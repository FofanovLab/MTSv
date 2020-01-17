import os
import hashlib
import logging
import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import namedtuple

from mtsv.scripts.summary import (
    get_chunked_reader, add_taxa_count,
    drop_rows_over_max_taxa, Lineage,
    get_signature_mask)
from mtsv.utils import (
    get_ete_ncbi, error)



BLOOM_SIZE = 500000


class BloomFilter:
    """
    VERY simple bloom filter with single
    hash function.
    """

    def __init__(self, array):
        self.array = array

    def get_index(self, item):
        return int(
            hashlib.sha256(
                bytes(item, 'utf-8')).hexdigest(), 16) % len(self.array)

    def add(self, item):
        """
        Add item to filter.
        """
        self.array[self.get_index(item)] = True

    def check(self, item):
        """
        Check if item in filter.
        """
        return self.array[self.get_index(item)]

    def get_bloom(self):
        """
        Return bloom array
        """
        return self.array


class BloomDataStore:
    """
    Bloom dataset object
    """

    def __init__(self, df, size=BLOOM_SIZE):
        """
        Initialize with a dataframe of bloom filters
        where each column is a bloom fitler for a 
        different id. New bloom filters will
        be the same size as existing or if none
        present, BLOOM_SIZE.
        """
        self.df = df
        self.size = df.shape[0] if df.shape[0] else size

    def add_new(self, name):
        """
        Add new id column initialized with all zeros.
        Raise KeyError if name is already in BloomDataStore.
        """
        if name in self.df:
            raise KeyError(
                "{0} already exists in BloomDataStore".format(name))
        self.df[name] = pd.Series(np.zeros(self.size), dtype=bool)

    def update_bloom(self, name, array):
        """
        Change array for id
        """
        self.df[name] = array

    def get_bloom(self, name):
        """
        Return array for id
        """
        return self.df[name]

    def remove_bloom(self, name):
        """
        Remove a bloom filter from dataset
        """
        self.df = self.df.drop(name, 1)

    def get_percent_filled(self, name):
        """
        Return how full a bloom filter is for id.
        """
        total = self.df[name].sum()
        return total/self.size

    def get_datastore(self):
        """
        Returns bloom dataset with any modifications
        """
        return self.df


class DataStore:
    def __init__(self, data_store_file):
        self.file = data_store_file
        self.open()

    def update(self, name, df, metadata):
        """
        Write dataframe and metadata to path name. 
        """
        self.store.put(name, df)
        self.store.get_storer(name).attrs.metadata = metadata

    def update_metadata(self, name, metadata):
        """
        Add metadata to existing dataset.
        Raises KeyError if not in datastore.
        """
        curr_meta = self.store.get_storer(name).attrs.metadata
        curr_meta.update(metadata)
        self.store.get_storer(name).attrs.metadata = curr_meta
    

    def get_metadata(self, name):
        """
        Retrieve metadata from dataset.
        Raises KeyError if not in datastore
        """
        return self.store.get_storer(name).attrs.metadata

    def open(self):
        """
        Opens datastore at path, creates path if
        it does not exist. 
        """
        # Avoid HDF5 file locking error
        os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
        if not os.path.isfile(self.file):
            os.makedirs(os.path.dirname(self.file), exist_ok=True)
        self.store = pd.HDFStore(self.file, complevel=5, complib='bzip2')


    def close(self):
        """
        Closes datastore
        """
        self.store.close()

    def get_dataset(self, name):
        """
        Returns dataset in datastore
        Raises KeyError if not in datastore
        """
        return self.store[name]

    def remove(self, name):
        """
        Removes dataset in datastore
        Raises KeyError if not in datastore
        """
        self.store.remove(name)

    def group_info(self):
        """
        Returns dictionary of information for
        groups in datastore including metadata,
        and approximate size.
        """
        groups = {}
        for group in self.store.keys():
            groups[group] = {
                'metadata': self.store.get_storer(group).attrs.metadata,
                'size': self.store[group].memory_usage().sum()}
        return groups


DataSetInfo = namedtuple('DataStoreInfo',
        ['artifact', 'main_ds_path',
        'exp_ds_path', 'bloom_ds_path'])

def get_dataset_info(artifact_path, main_ds_path):
    exp_ds_path = get_expected_data_path(main_ds_path)
    bloom_ds_path = get_bloom_data_path(main_ds_path)
    return DataSetInfo(
        artifact=artifact_path, main_ds_path=main_ds_path,
        exp_ds_path=exp_ds_path, bloom_ds_path=bloom_ds_path)


def get_hash_from_artifact(artifact_file):
    with open(artifact_file, 'rb') as hash_file:
        return hashlib.md5(hash_file.read()).hexdigest()



def init_expected_metadata(dataset_info):
    todays_date = str(datetime.datetime.now())
    return {"file_name": dataset_info.artifact,
            "creation_date": todays_date,
            "last_accessed": todays_date}

def init_expected_df():
    return pd.DataFrame([],
                      columns=["unique_signature", "unique"], dtype=np.float64)


def init_bloom_df():
    return pd.DataFrame([], dtype=np.bool)

def init_bloom_metadata():
    return {}


def update_datastore(
        clp_file_list, queries_list, taxid_list, datastore,
        dataset_info, max_taxa_per_query, taxdump):
    ds = DataStore(datastore)
    exp_df, exp_meta = get_expected_dataset(dataset_info, ds)
    exp_df = update_species_expected_dataframe(
    clp_file_list, taxid_list, exp_df, max_taxa_per_query)
    bloom_df, bloom_meta = get_bloom_dataset(
        dataset_info.bloom_ds_path, ds)
    exp_df, bloom_df, bloom_meta = update_genus_expected_dataframe(
        clp_file_list, taxid_list, queries_list,
        exp_df, bloom_df, bloom_meta, taxdump, max_taxa_per_query)
    exp_df.index = exp_df.index.astype(np.int64)
    ds.update(dataset_info.exp_ds_path, exp_df, exp_meta)
    ds.update(dataset_info.bloom_ds_path, bloom_df, bloom_meta)
    ds.close()


def update_genus_expected_dataframe(
        clp_file_list, taxid_list, queries_list, exp_df,
        bloom_df, bloom_meta, taxdump, max_taxa_per_query):
    bloom_store = BloomDataStore(bloom_df)
    ncbi = get_ete_ncbi(taxdump)
    lineage = Lineage(ncbi)
    genera = [
        lineage.get_lineage_levels(taxid)['genus'] for taxid in taxid_list]
    for genus, queries, clp_file in zip(genera, queries_list, clp_file_list):
        # Do not run species that have been rolled up to genus level
        # (aka species that don't have a genus)
        if ncbi.get_rank([genus])[int(genus)] != 'genus':
            continue
        # add new genus if not in metadata
        if genus not in bloom_meta:
            # add new empty zeros array column
            bloom_store.add_new(genus)
            bloom_meta.update({genus: "active"})
        # once a bloom filter is filled, the final estimate of expected
        # value is final. The genus is marked as inactive
        # and the bloom filter is removed.
        if bloom_meta[genus] == "inactive":
            continue
        genus_summary_data = get_genus_summary_data(
            clp_file, genus, max_taxa_per_query, queries,
            bloom_store, lineage)
        # add counts to expected dataframe 
        if genus in exp_df.index:
            exp_df.loc[genus] = np.add(exp_df.loc[genus], genus_summary_data)
        else:
            exp_df.loc[genus] = genus_summary_data
        # remove bloom filters that are nearly full and mark as inactive
        if bloom_store.get_percent_filled(genus) >= 0.9:
            bloom_store.remove_bloom(genus)
            bloom_meta.update({genus: "inactive"})
    return exp_df, bloom_df, bloom_meta


    

def map_query_sequence_to_ids(query_file, ids):
    """
    Ids is a list like of query ids
    Query_file is a fasta formated file that contains
    the same query ids. Pulls query sequences from
    queries file. Returns a dictionary mapping
    query id to sequence record.
    """
    # sort query id to make it easier to search for sequences
    # in sorted query_file
    sorted_ids = sort_query_ids_numerically(list(ids))
    seqs = {}
    with open(query_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if not sorted_ids:
                break
            if record.id == sorted_ids[0]:
                seqs[sorted_ids.pop(0)] = str(record.seq)
    return seqs




def sort_query_ids_numerically(ids):
    """
    Returns query ids in the format "R1_10" sorted
    ascending by the numerical part of the id
    (1 in this example).
    """
    return sorted(
        ids, key=lambda x: int(x.split("_")[0].lstrip("R")))


def run_bloom_filter(bloom, seqs):
    passed = []
    for name, seq in seqs.items():
        if not bloom.check(seq):
            bloom.add(seq)
            passed.append(name)
    return passed


            
def get_genus_summary_data(
        clp_file, genus,
        max_taxa_per_query, queries,
        bloom_store, lineage):
    # randomly sample 10% of remaining hits
    sampled_hits = get_raw_summary_dataframe(
        clp_file, max_taxa_per_query).sample(frac=0.1)
    # add taxa from sampled hits to lineage
    lineage.update_lineage(sampled_hits['taxa'])
    header, base_array = lineage.expand_taxa_array(
        sampled_hits['taxa'], "genus")
    expanded_df = pd.DataFrame(
        base_array, columns=header, index=sampled_hits.query_id)
    # get only queries that contain genus
    expanded_df = expanded_df[expanded_df[genus] == True]
    seqs = map_query_sequence_to_ids(queries, expanded_df.index)
    bloom = BloomFilter(bloom_store.get_bloom(genus))
    passed = run_bloom_filter(
        bloom, seqs)
    bloom_store.update_bloom(genus, bloom.get_bloom())
    expanded_df = expanded_df.loc[passed, :] # only rows that were not in filter
    unique_total = expanded_df[genus].sum()
    sig_mask = get_signature_mask(expanded_df)
    unique_sig = (expanded_df[genus] * sig_mask).sum()
    return np.array([unique_sig, unique_total], dtype=np.float64)


def get_bloom_dataset(bloom_path, ds):
    try:
        bloom_df = ds.get_dataset(bloom_path)
        bloom_meta = ds.get_metadata(bloom_path)
    except KeyError:
        bloom_df = init_bloom_df()
        bloom_meta = init_bloom_metadata()
    return bloom_df, bloom_meta


def get_expected_dataset(dataset_info, ds):
    try:
        exp_df = ds.get_dataset(dataset_info.exp_ds_path)
        exp_meta = ds.get_metadata(dataset_info.exp_ds_path)
    except KeyError:
        exp_df = init_expected_df()
        exp_meta = init_expected_metadata(dataset_info)
    return exp_df, exp_meta


def select_taxa_from_summary(df, taxa):
    """
    From dataframe with taxa column of tuples, filter
    rows that do not contain target taxa.
    """
    df['filter'] = df.taxa.apply(
        lambda x: True if int(taxa) in x else False)
    return df[df['filter'] == True]


def update_species_expected_dataframe(
        clp_file_list, taxid_list, exp_df,
        max_taxa_per_query):
    for clp_file, taxid in zip(clp_file_list, taxid_list):
        exp_df.loc[taxid] = get_species_summary_data(
            clp_file, taxid, max_taxa_per_query)
    return exp_df


def get_raw_summary_dataframe(clp_file, max_taxa_per_query):
    return drop_rows_over_max_taxa(
        add_taxa_count(
            get_chunked_reader(clp_file, chunksize=None)),
        max_taxa_per_query)


def get_species_summary_data(clp_file, taxid, max_taxa_per_query):
    df = select_taxa_from_summary(get_raw_summary_dataframe(
        clp_file, max_taxa_per_query), taxid)
    unique_sig = len(df[df.n_taxa == 1])
    unique_total = len(df)
    return np.array([unique_sig, unique_total], dtype=np.float64)


def format_group_info_into_table(datastore):
    info = datastore.group_info()
    data = {}
    header = ["Artifact File", "Creation Date", "Last Accessed", "Size (bytes)"]
    unique_groups = list(
        set(["/".join(p.split("/")[:-1]) for p in info.keys()]))
    for group in unique_groups:
        bloom = info.get("{}/bloom".format(group), {})
        expected = info.get("{}/expected".format(group), {})
        binning_params = [s.rsplit("_", 1)[-1] for s in group.split("/")[1:]]
        data_key = "-".join(binning_params)
        data[data_key] = [
            expected['metadata'].get('file_name', None),
            expected['metadata'].get('creation_date', None),
            expected['metadata'].get('last_accessed', None),
            expected.get('size', 0)]
        data[data_key][-1] += bloom.get('size', 0)
    df = pd.DataFrame.from_dict(
        data, orient="index", columns=header)
    df.index.name = "HASH-Kmer_size-Edits-Seed_Size-Seed_Gap-Min_Seeds"
    return df


def get_datastore_info_table(datastore):
    ds = DataStore(datastore)

    data_table = format_group_info_into_table(ds)
    pd.set_option('display.width',  200)
    pd.set_option('display.colheader_justify', 'center')
    ds.close()
    return data_table

def remove_datastore(datastore, ds_path_id):
    ds = DataStore(datastore)
    ds_path = "hash_{0}/kmer_size_{1}/edits_{2}/" \
        "seed_size_{3}/seed_gap_{4}/min_seeds_{5}".format(
            *ds_path_id.split("-"))
    try:
        ds.remove("{}/expected".format(ds_path))
        ds.remove("{}/bloom".format(ds_path))
    except KeyError:
        error("Dataset not found in datastore")
    finally:
        ds.close()

def get_expected_data_path(main_ds_path):
    return os.path.join(
        main_ds_path,
        "expected")


def get_bloom_data_path(main_ds_path):
    return os.path.join(
        main_ds_path,
        "bloom")


def get_main_datastore_path(artifact_path, binning_params):
    """
    Returns string indicating the path to
    dataset in datastore based on artifact hash and
    binning parameters used for run.
    """
    artifact_hash = get_hash_from_artifact(artifact_path)
    return "hash_{name}/kmer_size_{kmer_size}/edits_{edits}/"\
        "seed_size_{seed_size}/seed_gap_{seed_gap}"\
        "/min_seeds_{min_seeds}".format(
            name=artifact_hash, **binning_params)


def update_last_accessed_metadata(ds, dataset_path):
    metadata = {'last_accessed': str(datetime.datetime.now())}
    ds.update_metadata(dataset_path, metadata)


def candidates_not_in_database(
        candidates, datastore, main_ds_path):
    # get hash for artifacts file, this ensures that the
    # information in the file has not changed from previous
    # run
    ds = DataStore(datastore)
    # get data group path string
    ds_path = get_expected_data_path(main_ds_path)
    # get dataset if it exists
    try:
        taxa_in_ds = ds.get_dataset(ds_path).index
        # keep track of accessed dates so old files can be removed
        update_last_accessed_metadata(ds, ds_path)
        return get_candidates_not_in_database(taxa_in_ds, candidates)
    except KeyError:
        # No entry for this artifact or parameters so
        # all need to be added
        return list(map(str, candidates))
    finally:
        ds.close()


def get_candidates_not_in_database(taxa_in_ds, candidates):
    return list(map(str, list(np.setdiff1d(
        np.array(candidates, dtype=int),
        np.array(taxa_in_ds, dtype=int),
        assume_unique=True))))

def get_expected_data_for_analysis(datastore, main_ds_path):
    ds = DataStore(datastore)
    # get expected data path
    ds_path = get_expected_data_path(main_ds_path)
    df = ds.get_dataset(ds_path)
    ds.close()
    return df

    
