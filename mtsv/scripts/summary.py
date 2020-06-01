import os
import logging
from collections import defaultdict
import pandas as pd
import numpy as np
from itertools import chain
from ete3 import NCBITaxa
from Bio import SeqIO
import click
from multiprocessing import Pool
from functools import partial
from mtsv.utils import (get_ete_ncbi, config_logging)


LEVELS = ['species', 'genus', 'family',
                   'order', 'class', 'phylum', 'superkingdom']

class Lineage:
    def __init__(self, ncbi):
        self.ncbi = ncbi
        self.lineage = {}

    def get_lineage_levels(self, taxon):
        """
        For a given taxid, returns a lineage based on
        the levels in LEVELS. Returns a dictionary with
        LEVELS as keys.
        """
        lineage_map = {v: k for k, v in self.ncbi.get_rank(
            self.ncbi.get_lineage(taxon)).items()}
        lineage = {}
        prev_taxon = taxon
        # if a level is not present, fill with previous level
        # this avoids counting a signature hit if a taxon
        # drops out at a certain level.
        for level in LEVELS:
            prev_taxon = lineage_map.get(level, prev_taxon)
            lineage[level] = prev_taxon
        return lineage

    def update_lineage(self, new_taxa):
        """
        New_taxa is a series of lists of taxids.
        This method updates the lineage dict for any
        taxids that are not already included.
        """
        # get taxa that need to be added to the list
        unique_new_taxa = self._get_new_unique_taxa(new_taxa)
        new_lineage_dict = {
            taxid: self.get_lineage_levels(taxid)
            for taxid in unique_new_taxa}
        self.lineage.update(new_lineage_dict)

    def get_level_from_taxa(self, taxa, level):
        """
        Returns all unique taxids for a set of taxids
        at a given level.
        """
        return list(set([self.lineage[taxid][level] for taxid in taxa]))

    def _get_new_unique_taxa(self, new_taxa):
        """
        Given a new group of taxa, returns only the ones
        that are not already included in the lineage dict.
        """
        return list(np.setdiff1d(
            list(set(chain.from_iterable(new_taxa))),
            list(self.lineage.keys()), assume_unique=True))

    def expand_taxa_array(self, taxa, level):
        """
        Takes a series of lists of taxa, rolls up to the 
        provided level given self.lineage. Returns
        an expanded boolean array where each column is a unique
        taxid at the level and each row is a query with
        hits indicated by True. 
        """
        taxa_header = list(
            set(taxid[level] for taxid in self.lineage.values()))
        index_dict = {v: i for i, v in enumerate(taxa_header)}
        # get 2d boolean array of which taxids are hits
        # for each query
        base_array = np.zeros((len(taxa), len(taxa_header)))
        x_index = []
        y_index = []
        for i, taxa_set in enumerate(taxa):
            cols = [
                index_dict[t] for t in
                self.get_level_from_taxa(taxa_set, level)]
            x_index += [i] * len(cols)
            y_index += cols
        base_array[x_index, y_index] = 1
        return (taxa_header, base_array)


    def get_lineage_from_level(self, taxa, level):
        """
        Returns the lineage for all taxa in taxa from 
        the provied level on.
        """
        lineage_from_level = {}
        remaining_taxa = list(taxa.copy())
        for k, v in self.lineage.items():
            if v[level] in remaining_taxa:
                lineage_from_level[v[level]] = v
                remaining_taxa.remove(v[level])
            if len(remaining_taxa) == 0:
                break
        lineage_from_level = pd.DataFrame.from_dict(
            lineage_from_level, orient="index")
        return lineage_from_level[LEVELS[LEVELS.index(level) : ]]
        

def summary(
    infile, outfile, max_taxa_per_query,
    taxdump, chunksize=100000, threads=1):
    """
    Parses merged output from binning process into a csv with
    rows that have more than max_taxa_per_query hits removed.
    Returns a list of unique taxa that are present in the remaining rows.
    """
    logging.info(
        "Parsing merged binning file: {}".format(
            os.path.abspath(infile)))
    ncbi = get_ete_ncbi(taxdump)
    summary_table = parse_merged_file_in_chunks(
        taxdump, infile, chunksize, max_taxa_per_query, threads)
    summary_table = format_summary_table(summary_table, ncbi)
    summary_table.to_csv(outfile, index=False, float_format="%.0f")
    logging.info(
        "Finished writing summary to: {}".format(
            os.path.abspath(outfile)))



def format_summary_table(summary_table, ncbi):
    """
    Adds scientific names to summary table, reorders the columns
    and sorts the values. Returns a dataframe
    """
    # get scientific name from unique values
    summary_table = get_scientific_names(summary_table, ncbi)
    summary_table = validate_level(summary_table, ncbi)
    summary_table = summary_table[[
        'taxid', 'scientific_name', 'level', 'total', 'unique',
        'signature', 'unique_signature', 'weighted_support',
        'unique_weighted_support'] + LEVELS[::-1]]
    summary_table = summary_table.sort_values(
        by=["signature", "total"],
        na_position="first", ascending=False)    
    return summary_table



def process_chunk(chunk, taxdump, chunk_setting, max_taxa_per_query):
    chunk_index, chunk = chunk
    chunk_size = chunk.shape[0]
    chunk = process_merged_dataframe(chunk, max_taxa_per_query)
    ncbi = get_ete_ncbi(taxdump, quiet=True)
    lineage = Lineage(ncbi)
    lineage.update_lineage(chunk['taxa'])
    result = calculate_summary(chunk, lineage)
    logging.info(
        "Finished processing chunk: {0} to {1}".format(
            chunk_index * chunk_setting,
            chunk_index * chunk_setting + chunk_size
        ))
    return result
    

def parse_merged_file_in_chunks(taxdump, infile, chunksize,
    max_taxa_per_query, threads):
    """
    Main loop for parsing chunks of merged binning file.
    Returns a summary dataframe with counts and lineage information
    for each taxon.
    """

    summary_table = pd.DataFrame([])
    p = Pool(threads)
    process_chunk_func = partial(
        process_chunk, taxdump=taxdump,
        chunk_setting=chunksize,
        max_taxa_per_query=max_taxa_per_query)
    
    for result in p.imap(process_chunk_func,
            enumerate(get_chunked_reader(infile, chunksize))):
        summary_table = summary_table.append(result)
    logging.info("Aggregating summary data")
    # total up multiple entries for each taxid from chunks.
    summary_table = aggregate_rows(summary_table, 'taxid', AGG_FUNCT)
    summary_table['taxid'] = summary_table.index
    logging.info("Finished aggregating summary data")
    return summary_table

def agg_function_for_levels(x):
    return x.unique()[0]

AGG_FUNCT = {
    'unique': 'sum',
    'total': 'sum',
    'unique_signature': 'sum',
    'signature': 'sum',
    'weighted_support': 'sum',
    'unique_weighted_support': 'sum',
    'level': agg_function_for_levels,
'species': agg_function_for_levels,
    'genus': agg_function_for_levels,
    'family': agg_function_for_levels,
    'order': agg_function_for_levels,
    'class': agg_function_for_levels,
    'phylum': agg_function_for_levels,
    'superkingdom': agg_function_for_levels
}

def aggregate_rows(df, group, agg_func):
    return df.groupby(group).agg(agg_func)


def validate_level(summary_table, ncbi):
    """
    Some taxids in the database may not be at the species level initially.
    This will cause them to have an incorrect rank. Update to value 
    provided by ncbi database.
    """
    taxonomic_rank = ncbi.get_rank(summary_table['taxid'].unique())
    summary_table['level'] = summary_table['taxid'].apply(
        lambda x: taxonomic_rank.get(x, "Undefined"))
    return summary_table


def get_scientific_names(summary_table, ncbi):
    """
    Add 'scientific_name' column to summary_table dataframe
    from the 'taxid' column
    """
    scientific_name_map = ncbi.get_taxid_translator(
        summary_table['taxid'].unique())
    summary_table['scientific_name'] = summary_table['taxid'].apply(
        lambda x: scientific_name_map.get(x, "None"))
    return summary_table


def get_total_unique_counts(base_array):
    """
    Total unique counts are the number of query hits for each
    taxon.
    """
    return base_array.sum(axis=0)

def get_total_counts(base_array, counts):
    """
    Total counts are the number of query hits for each
    taxon times the number of times the query
    appeared in the data.
    """
    return (base_array.T * counts).T.sum(axis=0)

def get_signature_mask(base_array):
    """
    Finds query rows where there is only a single
    taxa hit (signature hit) and returns a mask that
    selects only these rows.
    """
    return np.equal(base_array.sum(axis=1), 1)

def get_total_signature_counts(base_array, sig_mask, counts):
    """
    Total signature counts are the number of signature
    query hits (where it was the only hit) for each
    taxon times the number of times the query appeared
    in the data.
    """
    return ((base_array.T * sig_mask) * counts).T.sum(axis=0)

def get_total_unique_signature_counts(base_array, sig_mask):
    """
    Total unique signature counts are the number of
    signature query hits (where it was the only hit)
    for each taxon.
    """
    return (base_array.T * sig_mask).T.sum(axis=0)

def get_weights(base_array):
    """
    Calculate weights for weighted support.
    """
    return np.divide(1, base_array.sum(axis=1))

def get_total_weighted_support(base_array, weights, counts):
    """
    Total weighted support is the number of query hits
    for each taxa, weighted by how many other taxa were
    also a hit for each query. This is multiplied by 
    how many times each query appeared in the data.
    """
    return ((base_array.T * weights) * counts).T.sum(axis=0)

def get_total_unique_weighted_support(base_array, weights):
    """
    Total unique weighted support is the number of query hits
    for each taxa, weighted by how many other taxa were
    also a hit for each query.
    """
    return (base_array.T * weights).T.sum(axis=0)


def calculate_summary(data, lineage):
    """
    Data is a data_frame with columns 'taxa' and 'query_counts' and
    lineage is a Lineage object.
    Function returns a dataframe of stats for each taxa
    in data at each rollup level in LEVELS. 
    """
    counts = data['query_counts']
    summary_table = pd.DataFrame([])
    for level in LEVELS:
        summary_table = summary_table.append(
            get_summary_table(lineage, data['taxa'], level, counts))
    return summary_table
        

def get_summary_table(lineage, taxa, level, counts):
    """
    Calculates count summaries for each taxa at a given level.
    returns a dataframe with columns: "unique",
    "total", "unique_signature",
    "signature", "unique_weighted_support",
    "weighted_support", "taxid", "level" and each taxonomic 
    level.
    """
    labels, base_array = lineage.expand_taxa_array(
        taxa, level)
    count_table = (
        calculate_counts(base_array, counts))
    count_table['taxid'] = labels
    count_table['level'] = level
    count_table = count_table.merge(
        get_lineage_for_summary_table(
            lineage, count_table['taxid'], level),
            left_on="taxid", right_index=True)
    return count_table


def get_lineage_for_summary_table(
    lineage, taxa, level):
    """
    For each taxonomic level, gets the lineage in a 
    from lineage and sets lower taxonmic levels to None.
    Returns a dataframe with columns for all
    taxonomic levels. 
    """
    lineage = lineage.get_lineage_from_level(taxa, level)
    missing = np.setdiff1d(LEVELS, lineage.columns)
    for miss in missing:
        lineage[miss] = None
    return lineage[LEVELS]


def calculate_counts(base_array, counts):
    """
    Using a boolean base_array, calculates summary counts
    for each taxon: "unique", "total", "unique_signature",
    "signature", "unique_weighted_support", and
    "weighted_support"
    """
    counts = counts.to_numpy(dtype=int)
    count_df = pd.DataFrame([])
    count_df['unique'] = get_total_unique_counts(base_array)
    count_df['total'] = get_total_counts(base_array, counts)
    sig_mask = get_signature_mask(base_array)
    count_df['unique_signature'] = get_total_unique_signature_counts(
        base_array, sig_mask)
    count_df['signature'] = get_total_signature_counts(
        base_array, sig_mask, counts)
    weights = get_weights(base_array)
    count_df['unique_weighted_support'] = get_total_unique_weighted_support(
        base_array, weights)
    count_df['weighted_support'] = get_total_weighted_support(
        base_array, weights, counts)
    return count_df
    




def get_chunked_reader(clp_infile, chunksize=100000):
    """
    Takes a clp_infile and returns a dataframe reader in chunksize
    chunks. The reader splits the input file on : and converts
    the comma separated list of taxa to the left of : into a tuple.
    The headers are "query_id" and "taxa". query_id is formated
    like "R100_10" where "R100" is a unique identifier and "10"
    is the number of times the query appeared in the data. Taxa
    is a column of taxid tuples of taxa hits for the query.
    """
    return pd.read_csv(clp_infile,
                       names=['query_id', 'taxa'], dtype={'query_id': str},
                       sep=':', chunksize=chunksize,
                       converters={
                           'taxa': lambda x: tuple(
                               [int(i) for i in x.split(',')])}
                       )


def process_merged_dataframe(chunk, max_taxa_per_query):
    """
    Transforms a dataframe with columns 'query_id', 'taxa'.
    A "taxa_count" column is added that provides the length
    of the values in "taxa". Rows that have greater than 
    max_taxa_per_query taxa are removed. The counts from
    query id are added to a new column "query_counts"
    """
    data = chunk.copy()
    return parse_query_counts(
        drop_rows_over_max_taxa(
            add_taxa_count(data),
            max_taxa_per_query))


def parse_query_counts(data):
    """
    Data is a dataframe with columns "query_id" and "taxa".
    Parses the count information from "query_id" into
    a separate "query_counts" column. Returns a dataframe
    with "query_id", "taxa" and "query_counts" columns.
    """
    data = data.copy()
    data["query_counts"] = data["query_id"].str.split(
        "_", expand=True)[1]
    return data


def drop_rows_over_max_taxa(
        data, max_taxa_per_query):
    """
    Drop rows where "n_taxa" is greater
    than max_taxa_per_query. This removes uninformative
    rows with too many taxa.
    """
    return data[data["n_taxa"] <= max_taxa_per_query]


def add_taxa_count(data, taxa_col="taxa", taxa_count_col="n_taxa"):
    """
    Data is a dataframe with column "taxa", the length of the
    values in taxa is calculated in a new column "n_taxa".
    A dataframe with added "n_taxa" column is returned.
    """
    data["n_taxa"] = data["taxa"].apply(
        lambda x: len(x))
    return data


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
                seqs[sorted_ids.pop(0)] = record
    return seqs


def sort_query_ids_numerically(ids):
    """
    Returns query ids in the format "R1_10" sorted
    ascending by the numerical part of the id
    (1 in this example).
    """
    return sorted(
        ids, key=lambda x: int(x.split("_")[0].lstrip("R")))

def parse_clp_file_in_chunks_for_extract(
        taxdump, clp_file, queries, taxid, outfile,
        chunksize, max_taxa_per_query):
    taxdump = get_ete_ncbi(taxdump)
    lineage = Lineage(taxdump)
    level = lineage.ncbi.get_rank([taxid])[int(taxid)]
    if level not in LEVELS:
        level = "species"
    with open(outfile, 'w') as out:
        for chunk in get_chunked_reader(clp_file, chunksize):
            chunk = process_merged_dataframe(chunk, max_taxa_per_query)
            lineage.update_lineage(chunk['taxa'])
            header, base_array = lineage.expand_taxa_array(
                chunk['taxa'], level
            )
            expanded_df = pd.DataFrame(
                base_array, columns=header, index=chunk.query_id)
            try:
                expanded_df = expanded_df[expanded_df[int(taxid)] == True]
            except KeyError:
                continue
            seqs = map_query_sequence_to_ids(queries, expanded_df.index)
            SeqIO.write(seqs.values(), out, 'fasta')
            


@click.command()
@click.option(
    '--max_taxa_per_query', default=100, type=click.IntRange(min=1, max=1000),
    help='Ignore queries that have more than max_taxa_per_query hits')
@click.option(
    '--chunksize', default=100000, type=click.IntRange(min=1),
    help="Set the number of data rows processed in each chunk."
)
@click.option(
    "--taxdump", default=None, type=click.Path(exists=True, dir_okay=False),
    help="""
    Provide path to taxdump file, otherwise default setting for 
    ete3 will be used.""")
@click.option(
    "--threads", "-t", default=1, type=click.IntRange(min=1, max=24),
    help="Set number of processes."
)
@click.argument('infile', type=click.Path(exists=True, dir_okay=False))
@click.argument('outfile', type=click.Path(writable=True))
def main(infile, outfile, taxdump, max_taxa_per_query, chunksize, threads):
    logging.basicConfig(
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=getattr(logging, "INFO"),
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')
    summary(infile, outfile, max_taxa_per_query, taxdump, chunksize, threads)


def main_snake(
    infile, outfile, max_taxa_per_query,
    taxdump, chunksize, log, threads):
    config_logging(log, "INFO")
    summary(
        infile, outfile, max_taxa_per_query,
        taxdump, chunksize, threads)


if __name__ == "__main__":
    try:
        main_snake(snakemake.input[0], snakemake.output[0],
        snakemake.params["max_taxa_per_query"],
        snakemake.params['taxdump'], snakemake.params['chunksize'],
        snakemake.log[0], snakemake.threads)
    except NameError:
        main()
