import logging
import argparse
import sys
import json
import os.path
import numpy as np
from itertools import chain
from multiprocessing import Pool
from Bio import SeqIO
from mtsv.utils import config_logging, error, line_generator#, fasta_generator
from mtsv.parsing import file_type, outpath_type, parse_query_id, outfile_type


def fasta_generator(file_name, n_records):
    with open(file_name, 'r') as handle:
        records = {}
        for record in SeqIO.parse(handle, 'fasta'):
            records[record.id] = record
            if len(records) == n_records:
                yield records
                records = {}
        if len(records):
            yield records
        return


class Result:
    def __init__(self, all_queries, unaligned_queries):
        self._all_queries = all_queries
        self._unaligned_queries = unaligned_queries
        self._all_query_counts = np.array(
            [parse_query_id(q.id) for q in all_queries], dtype=int)
        self._n_samples = len(self._all_queries[0])
        self._unaligned_query_counts = np.array(
            [parse_query_id(q.id) for q in unaligned_queries], dtype=int)
    
    @property
    def results(self):
        return {"total_unique_queries": self.total_unique_queries,
                "total_unique_unaligned_queries": 
                self.total_unique_unaligned_queries,
                "total_queries_by_sample": self.total_queries_by_sample,
                "total_queries": self.total_queries,
                "total_unaligned_queries": self.total_unaligned_queries,
                "total_unique_queries_by_sample": 
                self.total_unique_queries_by_sample,
                "total_unaligned_queries_by_sample": 
                self.total_unalianged_queries_by_sample,
                "total_unique_unaligned_queries_by_sample": 
                self.total_unique_unalianged_queries_by_sample,
                "unaligned_queries_by_sample": 
                self.unalianged_queries_by_sample}
    
    @property
    def unalianged_queries_by_sample(self):
        by_sample = [[] for sample in range(self._n_samples)]
        for query, counts in zip(
            self._unaligned_queries, self._unaligned_query_counts):
            for i, count in enumerate(counts):
                if count:
                    by_sample[i].append(query)
        return by_sample

    @property
    def total_unique_queries(self):
        return len(self._all_queries)

    @property
    def total_unique_unaligned_queries(self):
        return len(self._unaligned_queries)

    @property
    def total_queries(self):
        return np.sum(self._all_query_counts)

    @property
    def total_unaligned_queries(self):
        return np.sum(self._unaligned_query_counts)

    @property
    def total_queries_by_sample(self):
        return np.sum(self._all_query_counts, axis=0)

    @property
    def total_unique_queries_by_sample(self):
        return np.sum(np.array(self._all_query_counts, dtype=bool), axis=0)

    @property
    def total_unalianged_queries_by_sample(self):
        return np.sum(self._unaligned_query_counts, axis=0)

    @property
    def total_unique_unalianged_queries_by_sample(self):
        return np.sum(np.array(
            self._unaligned_query_counts, dtype=bool), axis=0)

# GLOBALS
hits = np.empty(0)


def parse_hits(line_list):
    return [l.split(":")[0] for l in line_list]

def read_hits(collapse_file, threads):
    global hits
    p = Pool(threads)
    file_gen = line_generator(collapse_file, 5000)
    hits = np.array(list(
        chain.from_iterable(p.imap(parse_hits, file_gen))), dtype=str)


def get_unaligned_queries_worker(fasta_queries):
    global hits
    result = Result(list(fasta_queries.values()),
        [fasta_queries[k] for k in np.setdiff1d(
            list(fasta_queries.keys()), hits)])
    return result.results


def write_sequences(unaligned_queries_by_sample, file_names):
    for queries, file_name in zip(unaligned_queries_by_sample, file_names):
        with open(file_name, 'a') as handle:
            SeqIO.write(queries, handle, "fasta")


def get_outfiles(sample_names, n_samples, outpath):
    file_names = []
    for sample in range(n_samples):
        try:
            name = sample_names[sample]
        except IndexError:
            LOGGER.warn(
                "Provided sample names are not equal to the number of samples.")
            name = "Sample_{}".format(sample + 1)
        except TypeError:
            name = "Sample_{}".format(sample + 1)
        file_names.append(
            os.path.join(
                outpath, "unaligned_queries_{}.fasta".format(name)))
    # clear outfiles
    [open(fn, 'w').close() for fn in file_names]
    return file_names
        

def peek_at_samples_in_fasta(fasta):
    with open(fasta, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return len(parse_query_id(record.id))


def add_summary_counts(total_dict, new_values):
    for key in new_values.keys():
        total_dict[key] += new_values[key]
    return total_dict

def init_total_dict(n_samples):
    return {"total_queries": 0,
            "total_unique_queries": 0,
            "total_hits": 0,
            "total_unaligned_queries": 0,
            "total_unique_unaligned_queries": 0,
            "total_unique_hits": 0,
            "total_queries_by_sample": np.zeros(n_samples, dtype=int),
            "total_unique_queries_by_sample": np.zeros(n_samples, dtype=int),
            "total_unaligned_queries_by_sample": np.zeros(n_samples, dtype=int),
            "total_unique_unaligned_queries_by_sample": np.zeros(
                n_samples, dtype=int),
            "total_hits_by_sample": np.zeros(n_samples, dtype=int),
            "total_unique_hits_by_sample": np.zeros(n_samples, dtype=int),
            }


def add_hit_counts(total_dict, hits, threads):
    p = Pool(threads)
    counts = p.map(parse_query_id, hits)
    p.close()
    p.join()
    counts = np.array(counts, dtype=int)
    total_dict['total_hits'] = np.sum(counts)
    total_dict['total_unique_hits'] = len(hits)
    total_dict['total_hits_by_sample'] = np.sum(counts, axis=0)
    total_dict['total_unique_hits_by_sample'] = np.sum(
        np.array(counts, dtype=bool), axis=0)
    return total_dict

def log_results(total_dict):
    for key, value in total_dict.items():
        LOGGER.info("{0}:\t{1}".format(
            key.replace("_", " "), value))

def write_summary_json(outpath, total_dict):
    summary_file = os.path.join(outpath, "query_stats.json")
    LOGGER.info("Query stats summary writing to {}".format(summary_file))
    # covert types for json encoding
    for key, value in total_dict.items():
        try:
            total_dict[key] = [int(v) for v in value]
        except TypeError:
            total_dict[key] = int(value)
    with open(summary_file, 'w') as json_out:
        json_out.write(json.dumps(total_dict))

def get_unaligned_queries(collapse_file, query_fasta,
        outpath, sample_names=None, threads=1):
    LOGGER.info("Reading query hits from: {}".format(collapse_file))
    read_hits(collapse_file, threads)

    n_samples = len(parse_query_id(hits[0]))
    LOGGER.info("{} sample(s) in file.".format(n_samples))

    if n_samples != peek_at_samples_in_fasta(query_fasta):
        msg = """
                The number of samples in query fasta does
                not match the number of sequences in merge file
              """
        LOGGER.error(msg)
        error(msg)
    outfiles = get_outfiles(sample_names, n_samples, outpath)
    LOGGER.info("Output will be written to: {}".format(", ".join(outfiles)))
    total_dict = init_total_dict(n_samples)
    total_dict = add_hit_counts(total_dict, hits, threads)

    with Pool(processes=threads) as pool:
        LOGGER.info("Reading queries from {}".format(query_fasta))
        fastas = fasta_generator(query_fasta, 100000)
        for chunk, result in enumerate(pool.imap(
                get_unaligned_queries_worker,
                fastas)):
            LOGGER.info("Writing Chunk {} to file".format(chunk))
            write_sequences(
                result['unaligned_queries_by_sample'], outfiles)
            del result['unaligned_queries_by_sample']
            add_summary_counts(total_dict, result)
    log_results(total_dict)
    write_summary_json(outpath, total_dict)
    LOGGER.info("FINISHED finding unaligned queries")

    




def parse_args():
    parser = argparse.ArgumentParser(
        prog="MTSv Unaligned Queries",
        description="Get queries that had no hits",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "clps", metavar="COLLAPSE_FILE", type=file_type,
        help="Path to merge output file."
    )

    parser.add_argument(
        "fasta", metavar="FASTA_FILE", type=file_type,
        help="Path to query fasta."
    )

    parser.add_argument(
        '-o', "--outpath", type=outpath_type,
        default="./",
        help="Path for unaligned queries output."
    )

    parser.add_argument(
        "--samples", nargs='+', default=None,
        help="Provide optional sample names for output files."
    )

    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads for multiprocessing"
    )

    parser.add_argument(
        "--log", type=outfile_type, default="./unaligned.log",
        help="Path of log file."
    )

    return parser.parse_args()
    
            

if __name__ == "__main__":

    try:
        config_logging(snakemake.log[0], "INFO")      
        LOGGER = logging.getLogger(__name__)
        get_unaligned_queries(snakemake.input[0],
            snakemake.input[1], snakemake.output[0],
            snakemake.params[0], snakemake.threads)
    except NameError:
        ARGS = parse_args()
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)
        get_unaligned_queries(
            ARGS.clps, ARGS.fasta, ARGS.outpath, ARGS.samples, ARGS.threads)




    

        




