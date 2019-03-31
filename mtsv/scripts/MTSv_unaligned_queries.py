import logging
import argparse
import sys
import json
import os.path
from os import popen
import numpy as np
from itertools import chain
from multiprocessing import Pool
from Bio import SeqIO
from mtsv.utils import config_logging, error, line_generator, fasta_generator
from mtsv.parsing import file_type, outpath_type, parse_query_id, outfile_type


strip = str.strip
split = str.split

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
    return [open(fn, 'w') for fn in file_names]

def close_outfiles(outfiles):
    [fn.close() for fn in outfiles]
        

def peek_at_samples_in_fasta(fasta):
    with open(fasta, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return len(parse_query_id(record.id))

def peek_at_samples_in_hits(hit_file):
    with open(hit_file, 'r') as handle:
        line = handle.readline()
        return len(line.split(":")[0].split("_")) - 1


def add_summary_counts(total_dict, new_values):
    for key in new_values.keys():
        total_dict[key] += new_values[key]
    return total_dict

def init_total_dict():
    return {"total_queries": 0,
            "total_unique_queries": 0,
            "total_hits": 0,
            "total_unaligned_queries": 0,
            "total_unique_unaligned_queries": 0,
            "total_unique_hits": 0,
            "total_queries_by_sample": 0,
            "total_unique_queries_by_sample": 0,
            "total_unaligned_queries_by_sample": 0,
            "total_unique_unaligned_queries_by_sample": 0,
            "total_hits_by_sample": 0,
            "total_unique_hits_by_sample": 0,
            }


def calculate_final_totals(total_dict):
    total_dict['total_unique_queries_by_sample'] = np.add(
        total_dict['total_unique_hits_by_sample'],
        total_dict['total_unique_unaligned_queries_by_sample'])
    total_dict['total_queries_by_sample'] = np.add(
        total_dict['total_hits_by_sample'],
        total_dict['total_unaligned_queries_by_sample'])
    total_dict['total_queries'] = np.sum(
        total_dict['total_queries_by_sample'])

def get_total_queries(query_fasta):
    n_lines = 5
    tail = ""
    while ">" not in tail:
        tail = popen("tail -n {0} {1}".format(n_lines, query_fasta)).read()
        n_lines *= 2
    return int(tail.split(">")[-1].split("\n")[0].split("_")[0][1:])


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




def process_line(line):
    vals = split(split(line, ":")[0], "_")
    return int(vals[0][1:]) - 1, np.array(vals[1:], dtype=int)

def get_unaligned_query_ids(hits, total_queries, total_dict):
    hits_idx = np.zeros(total_queries, dtype=bool)
    with open(hits, 'r') as infile:
        for line in infile:
            total_dict['total_unique_hits'] += 1
            idx, totals = process_line(line)
            hits_idx[idx] = True
            total_dict['total_hits_by_sample'] = np.add(
                total_dict['total_hits_by_sample'], totals)
            total_dict['total_unique_hits_by_sample'] = np.add(
                total_dict['total_unique_hits_by_sample'],
                np.array(totals, dtype=bool))
    total_dict['total_hits'] = np.sum(total_dict['total_hits_by_sample'])
    return np.where(hits_idx == False)[0]


def write_unaligned_query_seqs(queries, unaligned_query_ids, outfile):
    next_idx = 0
    current_record = 0
    adding_record = False
    with open(outfile, 'w') as outfile:
        with open(queries, 'r') as infile:
            for line in infile:
                is_header = (line[0] == ">")
                if is_header:
                    adding_record = False
                    try:
                        if current_record == unaligned_query_ids[next_idx]:
                            adding_record = True
                            next_idx += 1
                    except IndexError:
                        break
                    current_record += 1
                if adding_record:
                    outfile.write(line)

def write_unaligned_queries_by_sample(unaligned_fasta, outfiles, total_dict):
    with open(unaligned_fasta, 'r') as handle:
        for name, seq in read_fasta(handle):
            total_dict['total_unique_unaligned_queries'] += 1
            counts = [int(c) for c in split(name, ("_"))[1:]]
            total_dict['total_unaligned_queries'] += sum(counts)
            total_dict['total_unaligned_queries_by_sample'] = np.add(
                total_dict['total_unaligned_queries_by_sample'],
                counts)
            total_dict['total_unique_unaligned_queries_by_sample'] = np.add(
                total_dict['total_unique_unaligned_queries_by_sample'],
                np.array(counts, dtype=bool))
            for outfile, count in zip(outfiles, counts):
                write_by_sample(outfile, count, name, seq)
    close_outfiles(outfiles)

            
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
            
       
def write_by_sample(file_handle, count, name, seq):
    if count == 0:
        return
    file_handle.write(
        ">{0}_{1}\n{2}\n".format(split(name, "_")[0], count, seq))



def get_unaligned_queries(collapse_file, query_fasta,
        outpath, sample_names=None):
    LOGGER.info("Reading query hits from: {}".format(collapse_file))
    # read_hits(collapse_file, threads)
    n_samples = peek_at_samples_in_hits(collapse_file)
    # n_samples = len(parse_query_id(hits[0]))
    LOGGER.info("{} sample(s) in file.".format(n_samples))

    if n_samples != peek_at_samples_in_fasta(query_fasta):
        msg = """
                The number of samples in query fasta does
                not match the number of sequences in merge file
              """
        LOGGER.error(msg)
        error(msg)

    total_dict = init_total_dict()
    # total_dict = add_hit_counts(total_dict, hits, threads)
    total_queries = get_total_queries(query_fasta)
    total_dict['total_unique_queries'] = total_queries

    unaligned_query_ids = get_unaligned_query_ids(
        collapse_file, total_queries, total_dict)


    unaligned_queries_file = os.path.join(
        outpath, "unaligned_queries.fasta")
    LOGGER.info("Finding unaligned queries, writing to {}".format(
        unaligned_queries_file
    ))
    write_unaligned_query_seqs(
        query_fasta, unaligned_query_ids,
        unaligned_queries_file)
    LOGGER.info("Finished finding unaligned queries")

    LOGGER.info("Sorting unaligned queries by sample")
    outfiles = get_outfiles(sample_names, n_samples, outpath)
    LOGGER.info("Sorting unaligned queries by sample, writing to: {}".format(
        ", ".join([f.name for f in outfiles])))

    write_unaligned_queries_by_sample(
        unaligned_queries_file, outfiles, total_dict)
    LOGGER.info("Finished sorting unaligned queries by sample")
    calculate_final_totals(total_dict)
    log_results(total_dict)
    LOGGER.info("Writing summary json")
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
        "--log", type=outfile_type, default="./unaligned.log",
        help="Path of log file."
    )

    return parser.parse_args()
    
            

if __name__ == "__main__":

    try:
        config_logging(snakemake.log[0], "INFO")      
        LOGGER = logging.getLogger(__name__)
        get_unaligned_queries(snakemake.input[0],
            snakemake.input[1], snakemake.params[1],
            snakemake.params[0])
    except NameError:
        ARGS = parse_args()
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)
        get_unaligned_queries(
            ARGS.clps, ARGS.fasta, ARGS.outpath, ARGS.samples)




    

        




