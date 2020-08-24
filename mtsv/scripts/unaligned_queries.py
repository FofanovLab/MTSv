import logging
import numpy as np
import pandas as pd
import subprocess as sp
import click
from Bio import SeqIO

strip = str.strip
split = str.split
lstrip = str.lstrip


def get_last_query_id_in_fasta(query_file):
    tail = sp.check_output(['tail', '-n', '2', query_file])
    return int(tail.decode('utf-8').split("\n")[0].lstrip(">R").split("_")[0])


def get_unaligned_query_ids(max_query_id, hit_ids):
    """
    Returns an integer array of indices that are not in 
    hit_ids, assuming max_query_d is the max value.
    """
    queries = np.ones(max_query_id, dtype=np.bool)
    hit_ids = np.subtract(hit_ids, 1)
    queries[hit_ids] = False
    return np.where(queries)[0] + 1


def get_hit_ids(clp_file):
    ids = pd.read_csv(clp_file,
                      names=['query_id'],
                      converters={
                          'query_id': lambda x: split(lstrip(x, "R"), "_")[0]},
                      sep=':', usecols=[0])
    return ids['query_id'].to_numpy(dtype=np.int64)


def write_unaligned_queries(query_file, ids, outfile):
    """
    Writes query sequences in fasta format from
    fasta formatted query_file with
    ids in ids to outfile.
    query_file: (str) Path to query_file
    ids: (list like ints) ids that should be included
    outfile: (str) Path to file to write reads
    """
    sorted_ids = sorted(ids)
    count = 0
    with open(outfile, 'w') as out:
        with open(query_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if not sorted_ids:
                    break
                if int(lstrip(split(record.id, "_")[0], "R")) == sorted_ids[0]:
                    count += 1
                    sorted_ids.pop(0)
                    SeqIO.write(record, out, 'fasta')
    logging.info(
        "{0} unaligned queries have been written to {1}".format(
            count, outfile))


def get_unaligned_queries(query_file, clp_file, outfile):
    write_unaligned_queries(query_file,
                            get_unaligned_query_ids(
                                get_last_query_id_in_fasta(query_file),
                                get_hit_ids(clp_file)),
                            outfile)


@click.command()
@click.argument(
    'QUERY_FILE',
    type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument(
    'CLP_FILE',
    type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument(
    'OUTFILE', type=click.Path())
def main(query_file, clp_file, outfile):
    """
    Get unaligned queries (queries not present in CLP_FILE) from
    QUERY_FILE fasta. Write to OUTFILE.
    """
    get_unaligned_queries(query_file, clp_file, outfile)


if __name__ == "__main__":
    main()
