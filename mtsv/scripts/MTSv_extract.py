from Bio import SeqIO
from ete3 import NCBITaxa
import logging
import os
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from mtsv.utils import line_generator, config_logging 
from mtsv.parsing import parse_output_row
import gzip
import shutil


def get_decendants(taxids, des_flag):
    parents = []
    descendants = []
    for taxid in taxids:
        if des_flag:
            des = list(NCBI.get_descendant_taxa(
                taxid, collapse_subspecies=True))
            if taxid not in des:
                des.append(taxid)
        else:
            des = [taxid]
        descendants += des
        parents += [taxid for _ in range(len(des))]
    return (np.array(parents, dtype=int),
            np.array(descendants, dtype=int))


def get_targets(clps_data_chunk, parents, descendants):
    target_dict = defaultdict(list)
    for line in clps_data_chunk:
        data = parse_output_row(line)
        hits = np.where(np.in1d(descendants, data.taxa))[0]
        for hit in hits:
            target_dict[parents[hit]].append(data.read_name)
    return target_dict

def reduce_targets(taxids, targets):
    # ensure each taxid is in final dict even if there
    # were no hits
    all_targets = {tax:[] for tax in taxids}
    for target in targets:
        for k, v in target.items():
            all_targets[k].extend(v)
    return all_targets

def mtsv_extract(
    taxids, clps_file, query_fasta,
    descendants, by_sample, outpath, threads):
    # for higher level taxids, include all descendants in output
    parents, descendants = get_decendants(taxids, descendants)
    targets_partial = partial(
                        get_targets,
                        parents=parents,
                        descendants=descendants)
    get_lines = line_generator(clps_file, 5000)
    p = Pool(threads)
    target_dicts = p.imap(targets_partial, get_lines)
    p.close()
    p.join()
    targets = reduce_targets(taxids, target_dicts)
    # Don't build sequence dictionary if no
    # hits were found
    query_fasta_dict = {}
    if sum([len(t) for t in targets.values()]):
        query_fasta_dict = SeqIO.to_dict(
            SeqIO.parse(query_fasta, "fasta"))

    if by_sample:
        with open(clps_file, 'r') as infile:
            n_samples = len(
                parse_output_row(infile.readline()).counts)
            write_partial = partial(
                write_sequences_by_sample,
                query_fastas=query_fasta_dict,
                outpath=outpath,
                n_samples=n_samples)
    else:
        write_partial = partial(write_sequences,
                                query_fastas=query_fasta_dict,
                                outpath=outpath)
    p = Pool(threads)
    p.map(write_partial, targets.items())
    

def write_sequences_by_sample(
    targets, query_fastas, outpath, n_samples):                                                                                           
    if not len(targets[1]):
        logger.info(
            "There were no sequences for taxid: {}".format(
                targets[0]))
    sample_dict = {i:[] for i in range(n_samples)}
    for query in targets[1]:
        samples = np.where(
            np.array(query.split("_")[1:], dtype=int))[0]
        for sample in samples:
            sample_dict[sample].append(query)
    for sample, queries in sample_dict.items():
        fasta_file = os.path.join(
            outpath,
            "{0}_{1}.fasta".format(targets[0], sample + 1))
        fastq_file = os.path.splitext(fasta_file)[0] + ".fastq"
        with open(fastq_file, 'w') as qhandle, open(fasta_file, 'w') as ahandle:
            records = [query_fastas[query] for query in queries]
            SeqIO.write(records, ahandle, "fasta")
            for query, record in zip(queries, records):
                record.letter_annotations[
                    "phred_quality"] = [40] * len(record)
                # repeat for the number of times the query appeared
                # in the sample
                name = record.id
                for rep in range(int(query.split("_")[1:][sample])):
                    record.id = "{0}_{1}".format(name, rep)
                    SeqIO.write(record, qhandle, "fastq")
                    

def write_sequences(
    targets, query_fastas, outpath):
    fasta_file = os.path.join(
        outpath, str(targets[0]) + ".fasta")
    fastq_file = os.path.splitext(fasta_file)[0] + ".fastq"
    with open(fasta_file, 'w') as fhandle, open(fastq_file, 'w') as qhandle:
        if not len(targets[1]):
            logger.info(
                "There were no sequences for taxid: {}".format(
                    targets[0]))
        else:
            records = [query_fastas[target] for target in targets[1]]
            SeqIO.write(records, fhandle, "fasta")
            for query, record in zip(targets[1], records):
                record.letter_annotations[
                    "phred_quality"] = 40 * len(record)
                name = record.id
                for rep in range(
                    np.sum(np.array(query.split("_")[1:], dtype=int))):
                    record.id = "{0}_{1}".format(name, rep)
                    SeqIO.write(record, qhandle, 'fastq')


if __name__ == "__main__":
    config_logging(snakemake.log[0], "INFO")
    logger = logging.getLogger(__name__)    

    NCBI = NCBITaxa(taxdump_file=snakemake.params[0])

    mtsv_extract(
         snakemake.params[1],
         snakemake.input[1],
         snakemake.input[0],
         snakemake.params[3],
         snakemake.params[4],
         snakemake.params[2],
         snakemake.threads)
