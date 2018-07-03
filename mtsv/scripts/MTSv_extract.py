from Bio import SeqIO
from ete3 import NCBITaxa
import logging
import os
import argparse
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from mtsv.utils import config_logging 
from mtsv.parsing import parse_output_row, file_type, outfile_type, outpath_type
import gzip
import shutil



def line_generator(file_name, n_lines):
    go = True
    with open(file_name, 'r') as infile:
        while go:
            lines = []
            for _ in range(n_lines):
                l = infile.readline()
                if l == "":
                    go = False
                    break
                lines.append(l)
            yield lines
        return 


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
    taxids = [int(tax) for tax in taxids]
    # for higher level taxids, include all descendants in output
    parents, descendants = get_decendants(taxids, descendants)
    targets_partial = partial(
                        get_targets,
                        parents=parents,
                        descendants=descendants)
    get_lines = line_generator(clps_file, 5000)
    LOGGER.info("Extracting queries from {}".format(taxids))
    p = Pool(threads)
    target_dicts = p.imap(targets_partial, get_lines)
    p.close()
    p.join()
    targets = reduce_targets(taxids, target_dicts)
    LOGGER.info("Done extracting queries")

    if by_sample:
        with open(clps_file, 'r') as infile:
            n_samples = len(
                parse_output_row(infile.readline()).counts)
            write_partial = partial(
                write_sequences_by_sample,
                query_fastas=query_fasta,
                outpath=outpath,
                n_samples=n_samples)
    else:
        write_partial = partial(write_sequences,
                                query_fastas=query_fasta,
                                outpath=outpath)
    p = Pool(threads)
    p.map(write_partial, targets.items())
    LOGGER.info("Done writing to file")

    p.close()
    p.join()


def write_sequences_by_sample(
    targets, query_fastas, outpath, n_samples):
    LOGGER.info("Writing to file by sample")
                                                                                       
    if not len(targets[1]):
        LOGGER.info(
            "There were no sequences for taxid: {}".format(
                targets[0]))
    sample_dict = {i:[] for i in range(n_samples)}
    for query in targets[1]:
        samples = np.where(
            np.array(query.split("_")[1:], dtype=int))[0]
        for sample in samples:
            sample_dict[sample].append(query)
    for sample, queries in sample_dict.items():
        queries = set(queries)
        fasta_file = os.path.join(
            outpath,
            "{0}_{1}.fasta".format(targets[0], sample + 1))
        fastq_file = os.path.splitext(fasta_file)[0] + ".fastq"
        with open(fastq_file, 'w') as qhandle, open(fasta_file, 'w') as ahandle:
            with open(query_fastas, 'r') as qfasta:
                for record in SeqIO.parse(qfasta, 'fasta'):
                    if record.id in queries:
                        SeqIO.write(record, ahandle, "fasta")
                        record.letter_annotations[
                            "phred_quality"] = [40] * len(record)
                        name = record.id
                        # repeat for the number of times the query appeared
                        # in the sample
                        for rep in range(int(name.split("_")[1:][sample])):
                            record.id = "{0}_{1}".format(name, rep)
                            SeqIO.write(record, qhandle, "fastq")
                    

def write_sequences(
    targets, query_fastas, outpath):
    LOGGER.info("Writing to file")
    fasta_file = os.path.join(
        outpath, str(targets[0]) + ".fasta")
    fastq_file = os.path.splitext(fasta_file)[0] + ".fastq"
    tar = set(targets[1])
    with open(fasta_file, 'w') as fhandle, open(fastq_file, 'w') as qhandle:
        if not len(tar):
            LOGGER.info(
                "There were no sequences for taxid: {}".format(
                    targets[0]))
        else:
            with open(query_fastas, 'r') as qfasta:
                for record in SeqIO.parse(qfasta, 'fasta'):
                    if record.id in tar:
                        SeqIO.write(record, fhandle, "fasta")
                        record.letter_annotations[
                                "phred_quality"] = [40] * len(record)
                        name = record.id
                        for rep in range(
                            np.sum(np.array(name.split("_")[1:], dtype=int))):
                            record.id = "{0}_{1}".format(name, rep)
                            SeqIO.write(record, qhandle, 'fastq')


if __name__ == "__main__":
    try:
        config_logging(snakemake.log[0], "INFO")
        LOGGER = logging.getLogger(__name__)    

        NCBI = NCBITaxa(taxdump_file=snakemake.params[0])

        mtsv_extract(
            snakemake.params[1],
            snakemake.input[1],
            snakemake.input[0],
            snakemake.params[3],
            snakemake.params[4],
            snakemake.params[2],
            snakemake.threads)
    except NameError:
        
        PARSER = argparse.ArgumentParser(
            prog="MTSv Extract Queries for Taxa",
            description="Pull out reads for taxa.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )     
        PARSER.add_argument(
            "clp", metavar="CLP_FILE", type=file_type,
            help="Path to collapse file."
        )

        PARSER.add_argument(
            "query_fasta", metavar="QUERY_FASTA", type=file_type,
            help="Path to query fasta file."
        )

        PARSER.add_argument(
            "outpath", metavar="OUTPATH", type=outpath_type,
            help="Output directory"
        )

        PARSER.add_argument(
            "--taxdump", type=file_type, default=None,
            help="Path to taxdump file."
        )

        PARSER.add_argument(
            "--taxids", nargs="+", type=int,
            help="Taxids to pull."
        )

        PARSER.add_argument(
            "--descendants", action='store_true',
            help="Include all descendants in search"
        )

        PARSER.add_argument(
            "--by_sample", action='store_true',
            help="Break up sequences by sample."
        )

        PARSER.add_argument(
            "--log", type=outfile_type, default="extract.log",
            help="Name of log file."
        )

        PARSER.add_argument(
            "--threads", type=int, default=1,
            help="Number of threads to spawn."
        )

        ARGS = PARSER.parse_args()

        NCBI = NCBITaxa(taxdump_file=ARGS.taxdump)  
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)

        mtsv_extract(
            ARGS.taxids, ARGS.clp, ARGS.query_fasta,
            ARGS.descendants, ARGS.by_sample, ARGS.outpath, ARGS.threads)


