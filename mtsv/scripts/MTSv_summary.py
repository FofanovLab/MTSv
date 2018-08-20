import logging
import sys
import os
import re
import argparse
import pandas as pd
import numpy as np
from ete3 import NCBITaxa
from multiprocessing import Pool
from functools import partial
from mtsv.utils import config_logging, line_generator, get_ete_ncbi
from mtsv.parsing import parse_output_row, file_type, outfile_type


DIV_MAP = {2:"Bacteria", 10239: "Viruses (excluding environmental sample)",
           2157: "Archaea", 12884: "Viroids", 28384: "Other and synthetic sequences",
           2759: "Eukaryotes", 33090: "Green Plants", 4751: "Fungi",
           7742: "Vertebrates (excluding Primates and Chiroptera and Bos taurus and Canis lupus familiaris)",
           9443: "Primates (excluding Homo sapiens)",
           9397: "Chiroptera", 9913: "Bos Taurus", 9615: "Canis lupus familiaris",
           9606: "Homo sapiens"}


def tax2div(taxid):
    try:
        lineage = NCBI.get_lineage(taxid)
        for level in lineage[::-1]:
            if level in DIV_MAP:
                return DIV_MAP[level]
        return "Unknown"
    except ValueError:
        return "Unknown"


def parse_signature_hits(lines):
    data_dict = {}
    read_set = set()
    for line in lines:
        data = parse_output_row(line)
        read_set.add(data.read_id)
        for taxon in data.taxa:
            if taxon not in data_dict:
                data_dict[taxon] = {
                    i: np.zeros(4, dtype=int)
                    for i in range(len(data.counts))}
            for sample, count in enumerate(data.counts):
                data_dict[taxon][sample] += [count,
                    bool(count), count, bool(count)]
    return data_dict, read_set

def sig_reduce(iterator):
    dic, reads = zip(*iterator)
    read_set = set()
    for _read in reads:
        read_set.update(_read)
    data_dict = all_reduce(dic)
    return data_dict, read_set


def get_descendants(taxa):
    levels = NCBI.get_rank(taxa)
    descendants = {}
    for taxon in taxa:
        if taxon in levels and levels[taxon] != "species":
            try:
                for descendant in NCBI.get_descendant_taxa(
                        taxon, collapse_subspecies=True):
                    descendants[descendant] = taxon
            except AttributeError:
                # when rank does not exist
                for descendant in NCBI.get_descendant_taxa(
                        taxon, collapse_subspecies=False):
                    descendants[descendant] = taxon
    return descendants


def get_lineage_in_signature(taxon, signature_taxa):
    lineage = NCBI.get_lineage(taxon)
    for line in lineage[::-1][1:]:
        if line in signature_taxa:
            return line
    return False


def merge_dicts(sig_dict, all_dict):
    for k, v in all_dict.items():
        if k not in sig_dict:
            sig_dict[k] = v
        else:
            for kk, vv in v.items():
                sig_dict[k][kk] += vv
    return sig_dict


def parse_all_hits(lines, sig_taxa, sig_reads, descendants):
    all_data_dict = {}
    for line in lines:
        data = parse_output_row(line)
        if data.read_id in sig_reads:
            # already counted these in sig_hits
            continue
        # check if taxa are descendants of rolled up signature taxa
        # get those that overlap with signature
        sig_des = np.array([descendants[tax] for tax in np.intersect1d(
            data.taxa,
            list(descendants.keys()),
            assume_unique=True)], dtype=int)
        # get species taxa that overlap with sig taxa
        sig_spec = np.intersect1d(
            data.taxa, list(sig_taxa), assume_unique=True)
        # Unique taxa for read that are also signature
        # with roll up to signature taxa level
        taxa_iter = np.unique(np.append(sig_spec, sig_des))
        for taxon in taxa_iter:
            if taxon not in all_data_dict:
                all_data_dict[taxon] = {
                    i: np.zeros(4, dtype=int)
                    for i in range(len(data.counts))}
            for sample, count in enumerate(data.counts):
                all_data_dict[taxon][sample] += [count, bool(count), 0, 0]
    return all_data_dict

def all_reduce(iterator):
    data_dict = {}
    for dic in iterator:
        for taxon, samples in dic.items():
            if taxon in data_dict:
                for sample, counts in samples.items():
                    data_dict[taxon][sample] = np.add(
                    data_dict[taxon][sample], counts)
            else:
                data_dict[taxon] = samples
    return data_dict
        



def get_summary(all_file, sig_file, outfile, threads, header=None):
    LOGGER.info("Parsing Signature Hits")
    p = Pool(threads)
    get_lines = line_generator(sig_file, 5000)
    sig_results = p.imap(parse_signature_hits, get_lines)
    p.close()
    p.join()
    sig_dict, read_set = sig_reduce(sig_results)
    descendants = get_descendants(sig_dict.keys())
    LOGGER.info("Parsing All Hits")
    p = Pool(threads)
    get_lines = line_generator(all_file, 5000)
    parse_all_partial = partial(
        parse_all_hits,
        sig_taxa=list(sig_dict.keys()),
        sig_reads=read_set, descendants=descendants)
    all_results = p.imap(parse_all_partial, get_lines)
    p.close()
    p.join()
    all_dict = all_reduce(all_results)
    data_dict = merge_dicts(sig_dict, all_dict)
    taxid2name = NCBI.get_taxid_translator(data_dict.keys())
    LOGGER.info("Writing to File: {}".format(outfile))

    data_list = []
    for taxa, samples in data_dict.items():
        taxa_name = re.sub(
            '\W+', ' ', taxid2name[taxa]).strip() if taxa in taxid2name else "Undefined"
        row_list = [taxa, tax2div(taxa), taxa_name]
        for sample, value in samples.items():
            row_list += [value[0], value[1], value[2], value[3]]
        data_list.append(row_list)
    
    column_names = ["TaxID","Division", "Sci. Name"]
    n_cols = len(data_list[0])
    n_samples = int((n_cols - 3)/4)
    for c in range(n_samples):
        column_names += ["Total Hits (S{})".format(c+1),
                            "Unique Hits (S{})".format(c+1),
                            "Signature Hits (S{})".format(c+1),
                            "Unique Signature Hits (S{})".format(c+1)]

    data_frame = pd.DataFrame(
        data_list,
        columns=column_names)
    data_frame['total'] = data_frame[
        ["Unique Signature Hits (S{})".format(i)
        for i in range(1, n_samples + 1)]].sum(axis=1)
    data_frame = data_frame.sort_values('total', ascending=False)
    data_frame = data_frame.drop('total', axis=1)
    with open(outfile, 'w') as out:
        if header:
            out.write("#,,,{}\n".format(",,,,".join(header)))
    data_frame.to_csv(
            outfile, mode='a', index=False) 


if __name__ == "__main__":
    try:
        NCBI = get_ete_ncbi(snakemake.params[0])
        config_logging(snakemake.log[0], "INFO")      
        LOGGER = logging.getLogger(__name__)

        
        get_summary(
            snakemake.input[1],
            snakemake.input[0],
            snakemake.output[0],
            snakemake.threads,
            snakemake.params['header'])
    except NameError:
        PARSER = argparse.ArgumentParser(
        prog="MTSv Summary",
        description="Summarize number of total and signature "
                    "hits for each taxa.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        PARSER.add_argument(
            "all", metavar="COLLAPSE_FILE", type=file_type,
            help="Path to merge output file."
        )

        PARSER.add_argument(
            "sig", metavar="SIGNATURE_FILE", type=file_type,
            help="Path to signature output file."
        )
        
        PARSER.add_argument(
            "-o", "--output", type=outfile_type, default="./summary.csv",
            help="Summary output"
        )

        PARSER.add_argument(
            "--taxdump", type=file_type, default=None,
            help="Alternative path to taxdump. "
                "Default is home directory where ete3 "
                "automatically downloads the file."
        )

        PARSER.add_argument(
            "--threads", type=int, default=1,
            help="Number of threads for multiprocessing"
        )

        PARSER.add_argument(
            "--log", type=outfile_type, default="./summary.log",
            help="Path of log file."
        )

        ARGS = PARSER.parse_args()

        NCBI = NCBITaxa(
                taxdump_file=ARGS.taxdump)

        config_logging(ARGS.log, "INFO")      
        LOGGER = logging.getLogger(__name__)

        get_summary(
            ARGS.all,
            ARGS.sig,
            ARGS.output,
            ARGS.threads)
else:
    NCBI = NCBITaxa()



        




