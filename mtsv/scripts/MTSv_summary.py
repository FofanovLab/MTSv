import argparse
from ete3 import NCBITaxa
import os.path as path
from  os import getcwd
from os.path import expanduser
from collections import defaultdict, namedtuple
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial



div_map = {2:"Bacteria", 10239: "Viruses (excluding environmental sample)",
           2157: "Archaea", 12884: "Viroids", 28384: "Other and synthetic sequences",
           2759: "Eukaryotes", 33090: "Green Plants", 4751: "Fungi",
           7742: "Vertebrates (excluding Primates, Chiroptera, Bos taurus, Canis lupus familiaris)",
           9443: "Primates (excluding Homo sapiens)",
           9397: "Chiroptera", 9913: "Bos Taurus", 9615: "Canis lupus familiaris",
           9606: "Homo sapiens"}


split = str.split
strip = str.strip
rsplit = str.rsplit

Record = namedtuple('Record', ['read_id', 'counts', 'taxa'])

def path_type(input_path):
    if not path.isdir(input_path):
        raise argparse.ArgumentTypeError(
            "Not a valid path: {}".format(input_path))
    return path.abspath(input_path)

def file_type(input_file):
    if not path.isfile(input_file):
        raise argparse.ArgumentTypeError(
            "Not a valid file path: {}".format(input_file))
    return path.abspath(input_file)

def tax2div(taxid):
    try:
        lineage = NCBI.get_lineage(taxid)
        for level in lineage[::-1]:
            if level in div_map:
                return div_map[level]
        return "Unknown"
    except ValueError:
        return "Unknown"

def parse_line(line):
    line = rsplit(strip(line), ":", 1)
    taxa = [int(strip(tax)) for tax in split(line[1], ",")]
    counts = [int(c) for c in split(line[0], "_")[1:]]
    return taxa, counts, line[0]


def parse_signature_hits(sig_file):
    data_dict = {}
    read_set = set()
    descendant_lookup = {}
    with open(sig_file, 'r') as infile:
        for line in infile:
            data = parse_row(line)
            read_set.add(data.read_id)
            levels = NCBI.get_rank(data.taxa)
            for taxon in data.taxa:
                if taxon not in data_dict:
                    data_dict[taxon] = {
                        i: np.zeros(4, dtype=int)
                        for i in range(len(data.counts))}
                    if taxon in levels and levels[taxon] != 'species':
                        for descendant in NCBI.get_descendant_taxa(
                            taxon, collapse_subspecies=True):
                            descendant_lookup[descendant] = taxon
                for sample, count in enumerate(data.counts):
                    data_dict[taxon][sample] += [count,
                        bool(count), count, bool(count)]
    return data_dict, read_set, descendant_lookup


def get_lineage_in_signature(taxon, signature_taxa):
    lineage = NCBI.get_lineage(taxon)
    for line in lineage[::-1][1:]:
        if line in signature_taxa:
            return line
    return False
        

def parse_row(row):
    read_id, taxa = split(row, ":")
    taxa = np.array([tax for tax in split(taxa, ",")], dtype=int)
    counts = np.array([c for c in split(read_id, "_")[1:]], dtype=int)
    read_id = split(read_id, "_")[0]
    return Record(read_id=read_id, counts=counts, taxa=taxa)

def merge_dicts(sig_dict, all_dict):
    print("Merging Dictionary")
    for k, v in all_dict.items():
        if k not in sig_dict:
            sig_dict[k] = v
        else:
            for kk, vv in v.items():
                sig_dict[k][kk] += vv
    return sig_dict


def parse_all_hits(all_file, data_dict, sig_reads, descendants, verbose=False):
    sig_taxa = set(data_dict.keys())
    all_data_dict = {}
    with open(all_file, 'r') as infile:
        for count, line in enumerate(infile):
            print("LINE: {}".format(count))
            data = parse_row(line)
            if data.read_id in sig_reads:
                    # already counted these
                    continue
            sig_spec = data.taxa
            sig_des = np.intersect1d(
                    data.taxa,
                    list(descendants.keys()),
                    assume_unique=True)
            sig_des = np.array(
                    [descendants[tax] for tax in sig_des], dtype=int)
            if not verbose:
                sig_spec = np.intersect1d(
                    data.taxa, list(sig_taxa), assume_unique=True)
            taxa_iter = np.unique(np.append(sig_spec, sig_des))
            
            for taxon in taxa_iter:
                if taxon not in all_data_dict:
                    all_data_dict[taxon] = {
                            i: np.zeros(4, dtype=int)
                            for i in range(len(data.counts))}
                for sample, count in enumerate(data.counts):
                    all_data_dict[taxon][sample] += [count, bool(count), 0, 0]

                
    return merge_dicts(data_dict, all_data_dict)


def get_summary(all_file, sig_file, outpath, verbose=False):
    print("Parsing Signature Hits")
    data_dict, sig_reads, descendants = parse_signature_hits(sig_file)
    print("Parsing All Hits")
    data_dict = parse_all_hits(
        all_file, data_dict, sig_reads, descendants, verbose)
    taxid2name = NCBI.get_taxid_translator(data_dict.keys())
    print("Writing to File")

    data_list = []
    for taxa, samples in data_dict.items():
        taxa_name = taxid2name[taxa] if taxa in taxid2name else "Undefined"
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
    data_frame.to_csv(outfile, index=False) 



if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        prog="MTSv Summary",
        description="Summarize number of total and signature "
                    "hits for each taxa.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    PARSER.add_argument(
        "project_name", metavar='PROJECT_NAME', type=str,
        help="Project name and output file prefix"
    )

    PARSER.add_argument(
        "all", metavar="COLLAPSE_FILE", type=file_type,
        help="Path to MTSv-collapse output file"
    )

    PARSER.add_argument(
        "sig", metavar="SIGNATURE_FILE", type=file_type,
        help="Path to MTSv-inform output file"
    )
    
    PARSER.add_argument(
        "-o", "--out_path", type=path_type, default="./",
        help="Output directory"
    )

    PARSER.add_argument(
        "--update", action="store_true",
        help="Update taxdump"
    )

    PARSER.add_argument(
        "--taxdump", type=file_type, default=None,
        help="Alternative path to taxdump. "
             "Default is home directory where ete3 "
             "automatically downloads the file."
    )

    PARSER.add_argument(
        "--verbose", action="store_true",
        help="Report all taxa, not just signature hits"
    )

    ARGS = PARSER.parse_args()
    if ARGS.update:
        NCBI = NCBITaxa()
        NCBI.update_taxonomy_database()

    elif ARGS.taxdump is not None:
        NCBI = NCBITaxa(
            taxdump_file=path.abspath(ARGS.taxdump))
    
    else:
        NCBI = NCBITaxa()

    outfile = path.join(
        ARGS.out_path, "{0}_summary.csv".format(ARGS.project_name))

    get_summary(ARGS.all, ARGS.sig, outfile, ARGS.verbose)

