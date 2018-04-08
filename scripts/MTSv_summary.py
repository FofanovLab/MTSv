import argparse
from ete3 import NCBITaxa
import os.path as path
from  os import getcwd
from os.path import expanduser
from collections import defaultdict
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
import subprocess as sp 
import configparser
import sys
from Bio import SeqIO
from scipy.stats import binom_test

NCBI = NCBITaxa()

div_map = {2:'Bacteria', 10239: 'Viruses (excluding environmental sample',
           2157: 'Archaea', 12884: 'Viroids', 28384: "Other and synthetic sequences",
           2759: "Eukaryotes", 33090: "Green Plants", 4751: "Fungi",
           7742: "Vertebrates (excluding Primates, Chiroptera, Bos taurus, Canis lupus familiaris)",
           9443: "Primates (excluding Homo sapiens)",
           9397: "Chiroptera", 9913: "Bos Taurus", 9615: "Canis lupus familiaris",
           9606: "Homo sapiens"}


split = str.split
strip = str.strip
rsplit = str.rsplit


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
    lineage = NCBI.get_lineage(taxid)
    for level in lineage[::-1]:
        if level in div_map:
            return div_map[level]
    return "Unknown"

def parse_line(line):
    line = rsplit(strip(line), ":", 1)
    taxa = [int(strip(tax)) for tax in split(line[1], ",")]
    counts = [int(c) for c in split(line[0], "_")[1:]]
    return taxa, counts, line[0]


def parse_signature_hits(sig_file):
    data_dict = {}
    read_list = []
    with open(sig_file, 'r') as infile:
        for line in infile:
            taxa, counts, read = parse_line(line)
            read_list.append(read)
            for taxon in taxa:
                if taxon not in data_dict:
                    data_dict[taxon] = {
                        i: np.zeros(4, dtype=int)
                        for i in range(len(counts))}
                for sample, count in enumerate(counts):
                    data_dict[taxon][sample] += [count,
                        bool(count), count, bool(count)]
    return data_dict, read_list


def get_lineage_in_signature(taxon, signature_taxa):
    lineage = NCBI.get_lineage(taxon)
    for line in lineage[::-1][1:]:
        if line in signature_taxa:
            return line
    return False
        

def parse_row(counts, taxa):
    taxa = [int(strip(tax)) for tax in split(taxa, ",")]
    counts = [int(c) for c in split(counts, "_")[1:]]
    return counts, taxa


def parse_chunks_verbose(chunk, sig_reads):
    data_dict = {}
    chunk = chunk[~chunk[0].isin(sig_reads)]
    for _, (counts, taxa) in chunk.iterrows():
        counts, taxa = parse_row(counts, taxa) 
        for taxon in taxa:
            if taxon not in data_dict:
                data_dict[taxon] = {
                        i: np.zeros(4, dtype=int)
                        for i in range(len(counts))}
            for sample, count in enumerate(counts):
                data_dict[taxon][sample] += [count, bool(count), 0, 0]
    return data_dict
                

def parse_chunks(chunk, sig_reads, sig_taxa):
    data_dict = {}
    chunk = chunk[~chunk[0].isin(sig_reads)]
    for _, (counts, taxa) in chunk.iterrows():
        counts, taxa = parse_row(counts, taxa)    
        for taxon in np.intersect1d(taxa, sig_taxa, assume_unique=True):
                if taxon not in data_dict:
                    data_dict[taxon] = {
                        i: np.zeros(4, dtype=int)
                        for i in range(len(counts))}
                for sample, count in enumerate(counts):
                    data_dict[taxon][sample] += [count, bool(count), 0, 0]
    return data_dict     


def merge_dicts(sig_dict, dicts):
    print("Merging Dictionary")
    for dic in dicts:
        for k, v in dic.items():
            if k not in sig_dict:
                sig_dict[k] = v
            else:
                for kk, vv in v.items():
                    sig_dict[k][kk] += vv
    return sig_dict

def parse_all_hits(all_file, data_dict, sig_reads, threads=1, verbose=False):
    sig_taxa = list(data_dict.keys())
    n_rows = sum(1 for line in open(all_file))
    if verbose:
        func = partial(parse_chunks_verbose, sig_reads=sig_reads)
    else:
        func = partial(parse_chunks, sig_reads=sig_reads, sig_taxa=sig_taxa)
    p = Pool(threads)
    results = p.map(
        func, pd.read_csv(
            all_file, sep=":", header=None, chunksize=n_rows//threads))
    return merge_dicts(data_dict, results)


def get_seq_for_taxid(taxid, prune_config, outpath):
    fasta_path = "{0}_{1}.fasta".format(outpath, taxid)
    call = ["./MTSv_prune.py", "-c", "-cp",
            prune_config, "-txi", taxid, "-o",
            fasta_path]
    check_subprocess(call)
    return fasta_path


def get_sequence_from_taxa(prune_config, sig_taxa, threads):
    p = Pool(threads)
    p_get_seq_for_taxid = partial(get_seq_for_taxid,
                                  prune_config=prune_config)
    fastas = p.map(p_get_seq_for_taxid, sig_taxa)
    p.close()
    p.join()
    return fastas


def get_mask(seq_len, k_mer):
    mask = np.zeros(sum(seq_len))
    start = 0
    for l in seq_len:
        if l >= k_mer:
            mask[start + l-(k_mer-1): start + l] = 1
        else:
            mask[start: start + l] = 1
        start = start + l
    return mask


def get_random_kmers(fasta, kmer):
    seqs = list(SeqIO.parse(fasta, "fasta"))
    seq_lens = [len(s) for s in seqs]
    try:
        samples = np.random.choice(np.ma.compressed(np.ma.masked_array(
            np.arange(np.sum(seq_lens)), mask=get_mask(seq_lens, kmer))),
            100000, replace=False)
    except ValueError:
        samples = np.ma.compressed(np.ma.masked_array(
            np.arange(np.sum(seq_lens)), mask=get_mask(seq_lens, kmer)))
    seqs = "".join([str(s.seq) for s in seqs])
    kmers = Counter([seqs[position: position + kmer] for position in samples])
    return kmers


def get_sample_kmers(fastas, kmer):
    kmer_counters = []
    for fasta in fastas:
        kmer_counters.append(
            get_random_kmers(fasta, kmer))
    set_of_kmers=np.unique(
        np.array([list(k.keys()) for k in kmer_counters]).flatten())
    unique_kmers={k: np.zeros(len(kmer_counters), dtype=int)
                                for k in set_of_kmers}
    for taxa, kmer_counter in enumerate(kmer_counters):
        for kmer, count in kmer_counter.items():
            unique_kmers[kmer][taxa] = count
    return unique_kmers
            
def write_samples_to_file(outpath, unique_kmers):
    with open(outpath, 'w') as out:
        for i, (kmer, counts) in enumerate(unique_kmers.items()):
            out.write(">R{0}_{1}\n{2}\n".format(i, "_".join(counts), kmer))

def binner_helper(index_path, cmd):
    result_path = "{0}_binner.txt".format(index_path)
    this_cmd = cmd + ["--index", index_path,
        "--results", result_path]
    check_subprocess(this_cmd)
    return result_path


def run_binner(config, fasta_path, threads):
    cmd = [config["path"], "--fasta", fasta_path]
    for flag, value in config.items():
        if flag not in ["path", "index", "results"]:
            cmd += ["--{0}".format(flag), value]
    indices = np.genfromtxt(config['index'], dtype=str)
    p = Pool(threads)
    binner_p = partial(binner_helper, cmd=cmd)
    results = p.map(binner_p, indices)
    p.close()
    p.join()
    return results

def run_collapse(config, binner_results, collapse_results, threads):
    cmd = [config["path"]] + \
          [path for path in binner_results] + \
          ["--output", results_path, "--threads", threads]
    check_subprocess(cmd)

def run_signature(config, collapse_results, sig_results, threads):
    cmd = [config['path'], "--input",
           collapse_results, "--threads", threads,
           "--output", sig_results]
    for flag, value in config.items():
        if flag != "path":
            cmd += ["--{0}".format(flag), value]
    check_subprocess(cmd)

def check_subprocess(cmd):
    try:
        out= sp.check_output(cmd)
        print("{0} SUCCEEDED, result = {1}".format(" ".join(cmd), out))
    except sp.CalledProcessError as e:
        sys.exit("{0} FAILED, returned code {1}, output: {2}".format(
            " ".join(cmd), e.returncode, e.output))
    except OSError as e:
        sys.exit("{0} FAILED, returned code {1}, output: {2}".format(
            " ".join(cmd), e.returncode, e.output))
            

def get_expected_prop(config, sig_taxa, outpath, threads):
    #TODO: add a cutoff so only taxa with greater than X signature hits are analyzed
    config = configparser.ConfigParser()
    config.read(config)
    # Call prune with signature taxa
    fastas = get_sequence_from_taxa(config['prune']['config'], sig_taxa, threads)
    # Get random samples
    unique_kmers = get_sample_kmers(fastas, config['summary']['kmer'])
    # Format into fasta input file
    fasta_path = "{0}_summary.fasta".format(outpath)
    collapse_results = "{0}_collapse.txt".format(outpath)
    sig_results = "{0}_signature.txt".format(outpath)

    write_samples_to_file(fasta_path, unique_kmers)
    # Run binner with the same parameters
    binner_results = run_binner(config['binner'], fasta_path, threads)
    # Run collapse
    run_collapse(config['collapse'], binner_results, collapse_results)
    # Run signature
    run_signature(config['signature'],
        collapse_results, sig_results, threads)
    # Run summary, don't recursively call stats by stats=False
    summary = get_summary(collapse_results, sig_results,
        None, None, threads=threads, stats=False)
    # Get signature proportion for taxid
    prop_dict = {}
    for i, taxon in enumerate(sig_taxa):
        total = summary[taxon][i][1]
        sig = summary[taxon][i][3]
        prop_dict[taxon] = (sig, total)
    return prop_dict

def cohen_h(p1, p2):
    return 2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p2))

def get_summary(all_file, sig_file, outpath, config, threads=1, verbose=False, stats=True):
    print("Parsing Signature Hits")
    data_dict, sig_reads = parse_signature_hits(sig_file)
    print("Parsing All Hits")
    data_dict = parse_all_hits(all_file, data_dict, sig_reads, threads, verbose)
    taxid2name = NCBI.get_taxid_translator(data_dict.keys())
    if not stats:
        return data_dict
    else:
        print("Getting stats")
        exp_props = get_expected_prop(config, list(data_dict.keys()),
            path.join(outpath.rsplit("/", 1), "stats"), threads)
        data_list = []
        for taxa, samples in data_dict.items():
            row_list = [taxa, tax2div(taxa), taxid2name[taxa]]
            if taxa in exp_props:
                p = exp_props[taxa][0]/exp_props[taxa][1]
                row_list.append("{0}/{1}".format(
                    exp_props[taxa][0], exp_props[taxa][1]))
                row_list.append(str(p))
            else:
                p = 0
                row_list += [0, 0]
            for sample, value in samples.items():
                row_list += [value[0], value[1], value[2], value[3]]
                if value[3] != 0:
                    row_list += [binom_test(value[3], value[1], p), 
                                cohen_h(value[3]/value[1], p1)]
                else:
                    row_list += [0, 0]
            data_list.append(row_list)
        
        column_names = ["TaxID","Division", "Sci. Name", "Prop_frac", "Exp_Prop"]
        n_cols = len(data_list[0])
        n_samples = int((n_cols - 5)/5)
        for c in range(n_samples):
            column_names += ["Total Hits (S{})".format(c+1),
                                "Unique Hits (S{})".format(c+1),
                                "Signature Hits (S{})".format(c+1),
                                "Unique Signature Hits (S{})".format(c+1),
                                "P-value (S{})".format(c+1),
                                "Cohens_H (S{})".format(c+1)]

        data_frame = pd.DataFrame(
            data_list,
            columns=column_names)
        data_frame.to_csv(outpath+"_summary.csv")


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

    PARSER.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads for multiprocessing"
    )

    PARSER.add_argument(
        "-cf", "--config", type=file_type,
        help="Config file with parameters for binning, "
             "collapse, and signature"
    )

    ARGS = PARSER.parse_args()
    if ARGS.update:
        NCBI.update_taxonomy_database()

    if ARGS.taxdump is not None:
        NCBI.update_taxonomy_database(
            path.abspath(ARGS.taxdump))
    

    outpath=path.join(
        ARGS.out_path, ARGS.project_name)
    get_summary(ARGS.all, ARGS.sig, outpath, ARGS.config, ARGS.threads, ARGS.verbose)

