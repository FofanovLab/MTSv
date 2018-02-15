import argparse
from ete3 import NCBITaxa
import os.path as path
from  os import getcwd
from os.path import expanduser
from collections import defaultdict
import pandas as pd
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
    read = split(line[0], "_")[0]
    return taxa, counts, read

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
        

def parse_all_hits(all_file, data_dict, sig_reads):
    signature_taxa = data_dict.keys()
    with open(all_file, 'r') as infile:
        for line in infile:
            taxa, counts, read = parse_line(line)
            if read in sig_reads:
                continue
            for taxon in taxa:
                if taxon not in data_dict:
                    # check if rolled up to parent level
                    # lineage_in_signature = get_lineage_in_signature(
                    #     taxon, signature_taxa)
                    # if lineage_in_signature:
                    #     taxon = lineage_in_signature
                    # else:
                    data_dict[taxon] = {
                        i: np.zeros(4, dtype=int)
                        for i in range(len(counts))}
                for sample, count in enumerate(counts):
                    data_dict[taxon][sample] += [count, bool(count), 0, 0]

    return data_dict 


def get_summary(all_file, sig_file, outpath, cutoff=None):
    print("Parsing Signature Hits")
    data_dict, sig_reads = parse_signature_hits(sig_file)
    print("Parsing All Hits")
    data_dict = parse_all_hits(all_file, data_dict, sig_reads)
    taxid2name = NCBI.get_taxid_translator(data_dict.keys())
    print("Writing to File")
    data_list = []
    for taxa, samples in data_dict.items():
        row_list = [taxa, tax2div(taxa), taxid2name[taxa]]
        for sample, value in samples.items():
            row_list += [value[0], value[1], value[2], value[3]]
        data_list.append(row_list)
    
    column_names = ["TaxID","Division", "Sci. Name"]
    for c in range(int((len(data_list[0]) - 3)/4)):
        column_names += ["Total Hits (S{})".format(c+1),
                            "Unique Hits (S{})".format(c+1),
                            "Signature Hits (S{})".format(c+1),
                            "Unique Signature Hits (S{})".format(c+1)]

    data_frame = pd.DataFrame(
        data_list,
        columns=column_names)

    data_frame.sort_values(
        by='Unique Signature Hits', ascending=False)
    if cutoff is not None:
        data_frame = data_frame[
            data_frame['Unique Hits'] > cutoff]

    data_frame.to_csv(outfile, index=False) 
        

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        prog="MTSv Summary",
        description="Summarize number of hits for each taxa, "
                    "including signature hits.",
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
        "--count_cutoff", type=int, default=None,
        help="Only report taxa with more than "
             "COUNT_CUTOFF unique read hits."
    )

    ARGS = PARSER.parse_args()
    if ARGS.update:
        NCBI.update_taxonomy_database()

    if ARGS.taxdump is not None:
        NCBI.update_taxonomy_database(
            path.abspath(ARGS.taxdump))
    
    outfile = path.join(
        ARGS.out_path, "{0}_summary.csv".format(ARGS.project_name))

    get_summary(ARGS.all, ARGS.sig, outfile, ARGS.count_cutoff)


