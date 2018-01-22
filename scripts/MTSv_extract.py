import argparse
import os.path as path
from Bio import SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()

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

def get_outfile(_path, file_name_list, ext):
    return path.join(
        _path, "{0}.{1}".format(
            "_".join(file_name_list), ext))

def taxid_lookup(species_name):
    return ncbi.get_name_translator(
        [species_name])[species_name][0]

def species_lookup(taxid):
    taxid = int(taxid)
    return ncbi.get_taxid_translator(
        [taxid])[taxid]

def mtsv_extract(
        project_name, taxid, species, all_data,
        sig_data, read_fasta, out_path):
    all_out = get_outfile(
        out_path, [project_name, taxid, "all"], "fasta")
    sig_out = get_outfile(
        out_path, [project_name, taxid, "signature"], "fasta")
    all_read_ids = parse_data(all_data, taxid)
    sig_read_ids = parse_data(sig_data, taxid)
    print("Writing files for taxid {0} ({1})".format(taxid, species))
    write_sequences(
        read_fasta, all_read_ids, sig_read_ids,
        all_out, sig_out)

def write_sequences(
        read_fasta, all_read_ids,
        sig_read_ids, all_out, sig_out):
    with open(all_out, 'w') as all_handle, \
    open(sig_out, 'w') as sig_handle:
        with open(read_fasta, "rU") as read_handle:
            for record in SeqIO.parse(read_handle, "fasta"):
                if record.id in all_read_ids:
                    SeqIO.write(record, all_handle, "fasta")
                if record.id in sig_read_ids:
                    SeqIO.write(record, sig_handle, "fasta")

def parse_data(file_path, taxid):
    strip = str.strip
    split = str.split
    set_add = set.add
    reads = set()
    with open(file_path, "r") as data:
        for line in data:
            # save time by checking 
            # if id is present before 
            # str manipulation
            if taxid not in line:
                continue
            read_id, taxa = split(strip(line), ":")
            if taxid in split(taxa, ","):
                set_add(reads, read_id)
    return reads

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
    prog='MTSv Extract',
    description="Extracts all sequences, including signature "
                "hits, that aligned to a given taxid or species name.",
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
        "reads", metavar="READS_FASTA", type=file_type,
        help="Path to FASTA file from MTSv-readprep"
    )

    PARSER.add_argument(
        "-o", "--out_path", type=path_type, default="./",
        help="Output directory"
    )


    LOOK_UP_GROUP = PARSER.add_mutually_exclusive_group(required=True)
    LOOK_UP_GROUP.add_argument(
        '-t', '--taxid', type=str, default=None,
        help="Extract sequences by taxid" 
    )

    LOOK_UP_GROUP.add_argument(
        '-s', '--species', type=str, default=None,
        help="Extract sequences by species name"
    )
 

    ARGS = PARSER.parse_args()

    if ARGS.taxid is None:
        ARGS.taxid = taxid_lookup(ARGS.species)
    else:
        ARGS.species = species_lookup(ARGS.taxid)

    mtsv_extract(
        ARGS.project_name,
        ARGS.taxid, ARGS.species,
        ARGS.all,
        ARGS.sig, ARGS.reads,
        ARGS.out_path)
