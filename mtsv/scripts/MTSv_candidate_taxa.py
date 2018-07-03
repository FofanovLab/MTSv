import argparse
import pandas as pd
import numpy as np
from ete3 import NCBITaxa
from mtsv.parsing import file_type, outfile_type, positive_int

def get_sample_count(n_cols):
    n_cols -= 3
    return n_cols//4

def get_candidate_taxa(summary_file, outfile, signature_cutoff):
    df = pd.read_csv(summary_file)
    n_samples = get_sample_count(df.shape[1])
    taxa_set = set()
    for sample in range(n_samples):
        d = df[df[
            "Unique Signature Hits (S{})".format(
                sample + 1)] > signature_cutoff]
        taxa_set.update(list(d.TaxID))
    ranks = NCBI.get_rank(list(taxa_set))
    taxa = [k for k, v in ranks.items() if v == "species"]
    with open(outfile, 'w') as out:
        out.write("\n".join(taxa))

if __name__ == "__main__":
    try:
        NCBI = NCBITaxa(taxdump_file=snakemake.params[1])
        get_candidate_taxa(
            snakemake.input[0],
            snakemake.output[0],
            snakemake.params[0])
    except NameError:
        PARSER = argparse.ArgumentParser(
        prog="MTSv Candidate Taxa",
        description="Get list of candidate taxa from summary.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        PARSER.add_argument(
            "summary", metavar="SUMMARY_FILE", type=file_type,
            help="Path to summary output file."
        )

        PARSER.add_argument(
            "output", metavar="OUTPUT", type=outfile_type,
            help="Summary output"
        )

        PARSER.add_argument(
            "--signature_cutoff", type=positive_int, default=20,
            help="Minimum number of unique signature hits to "
                 "be considered a candidate taxa."
        )

        PARSER.add_argument(
            "--taxdump", type=file_type, default=None,
            help="Alternative path to taxdump. "
                "Default is home directory where ete3 "
                "automatically downloads the file."
        )

        ARGS = PARSER.parse_args()

        NCBI = NCBITaxa(
                taxdump_file=ARGS.taxdump)


        get_candidate_taxa(
            ARGS.summary,
            ARGS.output,
            ARGS.signature_cutoff)



        




