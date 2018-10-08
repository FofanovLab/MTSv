import argparse
import logging
import json
import pandas as pd
import numpy as np
from ete3 import NCBITaxa
from mtsv.utils import (
    config_logging, warn,
    get_precalculated_df, get_ete_ncbi)
from mtsv.parsing import file_type, outfile_type, positive_int

def get_sample_count(n_cols):
    n_cols -= 3
    return n_cols//4

def get_candidate_taxa(summary_file, signature_cutoff):
    df = pd.read_csv(summary_file, comment="#")
    n_samples = get_sample_count(df.shape[1])
    taxa_set = set()
    for sample in range(n_samples):
        d = df[df[
            "Unique Signature Hits (S{})".format(
                sample + 1)] > signature_cutoff]
        taxa_set.update(list(d.TaxID))
    ranks = NCBI.get_rank(list(taxa_set))
    if len(ranks) != len(taxa_set):
        LOGGER.info(
            "Some taxids were removed because they were not "
            "at the species level. Currently, mtsv analyze only "
            "supports species level analysis.")
    return [k for k, v in ranks.items() if v == "species"]

def check_supplied_taxids(can_taxa):
    ranks = NCBI.get_rank(can_taxa)
    if len(ranks) != len(can_taxa):
        LOGGER.info(
            "Some taxids have been removed either because they "
            "were not found in the database or because they were "
            "not at the species level. Currently, mtsv analyze only "
            "supports species level analysis.")
    return [k for k, v in ranks.items() if v == "species"]
    

def remove_prev_calculated(
    can_taxa, **params):
    try:
        df = get_precalculated_df()
    except FileNotFoundError:
        LOGGER.warn(
            "Expected data file not found, "
            "recalculating all expected values.")
        return can_taxa, '{}'
    # get taxids that match all params and return remaining
    # taxids that are not found
    matching_taxa = df.loc[
        (df["Database"] == params['db_path']) & 
        (df["Kmer_Size"] == params['kmer']) & 
        (df['Edits'] == params['edits']) & 
        (df['Seed_Size'] == params['seed_size']) & 
        (df['Seed_Gap'] == params['seed_gap']) & 
        (df['Min_Seeds'] == params['min_seeds']),
        ["Taxid", "Total_Hits", "Sig_Hits", "Ratio"]]
    return (np.setdiff1d(
        can_taxa, matching_taxa['Taxid']),
        
        json.dumps(matching_taxa.loc[(matching_taxa['Taxid'].isin(np.intersect1d(
            can_taxa, matching_taxa['Taxid']))), :].to_dict(orient="index"))
        )



def write_to_file(can_taxa, outfile):
    LOGGER.info("Writing to file {}".format(outfile))
    with open(outfile, 'w') as out:
        out.write("\n".join([str(t) for t in can_taxa]))

if __name__ == "__main__":
    try:
        config_logging(snakemake.log[0], "INFO")
        LOGGER = logging.getLogger(__name__)
        NCBI = get_ete_ncbi(snakemake.params['taxdump'])
        if snakemake.params['can_taxa_list'] != None:
            LOGGER.info(
                "Reading candidates from file {}".format(
                    snakemake.params['can_taxa_list']
                ))
            CANTAXA = check_supplied_taxids(np.array(np.genfromtxt(
                snakemake.params['can_taxa_list'], dtype=int), ndmin=1))
            LOGGER.info(
                "Read the following {0} taxa from file:\n{1}".format(
                    len(CANTAXA),
                    ",".join([str(c) for c in CANTAXA])
                ))
        else:
            LOGGER.info(
                "Retreiving candidate taxa from summary file {}".format(
                    snakemake.input[0]
                ))
            LOGGER.info(
                "Including only taxa that have greater than {} "
                "unique signature hits".format(
                    snakemake.params['cutoff']
                ))
            CANTAXA = get_candidate_taxa(
                snakemake.input[0],
                snakemake.params['cutoff'])
            LOGGER.info(
                "The following {0} taxa were retreived from file:\n{1}".format(
                    len(CANTAXA),
                    ",".join([str(c) for c in CANTAXA])
                ))
        PRECALC_CANTAXA = '{}'            
        if snakemake.params['use_data']:
            LOGGER.info("Removing taxa that already have estimates")
            CANTAXA, PRECALC_CANTAXA = remove_prev_calculated(
                CANTAXA, **snakemake.params['exp_db_params'])
            LOGGER.info(
                "The following taxa still need expected "
                "value estimates\n{}".format(
                    ",".join([str(c) for c in CANTAXA])
                ))
        write_to_file(CANTAXA, snakemake.output[1])
        with open(snakemake.output[0], 'w') as json_out:
            json_out.write(PRECALC_CANTAXA)
        if len(CANTAXA) > 300:
            warn(
                "The number of candidate taxa is very large "
                "which may result in a very large query fasta file "
                "that will take a long time to process and require a "
                "lot of memory. You may want to rerun with a stricter "
                "cutoff to reduce the size of the query fasta file "
                "or break up the candidate taxa into chunks and run them "
                "individually by passing each chunk into analyze using "
                "the --can_taxa_list option")
        LOGGER.info("Finished collecting candidate taxa")
            
            
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

        PARSER.add_argument(
            "--log", type=outfile_type, default="cantaxa.log",
            help="Path to log file"
        )

        ARGS = PARSER.parse_args()
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)
        NCBI = NCBITaxa(
                taxdump_file=ARGS.taxdump)

        write_to_file(get_candidate_taxa(
            ARGS.summary,
            ARGS.signature_cutoff),
            ARGS.outfile)



        




