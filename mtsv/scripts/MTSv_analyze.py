import argparse
import numpy as np
import pandas as pd
from collections import namedtuple
from functools import partial
from scipy.stats import binom_test
from mtsv.parsing import file_type, outfile_type


EXP = namedtuple('exp', ["total", "sig", "ratio"])

def get_sample_count(n_cols):
    n_cols -= 3
    return n_cols//4


def get_sorted_headers(n_cols):
    headers = ["TaxID", "Division", "Sci. Name", 'exp_prop']
    sample_headers = ["Total Hits", "Unique Hits", "Signature Hits",
                      "Unique Signature Hits", "obs_prop", "p_value", "cohen_h"]
    for idx in range(1, n_cols + 1):
        for sample_header in sample_headers:
            headers.append(format_header(sample_header, idx))
    return headers


def get_cols(idx):
    return ["Unique Hits (S{})".format(idx),
            "Unique Signature Hits (S{})".format(idx)]


def process_stats(obs_sum, db, adj=100):
    n_samples = get_sample_count(obs_sum.shape[1])
    exp_prop = partial(get_exp_prop, database=db)
    obs_sum['exp_prop'] = obs_sum.TaxID.apply(exp_prop)
    for sample in range(1, n_samples + 1):
        fmt = partial(format_header, sample=sample)
        obs_sum[fmt('obs_prop')] = np.divide(
            obs_sum[fmt("Unique Signature Hits")],
            obs_sum[fmt("Unique Hits")])
        adj_num = obs_sum[fmt('obs_prop')] * adj
        obs_sum[fmt('p_value')] = [
            get_binom([adj_n, adj, exp]) for
            adj_n, exp in zip(adj_num, obs_sum['exp_prop'])]
        obs_sum[fmt('cohen_h')] = [
            cohen_h([obs, exp])
            for (obs, exp) in zip(
                obs_sum[fmt("obs_prop")],
                obs_sum["exp_prop"])]
    return obs_sum[get_sorted_headers(n_samples)]
    

def cohen_h(row):
    p1, p2 = row
    try:
        return abs(
            2 * np.arcsin(np.sqrt(p1)) -
            2 * np.arcsin(np.sqrt(p2)))
    except ValueError:
        return np.nan


def get_binom(row):
    k, n, p = row
    try:
        return binom_test(int(k), int(n), p)
    except ValueError:
        return np.nan


def get_exp_prop(taxid, database):
    try:
        return database[taxid].ratio
    except KeyError:
        return np.nan


def format_header(string, sample):
    return string + " (S{})".format(sample)


def get_expected_db(exp_summary_file, taxa):
    exp_sum = pd.read_csv(exp_summary_file, index_col=0)
    taxid_map = {k: v for k, v in zip(taxa, range(1, len(taxa) + 1))}
    exp_database = {}
    for tax, idx in taxid_map.items():
        tot, sig = exp_sum.loc[tax, get_cols(idx)]
        exp_database[tax] = EXP(tot, sig, sig/tot)
    return exp_database

def run_analysis(
    exp_summary_file, obs_summary_file, taxa, outfile):
    taxa = np.array(np.genfromtxt(taxa, dtype=int), ndmin=1)
    exp_db = get_expected_db(exp_summary_file, taxa)
    obs_sum = pd.read_csv(obs_summary_file)
    # remove all rows with non-candidate taxa
    obs_sum = obs_sum[obs_sum['TaxID'].isin(taxa)]
    obs_sum = process_stats(obs_sum, exp_db)
    obs_sum.to_csv(outfile, index=None)

if __name__ == "__main__":
    try:
        run_analysis(
            snakemake.input[0],
            snakemake.params[0],
            snakemake.params[1],
            snakemake.output[0]
        )
    except NameError:
        PARSER = argparse.ArgumentParser(
            prog="MTSv Random Kmers",
            description="Get random kmers from taxids.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )     
        PARSER.add_argument(
            "exp", metavar="EXP_SUMMARY", type=file_type,
            help="Path to expected values summary."
        )

        PARSER.add_argument(
            "obs", metavar="OBS_SUMMARY", type=file_type,
            help="Path to observed values summary."
        )

        PARSER.add_argument(
            "taxa", metavar="TAXA", type=file_type,
            help="Path to list of candidate taxa"
        )

        PARSER.add_argument(
            "outpath", metavar="OUTPATH", type=outfile_type,
            help="Output file"
        )

        ARGS = PARSER.parse_args()
        run_analysis(
            ARGS.exp,
            ARGS.obs,
            ARGS.taxa,
            ARGS.outpath)
    