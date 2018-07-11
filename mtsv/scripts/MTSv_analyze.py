import argparse
import json
import logging
import numpy as np
import pandas as pd
from collections import namedtuple
from functools import partial
from scipy.stats import binom_test
from mtsv.parsing import file_type, outfile_type
from mtsv.utils import (
    get_precalculated_df, error, config_logging,
    write_to_precalculated_db)


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

def get_precalculated_values(precalc_file):
    precalc_dic = {}
    for _, values in json.loads(
        open(precalc_file, 'r').read()).items():
        precalc_dic[values['Taxid']] = EXP(
            values['Total_Hits'],
            values['Sig_Hits'],
            values['Ratio'])
    return precalc_dic


def get_expected_db(exp_summary_files, taxa):
    exp_database = {}
    all_taxa = [np.array([], dtype=int)]
    for exp_sum_file, taxa_file in zip(exp_summary_files, taxa):
        tax = np.array(np.genfromtxt(taxa_file, dtype=int), ndmin=1)
        exp_sum = pd.read_csv(exp_sum_file, index_col=0)
        taxid_map = {k: v for k, v in zip(tax, range(1, len(tax) + 1))}
        for t, idx in taxid_map.items():
            try:
                tot, sig = exp_sum.loc[t, get_cols(idx)]
                exp_database[t] = EXP(tot, sig, sig/tot)
            except KeyError as e:
                print(e)
        all_taxa = np.append(all_taxa, tax)
    return exp_database, all_taxa

def run_analysis(
    exp_summary_files, obs_summary_file, taxa, outfile, precalc_db=None):
    if len(exp_summary_files) != len(taxa):
        raise IOError(
            "Number of expected summary files must equal number of taxa files")
    exp_db, taxa = get_expected_db(exp_summary_files, taxa)
    new_db = {}
    if precalc_db:
        new_db = exp_db
        exp_db = {**exp_db, **precalc_db}
        taxa = np.append(taxa, list(precalc_db.keys()))
    obs_sum = pd.read_csv(obs_summary_file)
    # remove all rows with non-candidate taxa
    obs_sum = obs_sum[obs_sum['TaxID'].isin(taxa)]
    obs_sum = process_stats(obs_sum, exp_db)
    LOGGER.info("Writing to file {}".format(outfile))
    obs_sum.to_csv(outfile, index=None)
    return new_db

def return_non_nan(row):
    x, y = row
    if not np.isnan(x) and not np.isnan(y):
        return y
    elif np.isnan(x):
        return y
    else:
        return x


def update_db(new_db, **params):
    LOGGER.info("Updating expected value database")
    taxids = list(new_db.keys())
    total_hits, sig_hits, ratio = list(zip(*[
        (v.total, v.sig, v.ratio) for v in new_db.values()]))    
    new_df = pd.DataFrame({
        'Database': params['db_path'],
        'Kmer_Size': params['kmer'],
        'Edits': params['edits'],
        'Seed_Size': params['seed_size'],
        'Seed_Gap': params['seed_gap'],
        'Min_Seeds': params['min_seeds'],
        'Taxid': taxids,
        'Total_Hits': total_hits,
        'Sig_Hits': sig_hits,
        'Ratio': ratio }, columns=[
            'Database', 'Kmer_Size',
            'Edits', 'Seed_Size',
            'Seed_Gap', 'Min_Seeds',
            'Taxid', 'Total_Hits',
            'Sig_Hits', 'Ratio'], 
        dtype={
            'Database': str,
            'Kmer_Size': int,
            'Edits': int,
            'Seed_Size': int,
            'Seed_Gap': int,
            'Min_Seeds': int,
            'Taxid': int,
            'Total_Hits': int,
            'Sig_Hits': int,
            'Ratio': float
        })
    old_df = get_precalculated_df()
    header = new_df.columns
    res = old_df.merge(
        new_df, how="outer",
        on=['Database', 'Kmer_Size',
        'Edits', 'Seed_Size', 'Seed_Gap',
        'Min_Seeds', 'Taxid'])
    res['Total_Hits'] = res[["Total_Hits_x", "Total_Hits_y"]].apply(
        return_non_nan, axis=1)
    res['Sig_Hits'] = res[["Sig_Hits_x", "Sig_Hits_y"]].apply(
        return_non_nan, axis=1)
    res['Ratio'] = res[["Ratio_x", "Ratio_y"]].apply(
        return_non_nan, axis=1)
    res = res[header]
    write_to_precalculated_db(res)
    LOGGER.info("Finished updating expected value database")


    

if __name__ == "__main__":
    try:
        config_logging(snakemake.log[0], "INFO")
        LOGGER = logging.getLogger(__name__)
        PRECALC = get_precalculated_values(
            snakemake.params['taxa'])
        new_db = run_analysis(
            [snakemake.input[0]],
            snakemake.params[0],
            [snakemake.params['req_taxa']],
            snakemake.output[0],
            PRECALC
        )
        update_db(new_db, **snakemake.params['exp_db_params'])
        
    except NameError:
        PARSER = argparse.ArgumentParser(
            prog="MTSv Random Kmers",
            description="Get random kmers from taxids.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )     
        PARSER.add_argument(
            "--exp", metavar="EXP_SUMMARY", nargs='+', type=file_type,
            help="Path to expected values summary."
        )

        PARSER.add_argument(
            "obs", metavar="OBS_SUMMARY", type=file_type,
            help="Path to observed values summary."
        )

        PARSER.add_argument(
            "--taxa", metavar="TAXA", nargs='+', type=file_type,
            help="Path to list of candidate taxa"
        )

        PARSER.add_argument(
            "outpath", metavar="OUTPATH", type=outfile_type,
            help="Output file"
        )

        PARSER.add_argument(
            "--log", type=outfile_type, default="analyze.log",
            help="Path to log file"
        )

        ARGS = PARSER.parse_args()
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)
        run_analysis(
            ARGS.exp,
            ARGS.obs,
            ARGS.taxa,
            ARGS.outpath)
    