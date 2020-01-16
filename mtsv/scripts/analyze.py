import hashlib
import logging
import datetime
import os
import pandas as pd
import numpy as np

import json
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
from statsmodels.stats.proportion import proportions_ztost
from mtsv.scripts.summary import (
    LEVELS,
    get_chunked_reader, add_taxa_count,
    drop_rows_over_max_taxa, Lineage,
    get_signature_mask)
from mtsv.utils import (
    get_ete_ncbi, warn, template_path)


def get_species(df):
    return df[df['level'] == "species"]


def signature_taxa_cutoff(df, cutoff):
    """
    Return rows with unique signature hits greater
    than cutoff.
    """
    return df[df['unique_signature'] > cutoff]


def quantile_rule(col, value):
    def df_operation(df):
        return df[df[col] > df[col].quantile(value)]['taxid'].values
    return df_operation


def min_rule(col, value):
    def df_operation(df):
        return df[df[col] >= value]['taxid'].values
    return df_operation


def in_rule(col, value):
    if isinstance(value, int):
        value = [value]

    def df_operation(df):
        return df[df['taxid'].isin(value)]['taxid'].values
    return df_operation


def top_rule(col, value):
    def df_operation(df):
        return df.nlargest(value, col)['taxid'].values
    return df_operation


def filter_candidate_taxa(df, rule, column, value):
    """
    Input:
        df: (dataframe) from summary file 
        rule: (str) name of rule in {top, min, in, quant}
        col: (str) name of column to filter on
        value: (list or int) to apply to rule
    Output:
        List of taxids that pass filter.
    """
    rule_dict = {'top': top_rule,
                 'min': min_rule,
                 'in': in_rule,
                 'quant': quantile_rule}
    filter_func = rule_dict[rule](column, value)
    return list(filter_func(df))



def get_candidate_taxa(fn, rule, col, value):
    df = get_species(
            pd.read_csv(fn))
    return filter_candidate_taxa(df, rule, col, value)


def analyze(summary_data, expected_data, alpha, h, filter_params):
    """
    Performes statistical analysis on taxids that have expected data.
    Returns flag column to indicate significantly similar values
    for expected and observed proportions and a flag column that 
    indicates taxids that pass the initial user-defined filter_params.

    Input:
        summary_data: Dataframe with taxids as index and summary info
        expected_data: Dataframe with values for expected proportion for
                        taxids (index) that have been sampled.
        alpha: (float [0-1]) Cutoff value for statistical significance
        h: (float [0-1]) Cohen's h value for effect size. This determines
            the range for the equivalence test.
        filter_params: (dict) dictionary of filter params 
                        (keys: rule, column, value)
    Output:
        Dataframe that contains columns: "unique_signature_exp",
        "unique_signature_obs", "unique_exp", "unique_obs", "prop_exp",
        "prop_obs", "p-value", "significant", and "user_filter"

    """
    # Remove any non-signature taxa from summary
    summary_data = signature_taxa_cutoff(summary_data, 0)

    summary_data = summary_data.join(
        expected_data, lsuffix="", rsuffix="_exp")
    # Calculate proportions
    summary_data = get_proportions(
        get_proportions(summary_data, "_exp"), "")
    # Calculate p-value and perform equivalence test
    summary_data = hypothesis_test(get_pvalue(summary_data, h), alpha)
    summary_data = user_filter(summary_data, filter_params)
    return format_summary_data(summary_data)

def format_summary_data(summary_data):
    cols = ["taxid","scientific_name", "level", "total", "unique",
    "signature", "unique_signature",
    "weighted_support", "unique_weighted_support", "unique_exp",
            "unique_signature_exp", "prop_exp", "prop", "p-value",
    "significant", "user_filter"
    ] + LEVELS[::-1]
    summary_data = summary_data[cols]
    summary_data.columns = [
        "TaxID", "Scientific_Name", "Level", "Total_Hits", "Unique_Hits",
        "Signature_Hits", "Unique_Signature_Hits", "Weighted_Support",
        "Unique_Weighted_Support", "Exp_Unique_Hits",
        "Exp_Unique_Signature_Hits",
        "Exp_Prop(USH/UH)", "Obs_Prop(USH/UH)",
        "P_value", "Significant", "User_Filter"] + \
            [level.capitalize() for level in LEVELS[::-1]]
    return summary_data
    
def get_proportions(df, suffix):
    # To be in summary there must be a positive number of unique hits
    # so there should not be a divide by zero error
    # If there is they will be converted to nan
    df["prop{0}".format(suffix)] = \
        df['unique_signature{0}'.format(suffix)] / \
        df['unique{0}'.format(suffix)]
    return df


def get_upper_tost_bound(p1, h):
    upper = np.sin((h - (2 * np.arcsin(np.sqrt(p1))))/-2)**2
    if upper > p1:
        return p1
    return p1 - upper

def get_lower_tost_bound(p1, h):
    lower = np.sin((h + (2 * np.arcsin(np.sqrt(p1))))/2)**2
    if lower < p1:
        lower = 1
    return p1 - lower


def p_value_apply(row, h):
    """
    For row with ids: "unique_signature_exp", "unique_signature",
    "unique_exp", "unique", "prop_exp",
    calculate proportion_ztost, using upper and lower bounds determined
    using cohen's h value. The variance of the proportion estimate
    is calculated based on the sample proportion.
    """
    y1 = row['unique_signature_exp']
    y2 = row['unique_signature']
    n1 = row['unique_exp']
    n2 = row['unique']
    p1 = row['prop_exp']

    return proportions_ztost(
        [y1, y2], [n1, n2], get_lower_tost_bound(p1, h),
        get_upper_tost_bound(p1, h),
        prop_var="sample")[0]


def get_pvalue(df, h):
    df['p-value'] = df.apply(
        lambda row: p_value_apply(row, h), axis=1)
    return df

def hypothesis_test(df, alpha):
    """
    Adds "significant" column to df.
    String True if value in "p-value" column
    is <= alpha else False.
    Empty String if there is no value for
    "p-value"
    """
    df['significant'] = df['p-value'].apply(
        lambda x: "" if np.isnan(x) else str(x <= alpha) )
    return df

def user_filter(summary, filter_params):
    if 'taxid' not in summary.columns:
        summary['taxid'] = summary.index
    passed = filter_candidate_taxa(
        summary, **filter_params)
    summary['user_filter'] = False
    summary.loc[passed, 'user_filter'] = True
    return summary


def combine_species_and_genus(df):
    """
    Adds species that are significant. If none of the species
    in a genus are significant, a significant genus will be added.
    """
    species = df[(df['Level'] == "species") & (df['Significant']) & 
        (df['User_Filter'])].copy()
    genera_in_species = np.array(species['Genus'].unique(), dtype='int64')
    genera = df[(df['Level'] == "genus") & (df['Significant']) &
        (df['User_Filter']) & (
        ~df['TaxID'].isin(genera_in_species))].copy()
    return species.append(genera)


def calculate_abundance(df):
    df['abundance'] = df['Signature_Hits'] / \
        df['Unique_Signature_Hits']
    return df


def cat_analysis_files(analysis_files):
    df = pd.DataFrame([])
    for analysis_file in analysis_files:
        sample = os.path.basename(analysis_file).rstrip("_analysis.csv")
        next_df = pd.read_csv(analysis_file)
        next_df['sample'] = sample
        next_df = combine_species_and_genus(next_df)
        next_df = calculate_abundance(next_df)
        df = df.append(next_df)
    return df

def duplicate_column(df):
    dup_df = df.copy()
    dup_df['sample'] = df['sample'].apply(lambda x: x + "_dup")
    return df.append(dup_df)


def duplicate_row(df):
    dup_df = df.copy()
    dup_df['Scientific_Name'] = df['Scientific_Name'].apply(
        lambda x: x + "_dup")
    return df.append(dup_df)


def get_pivot_table(df):
    data = df.pivot(columns="sample", index="Scientific_Name",
                    values="abundance")
    data['avg'] = data.apply(lambda row: np.mean(row), axis=1)
    # place most abundant across all samples near top
    data = data.sort_values(by="avg", ascending=False)
    data = data.drop("avg", axis=1)
    data = data.replace(np.nan, 0)
    return data


def draw_figure(df, kwargs):
    sns.set_style("whitegrid")
    sns.set_context("paper")
    cmap = sns.cubehelix_palette(8, light=1, as_cmap=True)
    default_kwargs = {
        'cmap': cmap, "vmin": .8,
        "yticklabels": 1, "xticklabels": 1}
    default_kwargs.update(kwargs)
    g = sns.clustermap(df, **default_kwargs)
    ax = g.ax_heatmap
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    ax.set_xlabel("Sample")
    ax.set_ylabel("")
    return g



def heatmap_figure(analysis_files, output, table_output, kwargs):
    df = cat_analysis_files(analysis_files)
    if len(df['sample'].unique()) == 1:
        msg = """
        Only one sample, a duplicate sample is added for clustermap
        to work.
        """
        logging.warn(msg)
        warn(msg)
        df = duplicate_column(df) # need to have more than one column to work
    if len(df['Scientific_Name'].unique()) == 1:
        msg = """
        Only one taxa is significant, a duplicate of this taxa
        is added for clustermap to work.
        """
        logging.warn(msg)
        warn(msg)
        df = duplicate_row(df) # need to have more than one row to work
    df = get_pivot_table(df)
    df.to_csv(table_output)
    fig = draw_figure(df, kwargs)
    fig.savefig(output, bbox_inches="tight")


JSON_METADATA = {
"meta": {
    "ranks": {
        "Superkingdom": {"level": 1, "color": "blue"},
     			"Phylum": {"level": 2, "color": "blue lighten-1"},
     			"Class": {"level": 3, "color": "blue lighten-2"},
     			"Order": {"level": 4, "color": "blue lighten-3"},
     			"Family": {"level": 5, "color": "blue lighten-4"},
     			"Genus": {"level": 6, "color": "blue lighten-5"},
     			"Species": {"level": 7, "color": "white"},
    }
},
"name": "root",
"values": [
    {"header": "Scientific Name", "propName": "Scientific_Name",
        "colWidth": 2, "colIndex": 0, "tooltipped": False, "tooltip": ""},
    {"header": "Taxa ID", "propName": "TaxID", "colWidth": 1,
        "colIndex": 1, "tooltipped": False, "tooltip": ""},
    {"header": "Level", "propName": "Level", "colWidth": 1,
        "colIndex": 2, "tooltipped": False, "tooltip": ""},
    {"header": "Total Hits", "propName": "Total_Hits", "colWidth": 1,
        "colIndex": 3, "tooltipped": True,
        "tooltip": "Total query hits (includes counts for queries that occur multiple times)"},
  	{"header": "Unique Hits",
        "propName": "Unique_Hits", "colWidth": 1, "colIndex": 4,
        "tooltipped": True,
        "tooltip": "Number of unique query hits."},
  	{"header": "Signature Hits",
        "propName": "Signature_Hits", "colWidth": 1, "colIndex": 5,
        "tooltipped": True,
        "tooltip": "Total queries (including counts for duplicate queries) where taxid was the only hit."},
  	{"header": "Unique Signature Hits",
        "propName": "Unique_Signature_Hits", "colWidth": 1, "colIndex": 6,
        "tooltipped": True,
        "tooltip": "Unique queries where taxid was the only hit."},
    {"header": "Expected Proportion",
        "propName": "Exp_Prop(USH/UH)", "colWidth": 1, "colIndex": 7,
        "tooltipped": True,
        "tooltip": "Expected proportion (from simulated data) of unique signature hits to unique hits."},
  	{"header": "Observed Proportion",
        "propName": "Obs_Prop(USH/UH)", "colWidth": 1, "colIndex": 8,
        "tooltipped": True,
        "tooltip": "Observed sample proportion of unique signature hits to unique hits."},
   	{"header": "P Value", "propName": "P_value", "colWidth": 1, "colIndex": 9,
        "tooltipped": True,
        "tooltip": "P-value from equivalence test of observed and expected proportions with ranges determined by user-specified Cohen's H value."},
  	{"header": "Significant",
        "propName": "Significant", "colWidth": 1, "colIndex": 10,
        "tooltipped": True,
        "tooltip": "Result of hypothesis test based on user-specified alpha value."}
]}


def get_report_html(analysis_df, output, sample):
    
    df = pd.read_csv(analysis_df,
                     dtype={"TaxID": int, "P_values": str,
                     "Exp_Prop(USH/UH)": str, "Obs_Prop(USH/UH)": str})
    logging.info(df.head())
    df = df.replace(np.nan, '', regex=True)
    lineage_table = get_lineage_table(df)
    logging.info(df.head())
    logging.info(df.dtypes)
    logging.info(lineage_table)
    children = recursive_levels(df, lineage_table, "root", None)

    children['meta'] = JSON_METADATA['meta']
    children['name'] = JSON_METADATA['name']
    children['values'] = JSON_METADATA['values']
    children_json = json.dumps(children)
    logging.info(children_json)
    template = Template(
        open(template_path('report_template.html'), 'r').read())
    with open(output, 'w') as out:
        out.write(template.render(
            children_data=children_json,
            sample_name=sample, creation_date=str(datetime.datetime.now())))


    


def get_lineage_table(df):
    lineage_table = df[df.Level == "species"][[l.capitalize() for l in LEVELS]].copy()
    lineage_table.columns = LEVELS
    return lineage_table



def make_taxa_obj(df, taxid, children):

    data_cols = ["Scientific_Name", "TaxID", "Level",
                     "Total_Hits", "Unique_Hits",
                     "Signature_Hits", "Unique_Signature_Hits",
                     "Exp_Prop(USH/UH)", "Obs_Prop(USH/UH)",
                     "P_value", "Significant"]
    data = df[df.TaxID == int(taxid)].loc[:, data_cols].iloc[0].to_dict()
    logging.info(data)
    data['Level'] = data['Level'].capitalize()
    data['P_value'] = "" if data['P_value'] == "" else data['P_value']
    for k in [
        'TaxID', "Total_Hits", "Unique_Hits",
        "Signature_Hits", "Unique_Signature_Hits",
            "Exp_Prop(USH/UH)", "Obs_Prop(USH/UH)"]:
        data[k] = str(data[k])
    return {
        'values': data,
        'children': children
    }


def recursive_levels(df, lineage_table, level, taxid):
    tree_levels = ['root'] + LEVELS[::-1]
    if level == "species":
        return make_taxa_obj(df, taxid, children=[])
    next_level = tree_levels[tree_levels.index(level) + 1]
    if level == "root":
        unique_children = list(map(int, lineage_table['superkingdom'].unique()))
        return {
            "children": [
                recursive_levels(
                    df, lineage_table, next_level, tax)
                    for tax in unique_children]}
    else:
        unique_children = list(
            map(
                int,
                lineage_table[lineage_table[level] == taxid][
                    next_level].unique()))
        # Remove values that are carried up the tree
        if taxid in unique_children:
            unique_children.remove(taxid)
    return make_taxa_obj(
        df, taxid,
        children=[recursive_levels(df, lineage_table, next_level, tax)
        for tax in unique_children])
