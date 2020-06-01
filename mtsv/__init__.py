PROG_NAME = "MTSv"
VERSION = "2.0.0"

DEFAULT_CFG_FNAME = "mtsv.cfg"
DEFAULT_LOG_FNAME = "mtsv_{COMMAND}_{TIMESTAMP}.log"

CONFIG_STRING = """
# NOTE: Changes to the config file in the middle of the pipeline
# may force previous steps to be rerun.
#  
# ===============================================================
#
# READPREP: {readprep_description}
#
# ===============================================================

# {fastq_pattern_description}

 
fastq_pattern: {fastq_pattern_default}
 

# {fastp_params_description}

 
fastp_params: {fastp_params_default}
 

# {kmer_size_description}

 
kmer_size: {kmer_size_default}
 

# ===============================================================
#
# BINNING: {binning_description}
#
# ===============================================================

# {database_config_description}

 
database_config: {database_path}


# {edits_description}

 
edits: {edits_default}
 

# {binning_mode_description}

 
binning_mode: {binning_mode_default}


# {max_hits_description}


max_hits: {max_hits_default}

 

# {seed_size_description}
# Uncomment (remove "#" before name) to modify this parameter
# This will override the value set for the BINNING_MODE.

 
# seed_size: 
 

# {min_seeds_description}
# Uncomment (remove "#" before name) to modify this parameter
# This will override the value set for the BINNING_MODE.

 
# min_seeds:
 

# {seed_gap_description}
# Uncomment (remove "#" before name) to modify this parameter
# This will override the value set for the BINNING_MODE.

 
# seed_gap: 
 

# ===============================================================
# 
# ANALYSIS: {analyze_description}
#
# ===============================================================

# {filter_rule_description}

 
filter_rule: {filter_rule_default}
 

# {filter_value_description}

 
filter_value: {filter_value_default}
 

# {filter_column_description}

 
filter_column: {filter_column_default}
 

# {datastore_description}

 
datastore: {datastore}
 

# {sample_n_kmers_description}

 
sample_n_kmers: {sample_n_kmers_default}
 

# {alpha_description}

 
alpha: {alpha_default}
 

# {h_description}


h: {h_default}


# {figure_kwargs_description}


figure_kwargs: {figure_kwargs_default}
 

# ===============================================================
#
# EXTRACT: {extract_description}
#
# ===============================================================
    
# {extract_taxids_description}


extract_taxids: {extract_taxids_default}


"""


CLUSTER_CONFIG = """
__default__:
  cpus: '{threads}'
  mem: 5000
  log: '{log}.cluster'
  jobname: '{rule}'
  time: "30:00"

fastp:
  jobname: "fastp"

readprep:
  jobname: "readprep"

binning:
  jobname: "binning"
  mem: 32000
  cpus: 12
  time: "2:00:00"

collapse:
  jobname: "collapse"

init_taxdump:
  jobname: "init_taxdump"

summary:
  jobname: "summary"
  time: "1:00:00"
  mem: 20000
  cpus: 12

filter_candidate_taxa:
  jobname: "filter_candidate_taxa"

get_candidates_not_in_database:
  jobname: "get_candidates_not_in_database"

random_kmers:
  jobname: "random_kmers"

analyze_binning:
  jobname: "analyze_binning"
  mem: 30000
  time: "1:00:00"

analyze_collapse:
  jobname: "analyze_collapse"

update_datastore:
  jobname: "update_datastore"

analysis:
  jobname: "analysis"

analysis_figure:
  jobname: "analysis_figure"

analysis_html:
  jobname: "analysis_html"

extract:
  jobname: "extract"

unaligned_queries:
  jobname: "unaligned_queries"
  mem: 8000

"""
