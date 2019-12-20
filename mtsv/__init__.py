PROG_NAME = "MTSv"
VERSION = "1.1.0"

DEFAULT_CFG_FNAME = "mtsv.cfg"
DEFAULT_LOG_FNAME = "mtsv_{COMMAND}_{TIMESTAMP}.log"

CONFIG_STRING = """
# NOTE: Changes to the config file in the middle of the pipeline
# may force previous steps to be rerun. 
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
 

# {seed_size_description}
# Uncomment (remove "#" before name) to modify this parameter

 
# seed_size: 
 

# {min_seeds_description}
# Uncomment (remove "#" before name) to modify this parameter

 
# min_seeds:
 

# {seed_gap_description}
# Uncomment (remove "#" before name) to modify this parameter

 
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


# ===============================================================
#
# WGFAST: {wgfast_description}
#
# ===============================================================

# {wgfast_taxids_description}


wgfast_taxids: {wgfast_taxids_default}


# {wgfast_reference_dirs_description}


wgfast_reference_dirs: {wgfast_reference_dirs_default}


# {wgfast_params_description}


wgfast_params: {wgfast_params_default}


# ===============================================================
#
# CONCOCT: {concoct_description}
#
# ===============================================================

# {megahit_params_description}


megahit_params: {megahit_params_default}


# {cutup_params_description}


cutup_params: {cutup_params_default}


# {bwa_params_description}


bwa_params: {bwa_params_default}


# {concoct_params_description} 


concoct_params: {concoct_params_default}


"""


CLUSTER_CONFIG = """
__default__:
  cpus: '{threads}'
  mem: 5000
  log: 'logs/cluster_{log}'
  jobname: '{rule}'
  time: "30:00"

fastp:
  jobname: "fastp"

readprep:
  jobname: "readprep"

binning:
  jobname: "binning"
  mem: 30000
  time: "2:00:00"

collapse:
  jobname: "collapse"

summary:
  jobname: "summary"
  time: "1:00:00"
  mem: 8000

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

wgfast:
  jobname: "wgfast"
  time: "1:00:00"
  mem: 8000

wgfast_draw_tree:
  jobname: "wgfast_draw_tree"

wgfast_fasta_to_fastq:
  jobname: "wgfast_fasta_to_fastq"

wgfast_gzip:
  jobname: "wgfast_gzip"

megahit:
  jobname: "megahit"
  time: "1:00:00"
  mem: 8000

concoct_cut_up:
  jobname: "concoct_cut_up"

concoct_fasta2fastq:
  jobname: "concoct_fasta2fastq"

concoct_make_bam:
  jobname: "concoct_make_bam"

concoct_coverage_table:
  jobname: "concoct_coverage_table"

concoct:
  jobname: "concoct"

concoct_merge:
  jobname: "concoct_merge"

concoct_fasta_bins:
  jobname: "concoct_fasta_bins"

"""
