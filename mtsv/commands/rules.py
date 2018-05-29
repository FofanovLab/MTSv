import os
import glob
from mtsv.utils import bin_path
from snakemake.workflow import Workflow, Rules, expand
import snakemake.workflow
from snakemake import shell
from snakemake.logging import setup_logger




def setup_workflow(params):
    # setup_logger()
    workflow = Workflow(
        "__file__",
        overwrite_workdir=params['working_dir'])
    if params['cluster_cfg'] is not None:
        workflow.cluster_cfg = params['cluster_cfg']
    snakemake.workflow.rules = Rules()
    snakemake.workflow.config = dict()

    return workflow
    
def binner(cmd, workflow):
    binner_bin = bin_path('mtsv-binner')
    collapse_bin = bin_path('mtsv-collapse')
    #### TEMP #####
    p = "~/Desktop/Fofanov_Projects/Repos/MTSv/ext/target/release/"
    binner_bin = p + "mtsv-binner"
    collapse_bin = p + "mtsv-collapse"
    
    _indices = []
    index2path = {}
    for _path in cmd.params['fm_index_paths']:
        _dir = os.path.dirname(os.path.abspath(_path))
        for _p in glob.iglob(os.path.abspath(_path)):
            idx = os.path.basename(_p).split(".")[0]
            _indices.append(idx)
            index2path[idx] = _dir
    
    binned_files = [os.path.join(cmd.params['binning_outpath'], idx + ".bn") 
                for idx in _indices]
    print("BINNED", binned_files)
    print("PATH", list(index2path.values())[0])

    @workflow.rule(name='binning')
    @workflow.docstring("""Metagenomics binning""")
    @workflow.input(os.path.join(list(index2path.values())[0], "{index}.index"))
    @workflow.output(os.path.join(cmd.params['binning_outpath'], "{index}.bn"))
    @workflow.params(call=binner_bin, args=cmd.cml_args)
    @workflow.message("Executing Binner with {threads} threads on the following files {input}.")
    @workflow.log(os.path.join(cmd.params['binning_outpath'],"{index}.log"))
    @workflow.threads(cmd.params['threads'])
    @workflow.shellcmd("{params.call} --index {input} --threads {threads} --results {output} {params.args} > {log} 2>&1")
    @workflow.run
    def __rule_binning(input, output, params, wildcards,
                       threads, resources, log, version, rule,
                       conda_env, singularity_img, singularity_args,
                       use_singularity, bench_record, jobid, is_shell):
        shell(
            "{params.call} --index {input} --threads {threads} --results {output} {params.args} > {log} 2>&1")
    ####################
    # Collapse 
    ####################

    @workflow.rule(name='collapse')
    @workflow.docstring("""Combine the output of multiple separate mtsv runs. """)
    @workflow.message("Merging binning output files into {output}.")
    @workflow.input(binned_files)
    @workflow.output(cmd.params['merge_file'])
    @workflow.log(cmd.params['log_file'])
    @workflow.params(call=collapse_bin)
    @workflow.shellcmd("{params.call} {input} -o {output} > {log} 2>&1")
    @workflow.run
    def __rule_collapse(input, output, params, wildcards,
                        threads, resources, log, version, rule,
                        conda_env, singularity_img, singularity_args,
                        use_singularity, bench_record, jobid, is_shell):
        shell(
            "{params.call} {input} -o {output} > {log} 2>&1")

    return workflow


def readprep(cmd, workflow):
    readprep_bin = bin_path('mtsv-readprep')
    ##### TEMP ######
    p = "~/Desktop/Fofanov_Projects/Repos/MTSv/ext/target/release/"
    readprep_bin = p + "mtsv-readprep"

    @workflow.rule(name='readprep')
    @workflow.docstring("""Read fragment quality control and deduplication (FASTQ -> FASTA)""")
    @workflow.input(cmd.params['fastq'])
    @workflow.output(cmd.params['fasta_query_path'])
    @workflow.params(call=readprep_bin, args=cmd.cml_args)
    @workflow.log(cmd.params['log_file'])
    @workflow.threads(cmd.params['threads'])
    @workflow.shellcmd("{params.call} {input} -o {output} {params.args} >> {log} 2>&1")
    @workflow.run
    def __rule_readprep(input, output, params, wildcards,
                        threads, resources, log, version, rule,
                        conda_env, singularity_img, singularity_args,
                        use_singularity, bench_record, jobid, is_shell):
        shell(
            "{params.call} {input} -o {output} {params.args} >> {log} 2>&1")
    return workflow


def summary(cmd, workflow):
    signature_bin = bin_path('mtsv-signature')
    
    ###### TEMP ######
    p = "~/Desktop/Fofanov_Projects/Repos/MTSv/ext/target/release/"
    signature_bin = p + "mtsv-signature"

    @workflow.rule(name='signature')
    @workflow.docstring("""Find signature hits""")
    @workflow.input(
        infile=cmd.params['merge_file'],
        index=cmd.params['tree_index'])
    @workflow.output(cmd.params['signature_file'])
    @workflow.params(call=signature_bin, args=cmd.cml_args)
    @workflow.log(cmd.params['log_file'])
    @workflow.threads(cmd.params['threads'])
    @workflow.shellcmd(
        "{params.call} --index {input.index} --input {input.infile} --output {output} {params.args} >> {log} 2>&1")
    @workflow.run
    def __rule_signature(input, output, params, wildcards,
                        threads, resources, log, version, rule,
                        conda_env, singularity_img, singularity_args,
                        use_singularity, bench_record, jobid, is_shell):
        shell(
            "{params.call} --index {input.index} --input {input.infile} --output {output} {params.args} >> {log} 2>&1")

    return workflow

def analyze(cmd, workflow):
    return workflow


def extract(cmd, workflow):
    return workflow




    ####################
    # READPREP WORKFLOW
    ####################

    
    ####################
    # BINNING WORKFLOW
    ####################
    # # Make full list of index paths
    
