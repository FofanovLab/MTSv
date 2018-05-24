import configparser
import os
from functools import wraps
from snakemake.workflow import Workflow, Rules
import snakemake.workflow
from snakemake import shell
from snakemake.logging import setup_logger, Logger
import logging
from pkg_resources import resource_filename

from mtsv.parsing import (
    create_config_file,
    specfile_read,
    format_commands
    )


class Command:
    config_section = []
    targets = []

    def __init__(self, params):
        self._logger = logging.getLogger(__name__)
        self._params = params.params

    def __repr__(self):
        return str(self.__class__.__name__)

    @classmethod
    def get_config_line(cls):
        return cls.config_section

    @property
    def get_targets(self):
        return self.targets

    @property
    def params(self):
        return self._params

    def run(self):
        workflow = self.snek()
        workflow.check()
        print("EXECUTING")
        workflow.execute(nolock=True,
            targets=self.targets, updated_files=[],
            forceall=self.params['force'], resources={},
            dryrun=False)
        
    def snek(self):
        readprep = bin_path('mtsv-reaprep')
        binner = bin_path('mtsv-binner')
        collapse = bin_path('mtsv-collapse')
        signature = bin_path('mtsv-signature')
        # temp
        p = "~/Desktop/Fofanov_Projects/Repos/MTSv/ext/target/release/"
        readprep = p + "mtsv-readprep"
        binner = p + "mtsv-binner"
        collapse = p + "mtsv-collapse"
        signature = p + "mtsv-signature"

        setup_logger()
        workflow = Workflow(
            "__file__",
            overwrite_workdir=self.params['working_dir'])
        if self.params['cluster_cfg'] is not None:
            workflow.cluster_cfg = self.params['cluster_cfg']
        snakemake.workflow.rules = Rules()
        snakemake.workflow.config = dict()
        @workflow.rule(name='readprep')
        @workflow.docstring("""Read fragment quality control and deduplication (FASTQ -> FASTA)""")
        @workflow.input(self.params['fastq'])
        @workflow.output(self.params['query_fasta_path'])
        @workflow.params(cmd=readprep, params=self.cml_params)
        @workflow.log(self.params['log_file'].name)
        @workflow.threads(self.params['threads'])
        @workflow.shellcmd("{params.cmd} {input} -o {output} {params.params} >> {log} 2>&1")
        
        @workflow.run
        def __rule_readprep(input, output, params, wildcards, 
            threads, resources, log, version, rule,
            conda_env, singularity_img, singularity_args,
            use_singularity, bench_record, jobid, is_shell):
            shell(
                "{params.cmd} {input} -o {output} {params.params} >> {log} 2>&1")
        return workflow
# "{readprep} {input} -o {output} --threads {threads} {readprep_params}"
class Init(Command):
    config_section = []

    def __init__(self, params):
        print("Running init")
        super().__init__(params)

    def run(self):
        print("running init")
        create_config_file(
            self.params['config'])


class Analyze(Command):
    config_section = ["ANALYZE"]

    def __init__(self, params):
        print("Running analysis")
        super().__init__(params)
        self.targets = [self.params['analysis_file']]


class Binning(Command):
    config_section = ["BINNING"]

    def __init__(self, params):
        print("Running binning")
        super().__init__(params)
        self.targets = [join(
            self.params['binning_outpath'],
            '{index}_binned.bin')]


class Extract(Command):
    config_section=["EXTRACT"]

    def __init__(self, params):
        print("running extract")
        super().__init__(self, params)
        self.targets = []


class Readprep(Command):
    config_section=["READPREP"]

    def __init__(self, params):
        print("running readprep")
        super().__init__(params)
        self.targets = [self.params['query_fasta_path']]
        
        self.cml_params = format_commands(
            'READPREP',
            self.params,
            ['query_fasta_path', 'kmer', 'trim_mode'],
            ['--segment', str(self.params['kmer'])]
            if self.params['trim_mode'] == 'segment'
            else ['--{}'.format(self.params['trim_mode'])])
        


class Summary(Command):
    config_section=["SUMMARY"]
    def __init__(self, params):
        print("running summary")
        super().__init__(params)
        self.targets = []


class Pipeline(Command):
    config_section=["READPREP", "BINNING", "SUMMARY", "ANALYZE", "PIPELINE"]

    def __init__(self, params):
        print("Running Pipeline")
        super().__init__(params)
        self.targets=[self.params['analysis_file']]



def bin_path(cmd):
    """Return the specfile path for a given command name."""
    fp=os.path.join('ext', 'cmd')
    return resource_filename('mtsv', fp)
