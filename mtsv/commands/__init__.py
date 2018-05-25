import os
import logging
from mtsv.parsing import (
    create_config_file,
    format_commands
    )
from mtsv.commands.rules import (
    readprep,
    binner,
    summary,
    analyze,
    extract,
    setup_workflow
)


class Command:
    config_section = []

    def __init__(self, params):
        self._logger = logging.getLogger(__name__)
        self._params = params.params
        self._rules = []
        self._targets = []
        self._cml_args = []

    def __repr__(self):
        return str(self.__class__.__name__)

    @classmethod
    def get_config_line(cls):
        return cls.config_section

    @property
    def params(self):
        return self._params

    @property
    def rules(self):
        return self._rules

    @rules.setter
    def rules(self, rules):
        self._rules = rules

    @property
    def targets(self):
        return self._targets

    @targets.setter
    def targets(self, targets):
        self._targets = targets
        
    @property
    def cml_args(self):
        return self._cml_args

    @cml_args.setter
    def cml_args(self, cml_args):
        self._cml_args = cml_args

    def run(self):
        workflow = setup_workflow(self.params)
        for rule in self.rules:
            workflow = rule(self, workflow)
        workflow.check()
        print("EXECUTING")
        workflow.execute(nolock=True,
            targets=self.targets, updated_files=[],
            forceall=self.params['force'], resources={},
            dryrun=False)


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
        self.rules = [analyze]
        self.targets = [self.params['analysis_file']]


class Binning(Command):
    config_section = ["BINNING"]

    def __init__(self, params):
        print("Running binning")
        super().__init__(params)
        self.rules = [binner]
        index_names = [os.path.basename(p).split(
            ".")[0] + ".bn" for p in self.params['fm_index_paths']]
        # self.targets = [os.path.join(
        #     self.params['binning_outpath'], index_name) for index_name in index_names]
        self.targets = [os.path.join(
            self.params['binning_outpath'], "merged.clp")]
        self._mode_dict = {'fast': {'seed-size': 17, 'min-seeds': 5, 'seed-gap': 2},
                          'efficient': {'seed-size': 14, 'min-seeds': 4, 'seed-gap': 2},
                          'sensitive': {'seed-size': 11, 'min-seeds': 3, 'seed-gap': 1}}
        self._set_binning_mode()
        self.cml_args = format_commands(
            'BINNING',
            self.params,
            ['binning_mode', 'fm_index_paths', 'binning_outpath'], []
        )

    def _set_binning_mode(self):
        bin_mode = self._mode_dict[self.params['binning_mode']]
        self.params['seed-size'] = bin_mode['seed-size']
        self.params['min-seeds'] = bin_mode['min-seeds']
        self.params['seed-gap'] = bin_mode['seed-gap']

class Extract(Command):
    config_section=["EXTRACT"]

    def __init__(self, params):
        print("running extract")
        super().__init__(self, params)
        self.rules = [extract]
        self.targets = []


class Readprep(Command):
    config_section=["READPREP"]

    def __init__(self, params):
        print("running readprep")
        super().__init__(params)
        self.targets = [self.params['fasta_query_path']]
        self.rules = [readprep]     
        self.cml_args = format_commands(
            'READPREP',
            self.params,
            ['fasta_query_path', 'kmer', 'trim_mode'],
            ['--segment', str(self.params['kmer'])]
            if self.params['trim_mode'] == 'segment'
            else ['--{}'.format(self.params['trim_mode'])])


class Summary(Command):
    config_section=["SUMMARY"]
    def __init__(self, params):
        print("running summary")
        super().__init__(params)
        self.rules = [summary]
        self.targets = []


class Pipeline(Command):
    config_section=["READPREP", "BINNING", "SUMMARY", "ANALYZE", "PIPELINE"]

    def __init__(self, params):
        print("Running Pipeline")
        super().__init__(params)
        self.rules = [readprep, binner, summary, analyze]
        self.targets=[self.params['analysis_file']]




