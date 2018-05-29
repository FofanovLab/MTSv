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
        workflow.execute(nolock=True,
            targets=self.targets, updated_files=[],
            forceall=self.params['force'], resources={},
                         dryrun=False, printshellcmds=True)

class Database(Command):
    config_section = ["DATABASE"]

    def __init__(self, params):
        print("running Database")
        super().__init__(params)
        print(self.params)

class CustomDB(Command):
    
    config_section = ["CUSTOM_DB"]

    def __init__(self, params):
        print("running custom database")
        super().__init__(params)
        print(self.params)

class WGFast(Command):
    config_section = ["WGFast"]

    def __init__(self, params):
        print("Running WGFast")
        super().__init__(params)
        print(self.params)

class Init(Command):

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
        self.targets = [os.path.join(
            self.params['binning_outpath'], "merged.clp")]
        self._mode_dict = {'fast': {'seed_size': 17, 'min_seeds': 5, 'seed_gap': 2},
                          'efficient': {'seed_size': 14, 'min_seeds': 4, 'seed_gap': 2},
                          'sensitive': {'seed_size': 11, 'min_seeds': 3, 'seed_gap': 1}}
        self._set_binning_mode()
        self.cml_args = format_commands(
            'BINNING',
            self.params,
            ['binning_mode', 'fm_index_paths', 'merge_file', 'binning_outpath'], []
        )


    def _set_binning_mode(self):
        bin_mode = self._mode_dict[self.params['binning_mode']]
        for key, value in bin_mode.items():
            if self.params[key] == None:
                self.params[key] = value
            # change back to dash instead of underscore

            self.params[key.replace("_", "-")] = self.params[key]
            del self.params[key]


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
        # self.targets = [self.params['summary_file']]
        self.targets = [self.params['signature_file']]
        self.cmd_args = format_commands(
            'SUMMARY',
            self.params,
            ['taxdump_path', 'tax_level', 'merge_file', 'signature_file', 'summary_file', 'tree_index'],
            ["--{}".format(self.params['tax_level'])] if
                self.params['tax_level'] != 'species' else [])
        print("CMDS", self.cmd_args)


class Pipeline(Command):
    config_section=["READPREP", "BINNING", "SUMMARY", "ANALYZE", "PIPELINE"]

    def __init__(self, params):
        print("Running Pipeline")
        super().__init__(params)
        self.rules = [readprep, binner, summary, analyze]
        self.targets=[self.params['analysis_file']]




