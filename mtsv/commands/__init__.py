import os
import subprocess as sp
from mtsv.parsing import create_config_file
from mtsv.utils import snake_path

SNAKEFILES = {
    cmd:snake_path("{}_snek").format(cmd)
    for cmd in ["binning", "readprep", "summary", "analyze", "pipeline", "extract", "wgfast"]}

class Command:
    config_section = []

    def __init__(self, params):
        self._params = params
        self._rules = []

    def __repr__(self):
        return str(self.__class__.__name__)

    @classmethod
    def get_config_line(cls):
        return cls.config_section

    @property
    def params(self):
        return self._params.params
    
    @property
    def snake_params(self):
        return self._params.snake_params

    @property
    def rules(self):
        return self._rules

    @rules.setter
    def rules(self, rules):
        self._rules = rules

    def run(self):
        for rule in self.rules:
            cmd = ["snakemake", "--snakefile", rule, "--config"]
            config = [
                "{0}={1}".format(k,v)
                for k, v in self.params.items() if v is not None]
            cmd += config
            cmd += ["--restart-times", "3", "--nolock"]
            cmd += self.snake_params
            try:
                p = sp.run(cmd,
                        # stdout = sp.PIPE,
                        # stderr=sp.PIPE,
                        check=True)
                self._params.write_parameters()
            except sp.CalledProcessError as e:
                print(e)
            

class Init(Command):

    def __init__(self, params):
        super().__init__(params)

    def run(self):
        create_config_file(
            self.params['config'])


class Analyze(Command):
    config_section = ["ANALYZE"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['analyze']]


class Binning(Command):
    config_section = ["BINNING"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['binning']]

class Extract(Command):
    config_section=["EXTRACT"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['extract']]

class Readprep(Command):
    config_section=["READPREP"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['readprep']]     

class Summary(Command):
    config_section=["SUMMARY"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['summary']]


class Pipeline(Command):
    config_section=["READPREP", "BINNING",
        "SUMMARY", "ANALYZE", "PIPELINE"]

    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['pipeline']]

class WGFast(Command):
    def __init__(self, params):
        super().__init__(params)
        self.rules = [SNAKEFILES['wgfast']]


