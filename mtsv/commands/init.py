from __future__ import absolute_import
from .command import Command
from parsing import create_config_file

class Init(Command):
    config_section = []
    # def __init__(self, params):
    #     print("Running init")
    #     super().__init__(params)
    
    def run(self):
        print("running init")
        create_config_file(
            ["READPREP", "BINNING", "SUMMARY", "ANALYZE"],
            self.params['config'])
