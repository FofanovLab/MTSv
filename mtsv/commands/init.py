from __future__ import absolute_import
from .command import Command
from parsing import create_config_file, SECTIONS

class Init(Command):
    config_section = []
    # def __init__(self, params):
    #     print("Running init")
    #     super().__init__(params)
    
    def run(self):
        print("running init")
        print(self.params)
        create_config_file(
            SECTIONS,
            self.params['config'])
