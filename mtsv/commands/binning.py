from __future__ import absolute_import
from .command import Command

class Binning(Command):
    config_section = ["BINNING"]
    # def __init__(self):
    #     print("Running binning")
    #     super().__init__()

    def run(self):
        print("PARAMS", self.params)
