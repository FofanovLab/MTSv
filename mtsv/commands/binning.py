from __future__ import absolute_import
from .command import Command

class Binning(Command):
    config_line = ["BINNING"]
    def __init__(self):
        print("Running binning")
        super().__init__()

    def run(self):
        pass
