from __future__ import absolute_import
from .command import Command

class Pipeline(Command):
    config_line = ["READPREP", "BINNING", "SUMMARY", "ANALYSIS", "PIPELINE"]

    def __init__(self):
        print("Running Pipeline")
        super().__init__()

    def run(self):
        pass
