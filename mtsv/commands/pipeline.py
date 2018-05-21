from __future__ import absolute_import
from .command import Command

class Pipeline(Command):
    config_section = ["READPREP", "BINNING", "SUMMARY", "ANALYZE", "PIPELINE"]

    # def __init__(self):
    #     print("Running Pipeline")
    #     super().__init__()

    def run(self):
        pass
