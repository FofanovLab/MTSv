from __future__ import absolute_import
from .command import Command

class Analyze(Command):
    config_section = ["ANALYZE"]
    # def __init__(self):
    #     print("Running analysis")
    #     super().__init__()

    def run(self):
        pass
