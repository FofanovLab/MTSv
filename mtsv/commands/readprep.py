from __future__ import absolute_import
from .command import Command

class Readprep:
    config_line = ["READPREP"]
    def __init__(self):
        print("running readprep")
        super().__init__()

    def run(self):
        pass
