from __future__ import absolute_import
from .command import Command

class Summary:
    config_line = ["SUMMARY"]
    def __init__(self):
        print("running summary")
        super().__init__(self)
    
    def run(self):
        pass
