from __future__ import absolute_import
from .command import Command

class Init(Command):
    def __init__(self):
        print("Running init")
        super().__init__()
    
    def run(self):
        from .__init__ import create_config_file

        with open(self.config, 'w') as config:
            config.write(create_config_file())

