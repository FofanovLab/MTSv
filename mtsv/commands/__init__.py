import subprocess as sp
from mtsv.utils import snake_path, get_database_params, error, warn



class Command:
    config_section = []

    def __init__(self, params):
        self._params = params
        self._rules = []

    def __repr__(self):
        return str(self.__class__.__name__)

    @classmethod
    def get_config_line(cls):
        return cls.config_section

    @property
    def params(self):
        return self._params.params
    
    @property
    def snake_params(self):
        return self._params.snake_params
    
    @snake_params.setter
    def snake_params(self, params):
        self._params.snake_params = params

    @property
    def rules(self):
        return self._rules

    @rules.setter
    def rules(self, rules):
        self._rules = rules

    def run(self):
        for rule in self.rules:
            cmd = ["snakemake", "--snakefile", rule, "--config"]
            config = [
                "{0}={1}".format(k,v)
                for k, v in self.params.items() if v is not None]
            cmd += config
            cmd += self.snake_params
            try:
                p = sp.run(cmd,
                        check=True)
                self._params.write_parameters()
            except ( KeyboardInterrupt, sp.CalledProcessError) as e:
                warn("Unlocking directory after failed snakemake")
                sp.run(cmd + ["--unlock"], check=True )
                error(e)
            





