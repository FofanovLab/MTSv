import logging
from provenance import Parameters
class Command:
    config_line = []
    def __init__(self):
        self._logger = logging.getLogger(__name__)
    
    def __repr__(self):
        return str(self.__class__.__name__)
        
    def get_config_line(self):
        return self.config_line
    @property
    def parameters(self, **kwargs):
        parameters = Parameters(**kwargs)
        for key in kwargs:
            setattr(self, key, kwargs[key])

    
