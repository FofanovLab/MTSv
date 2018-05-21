import logging

class Command:
    config_section = []
    def __init__(self, params):
        self._logger = logging.getLogger(__name__)
        self._params = params.parameters
    
    def __repr__(self):
        return str(self.__class__.__name__)
        
    @classmethod
    def get_config_line(cls):
        return cls.config_section

    @property
    def params(self):
        return self._params


    

    


    
