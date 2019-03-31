from mtsv.parsing import outfile_type
import os

class Parameters:
    def __init__(self, params, snake_params):
        self._params = vars(params)
        self._snake_params = snake_params
    
    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, params):
        self._params = vars(params)

    @property
    def snake_params(self):
        return self._snake_params

    @snake_params.setter
    def snake_params(self, params):
        self._snake_params += params

    def write_parameters(self):
        file_name = outfile_type(os.path.join("./Params",
        "{cmd}_{timestamp}_params.txt".format(
            cmd=self.params['cmd_class'].__name__,
            timestamp=self.params['timestamp'])))
        with open(file_name, 'w') as out:
            for k, v in self.params.items():
                if hasattr(v, 'name'):
                    # file handles
                    v = v.name
                elif hasattr(v, '__name__'):
                    # class instance
                    v = v.__name__
                else:
                    # other
                    v = str(v)
                out.write("{0}: {1}\n".format(k, v))
            
        


    



        

   

    
