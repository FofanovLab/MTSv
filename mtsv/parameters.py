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

    def write_parameters(self, file_name):
        with open(file_name, 'w') as out:
            for k, v in self.params.items():
                try:
                    # get file name from file handles
                    v = v.name
                except AttributeError:
                    v = str(v)
                out.write("{0}: {1}\n".format(k, v))
            
        


    



        

   

    
