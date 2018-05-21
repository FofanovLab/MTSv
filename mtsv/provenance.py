class Parameters:
    def __init__(self, params):
        self._params = vars(params)
    
    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, params):
        self._params = vars(params)



    



        

   

    
