from parsing import parse_config_sections, set_types

class Parameters:
    def __init__(self, defaults, args):
        self._defaults = defaults
        self._args = args

    @property
    def defaults(self):
        return self._defaults

    @defaults.setter
    def defaults(self, defaults):
        self._defaults = defaults
    
    @property
    def args(self):
        return self._args

    @args.setter
    def args(self, args):
        self._args = args

    @property
    def config(self):
        return parse_config_sections(
            self.args.config,
            self.args.cmd_class.config_section)

    @property
    def parameters(self):
        return self._parse_params()

    def _parse_params(self):
        '''change default params to reflect modifications
        from the config file and from the command line'''
        params = {key: value['default']
            for key, value in self.defaults.items()}
        types = {key: value['type']
            for key, value in self.defaults.items()}
        print("DEFAULTS", params)
        params.update(get_modified_params(self.defaults, self.config))
        print("Params CONFIG", params, "CONFIG", self.config)
        params.update(get_modified_params(self.defaults, vars(self.args)))
        print("Params ARGS", params)
        print("ARGS", self.args)
        return set_types(params, types)


    
    
def get_modified_params(defaults, args):
    '''return dictionary of arguments that have been
    modified/different from defaults'''
    mod_params = {}
    for k, v in args.items():
        if v is None:
            continue
        if k in defaults and v != defaults[k]['default']:
            mod_params[k] = v
        if k not in defaults:
            mod_params[k] = v
    return mod_params


        

   

    
