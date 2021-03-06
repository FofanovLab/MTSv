import argparse
import configparser
import os
import sys
import ast
import logging
from glob import glob
from argutils import (read, export)

from mtsv.utils import(error, warn, specfile_read)
from mtsv import (DEFAULT_LOG_FNAME, DEFAULT_CFG_FNAME)

logger = logging.getLogger(__name__)

def make_sub_parser(subparser, cmd, cmd_class):
    global_defaults = get_global_config(cmd_class.config_section)
    help_str = global_defaults["_meta_{}".format(cmd)]["help"]
    p = subparser.add_parser(
        cmd, help=help_str,
        description="Additional Snakemake commands may also be provided")
    for arg, desc in global_defaults.items():
        if "_meta" in arg:
            continue
        if 'type' in desc:
            desc['type'] = TYPES[desc['type']]
        if 'choices' in desc and 'help' in desc:
            desc['help'] += " Choices are {}".format(desc['choices'])
        if 'default' in desc and 'help' in desc:
            desc['help'] += " (default: {})".format(desc['default'])
        if 'required' in desc and 'help' in desc:
            desc['help'] += " [REQUIRED]"
        if 'action' in desc and desc['action'] in ACTIONS:
            desc['action'] = getattr(
                sys.modules[__name__], desc['action'])
        arg = "--{}".format(arg)
        if 'positional' in desc:
            del desc['positional']
        p.add_argument(
            arg, **desc
        )
    add_default_arguments(p)
    p.set_defaults(cmd_class=cmd_class)



def get_global_config(include_cmds):
    args_dict = {}
    for cmd in include_cmds:
        spec = specfile_read(cmd.lower())
        args_from_yaml = read.from_yaml(spec)
        for key, value in args_from_yaml.items():
            if key in args_dict:
                # for full pipeline, ensure that
                # the outfile version replaces
                # any input files so argparse does
                # not assume they should be present
                if 'out' in args_dict[key]['type']:
                    continue
            args_dict[key] = value
    return args_dict




def add_default_arguments(parser):
    parser.add_argument(
        '--wd', '-wd', "--working_dir", type=str,
        default=os.getcwd(),
        help="Specify working directory to place output. "
             "(default: {})".format(os.getcwd())
    )
    parser.add_argument(
        '-c', "--config", type=outfile_type,
        help="Specify path to config file path, "
             "relative to working directory, "
             "not required if using default config. "
    )
    parser.add_argument(
        '-lf', "--log_file", type=outfile_type,
        default=os.path.join("./Logs", DEFAULT_LOG_FNAME),
        help="Set log file path, "
             "absolute or relative to working dir. "
             "(default: ./Logs/{})".format(DEFAULT_LOG_FNAME)
    )
    parser.add_argument(
        '-t', "--threads", type=positive_int,
        default=4,
        help="Number of worker threads to spawn. (default: 4)"
    )



def path_type(input_path):
    '''Path_type is for paths that should already exist.
    If path exists, returns absolute path, otherwise raise
    ArgumentTypeError'''
    if not os.path.isdir(input_path):
        raise argparse.ArgumentTypeError("Not a valid path")
    return os.path.abspath(input_path)

def flag_type(input_val):
    if input_val != "T" and input_val != "F":
        raise argparse.ArgumentTypeError(
            "Invalid flag, must be either T or F")
    if input_val == "T":
        return True
    if input_val == "F":
        return False

def project_dir_type(input_path):
    '''Creates a project directory if one does not exist and
    Throws PermissionError if there are no permissions to create
    directory. Returns absolute path to directory.'''
    input_path = os.path.abspath(input_path)
    if not os.path.isdir(input_path):
        try:
            os.mkdir(input_path)
        except PermissionError:
            error(
                "No permission to make directory: {}".format(input_path))
    if os.getcwd() != input_path:
        os.chdir(input_path)
    return input_path



def path_list_type(input_paths):
    '''Returns list of absolute paths using glob'''
    if not isinstance(input_paths, (list,)):
        input_paths = [input_paths]
    glob_list = []
    for path in input_paths:
        glob_list += glob(path)
    if not glob_list:
        raise argparse.ArgumentTypeError("Not a valid path")
    return list({os.path.abspath(path) for path in glob_list})

def arg_list_type(args):
    """ Returns a list of arguments """
    new_args = []
    for arg in args:
        arg = [a.replace("\'","").replace("\"", "") for a in arg.split()]
        new_args += arg
    return new_args



class Cmds(argparse.Action):
    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(Cmds, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        values = arg_list_type(values)
        setattr(namespace, self.dest, values)
    



class Glob(argparse.Action):
    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(Glob, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        values = path_list_type(values)
        setattr(namespace, self.dest, values)


def outpath_type(input_path):
    '''Outpath_type creates a directory if one does not exist.
    Throws PermissionError if there are no permissions to create
    directory. If path already exists and it is not empty, a warning
    is issued. Returns absolute path to directory.'''
    input_path = os.path.abspath(input_path)
    try:
        os.mkdir(input_path)
        logger.info("Creating directory: {}".format(input_path))
    except PermissionError:
        error(
            "No permission to make directory: {}".format(input_path))
    except OSError:
        logger.info(
            "Directory already exists: {}. ".format(input_path))
        if os.listdir(input_path):
            warn("Files in {} may be overwritten!".format(input_path))
    return input_path


def file_type(input_file):
    '''File type is for input files that
    should already exist.
    Check if input_file is a file.
    If yes, return absolute path to
    file else throw ArgumentTypeError
    exception'''

    input_file = os.path.abspath(input_file)
    if not os.path.isfile(input_file):
        raise argparse.ArgumentTypeError(
            "Not a valid file path: {}".format(input_file))
    return input_file


def outfile_type(input_file):
    '''Checks that path to file exists,
    if it doesn't, the path is created,
    returns abs path to file. PermissionError
    exception when there is no permission to 
    create directory'''
    input_file = os.path.abspath(input_file)
    path = os.path.dirname(input_file)
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except PermissionError:
            error(
                "No permission to create file: {}".format(input_file)
            )
    return input_file

def write_handle_type(input_file):
    '''Returns a handle to an outfile'''
    return open(outfile_type(input_file), 'w')


def read_handle_type(input_file):
    '''Returns a handle to an infile'''
    return open(file_type(input_file), 'r') 


def nonneg_int(input_val):
    ''' Make a non negative int type for argparse'''
    try:
        input_val = int(input_val)
        if input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Negative integer not valid")
    return input_val


def proportion(input_val):
    """make proportion type for argparse"""
    try:
        input_val = float(input_val)
        if input_val > 1 or input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a proportion")
    return input_val


def positive_int(input_val):
    '''Make a positive int type for argparse'''
    try:
        input_val = int(input_val)
        if input_val <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a positive integer")
    return input_val


TYPES = {
    'int': int,
    'str': str,
    'float': float,
    'bool': bool,
    'list': ast.literal_eval,
    'file_type': file_type,
    'outfile_type': outfile_type,
    'path_type': path_type,
    'outpath_type': outpath_type,
    'positive_int': positive_int,
    'nonneg_int': nonneg_int,
    'proportion': proportion,
    'project_dir_type': project_dir_type,
    'write_handle_type': write_handle_type,
    'read_handle_type': read_handle_type,
    'path_list_type': path_list_type,
    'flag_type': flag_type,
    'arg_list_type': arg_list_type
    }
ACTIONS = {'Glob': Glob, 'Cmds': Cmds}
