import argparse
import configparser
import os
import logging
import ast
from pkg_resources import resource_stream, resource_filename
from argutils import read, export
from utils import(
    error,
    warn
)

logger = logging.getLogger(__name__)


def create_config_file(commands, config_file):
    cfg_file_str = ""
    for cmd in commands:
        spec = specfile_read(cmd)
        opts = read.from_yaml(spec)
        meta = opts['_meta_{}'.format(cmd.lower())]
        del opts['_meta_{}'.format(cmd.lower())]
        opts['_meta'] = meta
        cfg_file_str += export.to_config(
            cmd, opts) + "\n"
    config_file.write(cfg_file_str)
    config_file.close()

def get_global_config(include_cmds):
    args_dict = {}
    for cmd in include_cmds:
        spec = specfile_read(cmd.lower())
        args_dict.update(read.from_yaml(spec))
    return args_dict

def get_global_defaults(cmd_class):
    return {k: {'default': v['default'],
        'type': v['type']} for 
        k, v in get_global_config(
            cmd_class.config_section).items()
        if 'default' in v and 'type' in v}


def set_types(params, types):
    for key, value in params.items():
        if key in types:
            params[key] = TYPES[types[key]](value)
    return params

def parse_config_sections(config_file, sections):
    # pass for init incase a config file already
    # exists
    if config_file.mode == 'w':
        return {}
    config = configparser.ConfigParser()
    config.read(config_file.name)
    config_for_sections = {}
    try:
        for section in sections:
            try:
                for cmd, val in config.items(section):
                    print(cmd, val)
                    config_for_sections[cmd] = val
            except configparser.NoSectionError:
                warn(
                    "{} section missing in config file, "
                    "using defaults".format(section))
    except configparser.ParsingError:
        error(
            "Cannot parse config file: {}".format(
                config_file.name))
    return config_for_sections

def path_type(input_path):
    '''Path_type is for paths that should already exist.
    If path exists, returns absolute path, otherwise raise
    ArgumentTypeError'''
    if not os.path.isdir(input_path):
        raise argparse.ArgumentTypeError("Not a valid path")
    return os.path.abspath(input_path)


def project_dir_type(input_path):
    '''Creates a project directory if one does not exist and
    changes the current working directory to this directory.
    Throws PermissionError if there are no permissions to create
    directory. If path already exists and it is not empty, a warning
    is issued. Returns absolute path to directory.'''
    print(input_path)
    try:
        os.mkdir(input_path)
        logger.info("Creating Working Directory: {}".format(input_path))
    except PermissionError:
        error(
            "No permission to make directory: {}".format(input_path))
    except OSError:
        logger.info(
            "Directory already exists: {}. ".format(input_path))
        if os.listdir(input_path):
            warn("Files in {} may be overwritten!".format(input_path))
    path = os.path.abspath(input_path)
    os.chdir(path)
    return path


def outpath_type(input_path):
    '''Outpath_type creates a directory if one does not exist.
    Throws PermissionError if there are no permissions to create
    directory. If path already exists and it is not empty, a warning
    is issued. Returns absolute path to directory.'''
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
    return os.path.abspath(input_path)


def file_type(input_file):
    '''File type is for input files that
    should already exist.
    Check if input_file is a file.
    If yes, return absolute path to
    file else throw ArgumentTypeError
    exception'''
    if not os.path.isfile(input_file):
        raise argparse.ArgumentTypeError("Not a valid file path")
    return os.path.abspath(input_file)


def outfile_type(input_file):
    '''Checks that path to file exists,
    if it doesn't, the path is created,
    returns abs path to file. PermissionError
    exception when there is no permission to 
    create directory'''
    path = os.path.dirname(os.path.abspath(input_file))
    print(os.path.abspath(input_file))
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except PermissionError:
            error(
                "No permission to create file: {}".format(input_file)
            )
    return os.path.join(path, input_file)

def write_handle_type(input_file):
    '''Returns a handle to an outfile'''
    return open(outfile_type(input_file), 'w')


def read_handle_type(input_file):
    '''Returns a handle to an infile'''
    
    return open(file_type(input_file), 'r') 


def positive_int(input_val):
    '''Make a positive int type for argparse'''
    try:
        input_val = int(input_val)
        if input_val <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Not a positive integer")
    return input_val


def nonneg_int(input_val):
    ''' Make a non negative int type for argparse'''
    try:
        input_val = int(input_val)
        if input_val < 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Negative integer not valid")
    return input_val



def specfile_read(name):
    """Return the specfile stream for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name + '.yml')
    return resource_stream(__name__, fp)


def specfile_path(name):
    """Return the specfile path for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name + '.yml')
    return resource_filename(__name__, fp)


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
    'project_dir_type': project_dir_type,
    'write_handle_type': write_handle_type,
    'read_handle_type': read_handle_type
}