import argparse
import configparser
import os
import sys
import ast
import logging
import numpy as np
from contextlib import suppress
from glob import glob
from collections import namedtuple


from mtsv.utils import(error, warn, specfile_read)
from mtsv.argutils import (read, export)
from mtsv import (DEFAULT_LOG_FNAME, DEFAULT_CFG_FNAME)

logger = logging.getLogger(__name__)
SECTIONS = ["READPREP", "BINNING", "SUMMARY", "ANALYZE", "EXTRACT", "WGFAST"]

split = str.split
strip = str.strip
rsplit = str.rsplit

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

def create_config_file(config_file):
    cfg_file_str = ""
    for cmd in SECTIONS:
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


def format_cml_params(include_cmd, args, ignore, include):
    command_list = []
    for arg, desc in get_global_config([include_cmd]).items():
        if '_meta' in arg:
            continue
        if 'positional' in desc:
            continue
        if arg in ignore:
            continue
        if arg in args and args[arg] is not None:
            try:
                val = args[arg].name
            except AttributeError:
                val = str(args[arg]) 
            command_list += ["--{}".format(arg), val]
    for incl in include:
        command_list.append(incl)
    return " ".join(command_list)


def add_default_arguments(parser):
    parser.add_argument(
        '-wd', "--working_dir", type=str,
        default=os.getcwd(),
        help="Specify working directory to place output. "
             "(default: {})".format(os.getcwd())
    )
    parser.add_argument(
        '-c', "--config", type=outfile_type,
        help="Specify path to config file path, "
             "relative to working directory, "
             "not required if using default config. "
             "(default: {})".format(DEFAULT_CFG_FNAME)
    )
    parser.add_argument(
        '-lf', "--log_file", type=outfile_type,
        default=DEFAULT_LOG_FNAME,
        help="Set log file path, "
             "absolute or relative to working dir. "
             "(default: ./{})".format(DEFAULT_LOG_FNAME)
    )
    parser.add_argument(
        '-t', "--threads", type=positive_int,
        default=4,
        help="Number of worker threads to spawn. (default: 4)"
    )

def get_missing_sections(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    return [s for s in SECTIONS if s not in config.sections()]

def parse_config_sections(config_file, sections):
    config = configparser.ConfigParser()
    config.read(config_file)
    config_for_sections = {}
    try:
        for section in sections:
            with suppress(configparser.NoSectionError):
                for cmd, val in config.items(section):
                    # Only add if there is a value
                    # avoids adding optional params
                    if val:
                        config_for_sections[cmd] = val                    
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

def flag_type(input_val):
    if not set(input_val).intersection(set("TF")):
        raise argparse.ArgumentTypeError(
            "Invalid flag, must be either T or F")
    return True if input_val == "T" else False

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
    glob_list = []
    for path in input_paths:
        glob_list += glob(path)
    if not glob_list:
        raise argparse.ArgumentTypeError("Not a valid path")
    return list({os.path.abspath(path) for path in glob_list})

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


Record = namedtuple('Record', ['read_id', 'counts', 'taxa', 'read_name'])

def parse_output_row(row):
    read_name, taxa = split(row, ":")
    taxa = np.array([tax for tax in split(taxa, ",")], dtype=int)
    counts = np.array([c for c in split(read_name, "_")[1:]], dtype=int)
    read_id = split(read_name, "_")[0]
    return Record(
        read_id=read_id, counts=counts, taxa=taxa, read_name=read_name)



def parse_query_id(query_id):
    return np.array([int(q) for q in query_id.split("_")[1:]])


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
    'flag_type': flag_type
    }
ACTIONS = {'Glob': Glob}
