import argparse
import logging
import sys
import datetime
import os
import configparser

from mtsv.argutils import (read, export)
from mtsv.parameters import Parameters
from mtsv.commands import (
    Init,
    Readprep,
    Binning,
    Summary,
    Analyze,
    Extract,
    Pipeline
)

from mtsv.parsing import (
    TYPES,
    make_sub_parser,
    parse_config_sections,
    get_missing_sections,
    add_default_arguments)

from mtsv.utils import(
    error,
    warn,
    config_logging,
    set_log_file
)

from mtsv import (
    DEFAULT_LOG_FNAME,
    DEFAULT_CFG_FNAME)

COMMANDS = {
    "analyze": Analyze,
    "binning": Binning,
    "readprep": Readprep,
    "summary": Summary,
    "extract": Extract,
    "pipeline": Pipeline
}

def add_cfg_to_args(argv, parser):
    '''treat config arguments as command line
    arguments to catch argparse errors'''
    config = get_config_from_argv(argv)
    config_args = parse_config_sections(
        config,
        get_command_from_argv(argv).config_section)
    for k, v in config_args.items():
        fmt_k = "--{}".format(k)
        if fmt_k not in argv and v != None:
            argv += [fmt_k, v]
    missing = get_missing_sections(config)
    args, snake_args = parser.parse_known_args(argv[1:])
    return args, snake_args, missing



def get_command_from_argv(argv):
    return COMMANDS[argv[1]]

def get_config_from_argv(argv):
    index = -1
    opts = ['-c', '--config']
    for opt in opts:
        if opt in argv:
            index = argv.index(opt)
    if index != -1:
        return argv[index + 1]
        


def change_wkdir(argv):
    index = -1
    opts = ['--working_dir', '-wd']
    for opt in opts:
        if opt in argv:
            index = argv.index(opt)
    if index != -1:
        argv[index + 1] = TYPES['project_dir_type'](argv[index + 1])

def setup_and_run(argv, parser):
    """Setup and run a command."""
    change_wkdir(argv)
    if argv[1] != "init":
        if '--config' in argv or '-c' in argv:
            args, snake_args, missing = add_cfg_to_args(argv, parser)
            if missing:
                warn(
                "Section(s) missing in config file, "
                "using defaults: {}".format(", ".join(missing)))
        else:
            args, snake_args = parser.parse_known_args()
        args.log_file = set_log_file(
            args.log_file,
            args.cmd_class.__name__,
            args.timestamp)
        config_logging(args.log_file, args.log)
    else:
        args, snake_args = parser.parse_args(), []

    params = Parameters(args, snake_args)

    cmd = args.cmd_class(params)

    cmd.run()



def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        prog="mtsv",
        description="Metagenomic analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.set_defaults(timestamp=datetime.datetime.now().strftime(
        '%Y-%m-%d_%H-%M-%S'))

    subparsers = parser.add_subparsers(
        title="commands", metavar="COMMAND",
        help="Pipeline Commands")

    parser_init = subparsers.add_parser(
        'init',
        help="initializes a directory with a pre-filled parameters file"
    )
    parser_init.add_argument(
        "-c", "--config", type=TYPES['write_handle_type'],
        default=DEFAULT_CFG_FNAME,
        help="Specify path to write config file, "
        "not required if using default config"
    )
    parser.add_argument(
        '-wd', "--working_dir", type=str,
        default=os.getcwd(),
        help="Specify working directory to place output. "
        "(default: {})".format(os.getcwd())
    )
    parser_init.set_defaults(cmd_class=Init)

    for command, cmd_class in COMMANDS.items():
        make_sub_parser(
            subparsers, command, cmd_class)

    # Return help if no command is passed
    if len(argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    try:
        setup_and_run(argv, parser)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)


if __name__ == "__main__":
    main()
