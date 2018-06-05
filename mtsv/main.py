from __future__ import absolute_import
import argparse
import logging
import sys
import datetime
import os
import configparser

from mtsv.argutils import (read, export)
from mtsv.provenance import Parameters
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


def add_cfg_to_args(argv, args, parser):
    '''treat config arguments as command line
    arguments to catch argparse errors'''
    config_args = parse_config_sections(
        args.config,
        args.cmd_class.config_section)
    for k, v in config_args.items():
        fmt_k = "--{}".format(k)
        if fmt_k not in argv and v != None:
            argv += [fmt_k, v]
    missing = get_missing_sections(args.config)
    return parser.parse_args(argv[1:]), missing


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
    args = parser.parse_args()
    if args.cmd_class.__name__ != "Init":
        if args.config is not None:
            args, missing = add_cfg_to_args(argv, args, parser)
            if missing:
                warn(
                "Section(s) missing in config file, "
                "using defaults: {}".format(", ".join(missing)))
        args.log_file = set_log_file(
            args.log_file,
            args.cmd_class.__name__,
            args.timestamp)
        config_logging(args.log_file, args.log)


    params = Parameters(args)
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
