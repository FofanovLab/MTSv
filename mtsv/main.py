from __future__ import absolute_import
import argparse
import logging
import sys
import datetime
import os

from mtsv import (
    DEFAULT_CFG_FNAME,
    DEFAULT_LOG_FNAME)

from commands import (
    Init,
    Readprep,
    Binning,
    Summary,
    Analyze,
    Extract,
    Pipeline
)

from parsing import (
    TYPES,
    specfile_path,
    specfile_read,
    get_global_config)

from utils import(
    error,
    warn,
    config_logging,
    set_log_file
)

# command_name: (cmd_class, include)
COMMANDS = {"analyze": Analyze,
            "binning": Binning,
            "readprep": Readprep,
            "summary": Summary,
            "extract": Extract,
            "pipeline": Pipeline,
            }

def make_sub_parser(subparser, cmd, args_dict, cmd_class):
    help_str = args_dict["_meta_{}".format(cmd)]["help"]
    p = subparser.add_parser(cmd, help=help_str)
    for arg, desc in args_dict.items():
        if "_meta" in arg:
            continue
        if 'type' in desc:
            desc['type'] = TYPES[desc['type']]
        if 'default' in desc and 'help' in desc:
            desc['help'] += " (default: {})".format(desc['default'])
            del desc['default']
        if 'positional' not in desc:
            arg = "--{}".format(arg)
        else:
            del desc['positional']
        p.add_argument(
            arg, **desc
        )
    p.set_defaults(cmd_class=cmd_class)



def setup_and_run(args):
    """Setup and run a command."""
    cmd = args.cmd_class()
    logger = logging.getLogger(__name__)
    args.log_file = set_log_file(
        args.log_file, str(cmd), args.timestamp)
    config_logging(args.log_file, args.log)
    logger.info("Starting {}".format(cmd))
    cmd.set_parameters(**vars(args))
    cmd.run()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        prog="MTSV",
        description="Metagenomic analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-c", "--config", type=TYPES['file_type'], default=DEFAULT_CFG_FNAME,
        help="Specify path to config file, "
        "not required if using default config")
    parser.add_argument(
        '-wd', "--working_dir", type=TYPES['project_dir_type'], default=os.getcwd(),
        help="Specify working directory to place output"
    )
    parser.add_argument(
        '-f', "--force", action="store_true",
        help="Force rerun of steps"
    )
    parser.add_argument(
        '-cls', "--cluster_cfg", type=TYPES['file_type'],
        help="Cluster configuration file"
    )
    parser.add_argument(
        '-l', "--log", type=str, default="INFO",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help="Set logging level"
    )
    parser.add_argument(
        '-lf', "--log_file", type=str,
        default=DEFAULT_LOG_FNAME,
        help="Log file"
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
        "-c", "--config", type=TYPES['outfile_type'], default=DEFAULT_CFG_FNAME,
        help="Specify path to write config file, "
        "not required if using default config"
    )
    parser_init.set_defaults(cmd_class=Init)

            
    for command, cmd_class in COMMANDS.items():
        args_dict = get_global_config(cmd_class)
        make_sub_parser(subparsers, command, args_dict, cmd_class)
    
    # Return help if no command is passed
    if len(sys.argv)==1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    args = parser.parse_args()
    try:
        setup_and_run(args)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)


if __name__ == "__main__":
    main()
