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
    specfile_path,
    specfile_read,
    parse_config_sections,
    get_global_config,
    get_missing_sections
)

from mtsv.utils import(
    error,
    warn,
    config_logging,
    set_log_file
)

from mtsv import (
    DEFAULT_CFG_FNAME,
    DEFAULT_LOG_FNAME)


COMMANDS = {"analyze": Analyze,
            "binning": Binning,
            "readprep": Readprep,
            "summary": Summary,
            "extract": Extract,
            "pipeline": Pipeline,
            }

def make_sub_parser(subparser, config, cmd, cmd_class):
    global_defaults = get_global_config(cmd_class.config_section)
    help_str = global_defaults["_meta_{}".format(cmd)]["help"]
    p = subparser.add_parser(cmd, help=help_str)
    for arg, desc in global_defaults.items():
        if "_meta" in arg:
            continue
        if 'type' in desc:
            desc['type'] = TYPES[desc['type']]
        if 'default' in desc and 'help' in desc:
            desc['help'] += " (default: {})".format(desc['default'])
            del desc['default']
        arg = "--{}".format(arg)
        if 'positional' in desc:
            del desc['positional']
        p.add_argument(
            arg, **desc
        )
    config_args = parse_config_sections(config, cmd_class.config_section)
    p.set_defaults(cmd_class=cmd_class, **config_args)


def add_cfg_to_args(argv, parser):
    '''treat config arguments as command line
    arguments to catch argparse errors'''
    args = parser.parse_args()
    d = vars(args)
    # change workingdir to abspath
    d['working_dir'] = os.getcwd()
    if "init" in argv:
        return args, []
    config_args = parse_config_sections(
        args.config, args.cmd_class.config_section)
    for k, v in config_args.items():
        fmt_k = "--{}".format(k)
        if fmt_k not in argv and v != None:
            argv += [fmt_k, v]
    missing = get_missing_sections(args.config)
    return parser.parse_args(argv[1:]), missing

    

def setup_and_run(argv, parser):
    """Setup and run a command."""
    args, missing = add_cfg_to_args(argv, parser)
    logger = logging.getLogger(__name__)
    args.log_file = set_log_file(
        args.log_file,
        args.cmd_class.__name__,
        args.timestamp)
    config_logging(args.log_file, args.log)
    if missing:
        warn(
            "Section(s) missing in config file, "
            "using defaults: {}".format(", ".join(missing)))
    params = Parameters(args)
    cmd = args.cmd_class(params)
    print(type(cmd))
    logger.info("Starting {}".format(cmd))

    cmd.run()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    dirparser = argparse.ArgumentParser(add_help=False)

    dirparser.add_argument(
        '-wd', "--working_dir", type=str,
        default=os.getcwd(),
        help="Specify working directory to place output"
    )

    dir_args, _ = dirparser.parse_known_args()
    working_dir = TYPES['project_dir_type'](dir_args.working_dir)
    os.chdir(working_dir)

    cfgparser = argparse.ArgumentParser(
        add_help=False
    )

    
    cfgparser.add_argument(
        '-c', "--config", type=TYPES['read_handle_type'],
        default=DEFAULT_CFG_FNAME,
        help="Specify path to config file, "
             "not required if using default config"
    )

    pre_args, _ = cfgparser.parse_known_args()

 
    parser = argparse.ArgumentParser(
        prog="MTSV",
        description="Metagenomic analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[dirparser, cfgparser]
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
    parser.add_argument(
        '-t', "--threads", type=TYPES['positive_int'],
        default=4,
        help="Number of worker threads to spawn."
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
            subparsers, pre_args.config,
            command, cmd_class)
    
    # Return help if no command is passed
    if len(argv)==1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    try:
        setup_and_run(argv, parser)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)


if __name__ == "__main__":
    main()
