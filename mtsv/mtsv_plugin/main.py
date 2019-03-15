import argparse
import sys
import datetime
from mtsv.commands import WGFast, Concoct
from mtsv.parameters import Parameters

from mtsv.parsing import (
    make_sub_parser,
    TYPES,
    parse_config_sections,
    get_missing_sections
)

from mtsv.utils import(
    error,
    warn,
    set_log_file
)

COMMANDS = {
    "wgfast": WGFast,
    "concoct": Concoct
}

def add_cfg_to_args(argv, parser):
    '''treat config arguments as command line
    arguments to catch argparse errors'''
    config = get_config_from_argv(argv)
    cmd_cfg_section = get_command_from_argv(argv).config_section
    config_args = parse_config_sections(
        config,
        cmd_cfg_section)
    for k, v in config_args.items():
        fmt_k = "--{}".format(k)
        if fmt_k not in argv and v != None:
            argv += [fmt_k] + v.split(" ")
    missing = set(cmd_cfg_section).intersection(
        set(get_missing_sections(config)))
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
    """Setup and run a command"""
    change_wkdir(argv)
    if '--config' in argv or '-c' in argv:
        args, snake_args, missing = add_cfg_to_args(argv, parser)
        if missing:
            warn("Section(s) missing in config file, "
                "using defaults: {}".format(", ".join(missing)))
    else:
        args, snake_args = parser.parse_known_args()
    args.log_file = set_log_file(
        args.log_file, 
        args.cmd_class.__name__,
        args.timestamp
    )
    params = Parameters(args, snake_args)
    cmd = args.cmd_class(params)
    cmd.run()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        prog="mtsv-plugin",
        description="Plugins and extensions to MTSv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.set_defaults(timestamp=datetime.datetime.now().strftime(
        '%Y-%m-%d_%H-%M-%S'))
    
    subparsers = parser.add_subparsers(
        title="commands", metavar="COMMANDS",
        help="Plugin Commands"
    )

    for command, cmd_class in COMMANDS.items():
        make_sub_parser(
            subparsers, command, cmd_class
        )
    
    # Return help if no command is passed
    if len(argv)==1:
        parser.print_help(sys.stdout)
        sys.exit(0)
    try:
        setup_and_run(argv, parser)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)


if __name__ == "__main__":
    main()

