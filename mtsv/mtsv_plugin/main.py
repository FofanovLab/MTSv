import argparse
import sys
from mtsv.commands import WGFast
from mtsv.parsing import make_sub_parser

COMMANDS = {
    "wgfast": WGFast
}


def main(argv=None):
    if argv is None:
        argv = sys.argv

    cfgparser = argparse.ArgumentParser(
        add_help=False
    )

    cfgparser.add_argument(
        '-c', "--config", type=str,
        help="Specify path to config file, "
             "not required if using default config."
    )
    
    pre_args, _ = cfgparser.parse_known_args()


    parser = argparse.ArgumentParser(
        prog="mtsv-plugin",
        description="Plugins and extensions to MTSv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[cfgparser]
    )
    
    
    subparsers = parser.add_subparsers(
        title="commands", metavar="COMMANDS",
        help="Plugin Commands"
    )

    for command, cmd_class in COMMANDS.items():
        make_sub_parser(
            subparsers, command, cmd_class
        )
    if len(argv)==1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    try:
        setup_and_run(parser)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)

