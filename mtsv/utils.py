import textwrap
import sys
import click
import logging
import os
from pkg_resources import resource_stream, resource_filename


logger = logging.getLogger(__name__)


class MTSVError(Exception):
    pass

def error(msg, exception=True, wrap=True):
    '''Prints an error message to logger and stderr and exits.'''
    if exception:
        raise MTSVError(msg)
    else:
        msg = quote(msg, width=75) if wrap else msg
        click.secho(msg, fg='red', err=True)
        sys.exit(1)


def warn(msg):
    '''Prints a warning message to stderr.'''
    text = quote(msg, quote="!> ", nl=False)
    click.secho(text, err=True, fg='red')


def quote(text, width=72, quote="", nl=True):
    if not text:
        return ""
    out = ""
    for line in text.split("\n"):
        sublines = textwrap.wrap(line, width=width, replace_whitespace=False)
        sublines = [quote + l for l in sublines]  # if l.strip()]
        out += "\n".join(sublines) + "\n"
    return out


def config_logging(file_name, level):
    """Configure logging

    Arguments:
        level (str): Logging level
    """
    logging.basicConfig(
        filename=file_name,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=getattr(logging, level),
        filemode='w',
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')


def set_log_file(log_file, command_name, timestamp):
    return log_file.format(
        COMMAND=str(command_name).lower(),
        TIMESTAMP=timestamp)


def specfile_read(name):
    """Return the specfile stream for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name.lower() + '.yml')
    return resource_stream(__name__, fp)


def specfile_path(name):
    """Return the specfile path for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name.lower() + '.yml')
    return resource_filename(__name__, fp)


def bin_path(cmd):
    """Return the binary path for a given command name."""
    fp = os.path.join('ext', cmd)
    return resource_filename('mtsv', fp)

def script_path(script_name):
    """ Return the script path for given script name."""
    fp = os.path.join('scripts', script_name)
    return resource_filename('mtsv', fp)
