import textwrap
import sys
import click
import logging

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
        logger.error(msg)
        sys.exit(1)


def warn(msg):
    '''Prints a warning message to stderr.'''
    text = quote(msg, quote="!> ", nl=False)
    click.secho(text, err=True, fg='red')
    logger.warning(msg)


def quote(text, width=72, quote="", nl=True):
    if not text:
        return ""
    out = ""
    for line in text.split("\n"):
        sublines = textwrap.wrap(line, width=width, replace_whitespace=False)
        sublines = [quote + l for l in sublines]  # if l.strip()]
        out += "\n".join(sublines) + "\n"
    return out


def config_logging(handle, level):
    """Configure logging

    Arguments:
        log_file_name (handle): log file handle
        level (str): Logging level
    """
    logging.basicConfig(
        stream=handle,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=getattr(logging, level),
        filemode='a',
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')


def set_log_file(log_file, command_name, timestamp):
    log_file = log_file.format(
        COMMAND=str(command_name).lower(),
        TIMESTAMP=timestamp)
    return open(log_file, 'w')
