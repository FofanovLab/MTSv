import textwrap
try:
    from ConfigParser import NoSectionError
except:
    from configparser import NoSectionError
from warnings import warn

META_KEY = "_meta"
DESC_KEY = "help"
EXCLUDE_FLAG = "_exclude"
FILE_W = "File-w"
FILE_R = "File-r"

def set_parser_defaults(parser, config):
    """Sets the defaults for an ArgumentParser from a ConfigParser object.

    :param parser: the ArgumentParser to set
    :param config: the ConfigParser to read the defaults from. Should have a 
    section matching the ArgumentParser's `prog` field.
    :returns: the updated ArgumentParser
    """

    try:
        defaults = dict(config.items(parser.prog))
        parser.set_defaults(**defaults)
        return parser
    except NoSectionError:
        warn("Section [{0}] not found in config file".format(parser.prog))
        return parser

def format_comment(text, width=72, quote="# "):
    """Line-wraps and pads text to write as a comment.

    :param text: the text to format
    :param width: number of characters to wrap text to
    :param quote: the string to prefix each line (include a space!)
    """

    if not text: 
        return ""
    lines = textwrap.wrap(text, width=width)
    lines = [quote + line.strip() for line in lines if line.strip()]
    return "\n".join(lines) + "\n"