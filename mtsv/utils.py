import textwrap
import sys
import click
import logging
import os
import json
import hashlib
from ete3 import NCBITaxa
from pkg_resources import resource_stream, resource_filename


LOGGER = logging.getLogger(__name__)


class MTSvError(Exception):
    pass

def error(msg, exception=False, wrap=True, trunc=100):
    """Prints an error message to logger and stderr and exits."""
    if exception:
        raise MTSvError(quote(msg, trunc=trunc, width=75) if wrap else msg)
    else:
        msg = quote(msg, trunc=trunc, width=75) if wrap else msg
        click.secho(msg, fg='red', err=True, bold=True)
        sys.exit(1)


def warn(msg):
    """Prints a warning message to stderr."""
    text = quote(msg, quote="!> ", nl=False)
    click.secho(text, err=True, fg='magenta', bold=True)

def info(msg, color="green"):
    """
    Prints info to stdout
    """
    text = quote(msg, nl=False)
    click.secho(text, fg=color, err=True)


def quote(text, width=80, quote="", nl=True, trunc=100):
    if not text:
        return ""
    out = ""
    for line in text.split("\n"):
        sublines = textwrap.wrap(line, width=width, replace_whitespace=False)
        sublines = [quote + l for l in sublines]  # if l.strip()]
        out += "\n".join(sublines) + "\n"
    if trunc is not None:
        out_list = out.split("\n")
        dots = "\n..." if len(out_list) > trunc else ""
        out = "\n".join(out_list[:trunc]) + dots 
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
        filemode='a',
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')



def specfile_read(name):
    """Return the specfile stream for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name.lower() + '.yml')
    return resource_stream(__name__, fp)


def specfile_path(name):
    """Return the specfile path for a given command name."""
    fp = os.path.join('commands', 'cmd_specs', name.lower() + '.yml')
    return resource_filename(__name__, fp)


def script_path(script_name):
    """ Return the script path for given script name."""
    fp = os.path.join('scripts', script_name)
    return resource_filename('mtsv', fp)

def template_path(template_name):
    """
    Return path to template.
    """
    fp = os.path.join('templates', template_name)
    return resource_filename(__name__, fp)

def snake_path(file_name):
    """ Return the script path for given script name."""
    fp = os.path.join('commands', 'snakefiles', file_name)
    return resource_filename('mtsv', fp)


def get_user_mtsv_dir():
    user = os.environ.get('HOME', '/')
    mtsv_path = os.path.join(user, ".mtsv")
    os.makedirs(mtsv_path, exist_ok=True)
    return mtsv_path

def data_path():
    """ Return expected values database path """
    mtsv_path = get_user_mtsv_dir()
    data_file = os.path.join(mtsv_path, "expected_datastore.hdf5")
    return data_file

def ete_database_data():
    """ Return path to ete3 database json """
    mtsv_path = get_user_mtsv_dir()
    fp = os.path.join(mtsv_path, "ete_databases.json")
    if not os.path.isfile(fp):
        with open(fp, 'w') as outfile:
            outfile.write("{}")
    return fp

def get_ete_ncbi(taxdump, quiet=False):
    #sqlite3.OperationalError
    ete_json = json.loads(open(ete_database_data(), 'r').read())
    mod_time = os.path.getmtime(taxdump)
    user = os.environ.get('HOME', '/')
    if user in ete_json and taxdump in ete_json[user]:
        entry = ete_json[user][taxdump]
        db_path = entry['db_path']
        if entry['modtime'] == mod_time:
            if os.path.isfile(db_path):
                if not quiet:
                    LOGGER.info("Ete taxdump database already exists")
                ncbi = NCBITaxa(dbfile=db_path)
                return ncbi
            else:
                if not quiet:
                    LOGGER.info(
                        "Ete taxdump database has been deleted, rebuilding")
        else:
            if not quiet:
                LOGGER.info(
                    "Old version of ete taxdump database exists, updating")
            ete_json[user][taxdump]['modtime'] == mod_time
    else:
        new_name = "mtsv_{}_taxa.sqlite".format(
            hashlib.md5(bytes(taxdump, 'utf8')).hexdigest())
        db_path = os.path.join(user, '.etetoolkit')
        if not os.path.isdir(db_path):
            os.makedirs(db_path)
        db_path = os.path.join(db_path, new_name)
        ete_json[user]= {taxdump: {'modtime': mod_time, 'db_path': db_path}}
        if not quiet:
            LOGGER.info(
                "New ete taxdump database created for taxdump {}".format(
                    taxdump))
    ncbi = NCBITaxa(dbfile=db_path, taxdump_file=taxdump)
    with open(ete_database_data(), 'w') as out:
        out.write(json.dumps(ete_json))
    return ncbi


def get_item_from_artifact(artifact_path, item_key):
    artifact = json.loads(open(artifact_path, 'r').read())
    return artifact[item_key]

    
    



