import textwrap
import sys
import click
import logging
import os
import json
import hashlib
import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO
from pkg_resources import resource_stream, resource_filename


LOGGER = logging.getLogger(__name__)


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
        filemode='a',
        format='%(asctime)s %(levelname)s: [%(name)s] %(message)s')


def set_log_file(log_file, command_name, timestamp):
    log_file = log_file.format(
        COMMAND=str(command_name).lower(),
        TIMESTAMP=timestamp)
    # clear out previous log
    if os.path.isfile(log_file):
        open(log_file, 'w').close()
    return log_file


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


def snake_path(rule_name):
    """ Return the script path for given script name."""
    fp = os.path.join('commands', 'snakefiles', rule_name)
    return resource_filename('mtsv', fp)

def data_path():
    """ Return expected values database path """
    user = os.environ.get('HOME', '/')
    data_file = os.path.join(user, ".mtsv/expected.csv")
    try:
        open(data_file, 'r')
    except FileNotFoundError:
        os.makedirs(os.path.join(user, ".mtsv"))
        init_precalculated_df(data_file)
    return data_file

def ete_database_data():
    """ Return path to ete3 database json """
    fp = os.path.join('data', 'ete_databases.json')
    return resource_filename('mtsv', fp)

def get_ete_ncbi(taxdump):
    #sqlite3.OperationalError
    ete_json = json.loads(open(ete_database_data(), 'r').read())
    mod_time = os.path.getmtime(taxdump)
    user = os.environ.get('HOME', '/')
    if user in ete_json and taxdump in ete_json[user]:
        entry = ete_json[user][taxdump]
        db_path = entry['db_path']
        if entry['modtime'] == mod_time:
            if os.path.isfile(db_path):
                LOGGER.info("Ete taxdump database already exists")
                ncbi = NCBITaxa(dbfile=db_path)
                return ncbi
            else:
                LOGGER.info(
                    "Ete taxdump database has been deleted, rebuilding")
        else:
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
        LOGGER.info(
                "New ete taxdump database created for taxdump {}".format(
                    taxdump))
    ncbi = NCBITaxa(dbfile=db_path, taxdump_file=taxdump)
    with open(ete_database_data(), 'w') as out:
        out.write(json.dumps(ete_json))
    return ncbi

    
def init_precalculated_df(data_file):
    with open(data_file, 'w') as out:
        out.write("Database,Kmer_Size,Edits,Seed_Size,"\
                  "Seed_Gap,Min_Seeds,Taxid,Total_Hits,Sig_Hits,Ratio\n")


def get_precalculated_df():
    return pd.read_csv(
        data_path(),
        dtype={
            'Database': str,
            'Kmer_Size': int,
            'Edits': int,
            'Seed_Size': int,
            'Seed_Gap': int,
            'Min_Seeds': int,
            'Taxid': int,
            'Total_Hits': int,
            'Sig_Hits': int,
            'Ratio': float
        })

def write_to_precalculated_db(df):
    df.to_csv(data_path(), index=False)

def line_generator(file_name, n_lines):
    go = True
    with open(file_name, 'r') as infile:
        while go:
            lines = []
            for _ in range(n_lines):
                l = infile.readline()
                if l == "":
                    go = False
                    break
                lines.append(l)
            yield lines
        return 

def fasta_generator(file_name, n_records):
    with open(file_name, 'r') as handle:
        records = {}
        for record in SeqIO.parse(handle, 'fasta'):
            records[record.id] = record
            if len(records) == n_records:
                yield records
                records = {}
        if len(records):
            yield records
        return


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
               


#    while go:
#        while dat.count("\n") < n_lines:
#            if f.tell() < file_size:
#                dat += f.read(n_lines*500)
#            else:
#                go = False
#                break
#        lines = [l for l in dat.split("\n") if l != ""]
#        dat = "\n".join(lines[n_lines:])
#        print("GEN", lines[:n_lines][-1])
#        yield lines[:n_lines]
#    return

def track_file_params(
    file_type, file_path, params):
    track_file = os.path.join(
        os.path.dirname(file_path), ".params")
    track_dict = {file_type:
                    {file_path: params}}
    if os.path.isfile(track_file):
        json_data = open(track_file).read()
        record = json.loads(json_data)
        if file_type in record:
            record[file_type][file_path] = track_dict[file_type][file_path]
        else:
            record[file_type] = track_dict[file_type]
    else:
        record = track_dict
    record = json.dumps(record)
    with open(track_file, 'w') as out:
        out.write(record)

def get_database_params(filepath, value):
    try:
        return json.loads(open(filepath, 'r').read())[value]
    except IOError:
        error("{} file is not present. Avoid moving datafiles because "
              "the directory includes required metadata".format(filepath))
    except ValueError:
        error("Cound not parse json: {}".format(filepath))
    



