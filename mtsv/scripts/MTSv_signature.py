import logging
import argparse
from ete3 import NCBITaxa
from collections import defaultdict
import numpy as np
from functools import partial
from multiprocessing import Pool, Process
from mtsv.utils import config_logging, get_ete_ncbi
from mtsv.parsing import (
    parse_output_row_bytes, 
    file_type, outfile_type,
    positive_int)



def separate_multiples_from_singletons(input_file):
    """The format of the input file should be:
            R1_1_0_1:123
            R2_0_10_2:123,222
            ...
        where R1_1_0_1 and R2_0_10_2 are query IDs and the numbers 
        following the colon are taxids to which the query aligned.
        singletons are queries that have only a single taxid and 
        are automatically assumed to be signature hits. Multiple
        taxids need to be processed to see if they are signature 
        at the provided roll up rank. The different combinations 
        of taxids are combined into a set (dictionary keys) to
        avoid redundant roll up calculations. The function returns
        a list of singleton records (named tuples) and a dictionary
        of multiples where the keys are the taxids for the roll up 
        and the values are a list of query records.
    """
    singleton_records = []
    multiples_records = defaultdict(list)
    with open(input_file, 'rb') as infile:
        append_single = singleton_records.append
        for line in infile:
            record = parse_output_row_bytes(line)
            # already signature append to signature list
            if len(record.taxa) == 1:
                append_single(record)
            else:
                multiples_records[tuple(record.taxa)].append(record)
    return singleton_records, multiples_records


def get_common_ancestor(taxids, rank):
    """The lineage for each taxid in the list of hits is determind and
    the intersection of each of these lineages is calculated to determine
    the common ancestors. This process is stopped early as soon as there
    are no common ancestors at the provided rank, indicating a non-signature
    hit. Otherwise the taxid of the common ancestor of the given rank 
    will be returned.
    """
    get_lineage = NCBI.get_lineage
    get_rank = NCBI.get_rank
    common_ancestor = np.array(get_lineage(taxids[0]), dtype=int)
    for taxid in taxids[1:]:
        common_ancestor = np.intersect1d(
            common_ancestor, get_lineage(taxid))
        if rank not in set(get_rank(common_ancestor).values()):
            return None
    return [key for key, value 
                in get_rank(common_ancestor).items()
                if value == rank][0]

def get_byte_stream(name, taxid):
    return name + b":" + bytes(str(taxid), 'ascii') + b"\n"


def write_signature(data, outfile):
    LOGGER.info("Writing singletons to file")
    with open(outfile, 'wb') as out:
        for query in data:
            out.write(
                get_byte_stream(
                    query.read_name, query.taxa[0]))


def write_signature_rollups(multiples, common_ancestors, outfile):
    LOGGER.info("Writing signature hits from roll up")
    with open(outfile, 'ab') as out:
        for ancestor, queries in zip(
                common_ancestors, multiples.values()):
            if ancestor is None:
                continue
            else:
                for query in queries:
                    out.write(
                        get_byte_stream(query.read_name, ancestor))
                

def signature(clps_file, rank, outfile, threads):
    LOGGER.info("Reading from file: {}".format(clps_file))
    singletons, multiples = separate_multiples_from_singletons(clps_file)
    LOGGER.info("Found {} queries with a single hit.".format(len(singletons)))
    proc = Process(target=write_signature, args=(singletons, outfile,))
    proc.start()
    LOGGER.info("Roll up level set to {}.".format(rank))
    if rank != 'species':
        LOGGER.info("Rolling up queries with multiple hits.")
        p = Pool(threads)
        common_ancestor_partial = partial(get_common_ancestor, rank=rank)
        common_ancestors = p.imap(
            common_ancestor_partial, multiples.keys(), chunksize=5000)
        p.close()
        p.join()
        proc.join()
        write_signature_rollups(multiples, common_ancestors, outfile)
    else:
        proc.join()
    LOGGER.info("MTSv Signature Completed")


if __name__ == "__main__":
    try:
        NCBI = get_ete_ncbi(snakemake.params[0])
        config_logging(snakemake.log[0], "INFO")
        LOGGER = logging.getLogger(__name__)
        signature(
            snakemake.input[0],
            snakemake.params[1],
            snakemake.output[0],
            snakemake.threads)
    except NameError:
        PARSER = argparse.ArgumentParser(
            prog="MTSv Signature",
            description="Find hits that are signature for a single taxon.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        PARSER.add_argument(
            "clp", metavar="COLLAPSE_FILE", type=file_type,
            help="Path to merged file, the output from mtsv-binning."
        )
        PARSER.add_argument(
            "sig", metavar="SIGNATURE_FILE", type=outfile_type,
            help="Path to output file."
        )
        PARSER.add_argument(
            '-r', "--rank", type=str, default="species",
            choices=["species", "genus", "family"],
            help="Consider hits to be signature up to a certain "
                 "taxanomic level. Ex: If genus is selected, a hit with "
                 "multiple species from the same genus will be considered to "
                 "be a signature hit for that genus, the same hit would not "
                 "be considered to be signature if the species level was "
                 "selected." 
        )
        PARSER.add_argument(
            "--taxdump", type=file_type, default=None,
            help="Alternative path to taxdump. "
                "Default is home directory where ete3 "
                "automatically downloads the file."
        )
        PARSER.add_argument(
            '-t', "--threads", type=positive_int, default=1,
            help="Number of threads for multiprocessing."
        )
        PARSER.add_argument(
            '-l', "--log", type=outfile_type, default="./signature.log",
            help="Path of log file."
        )

        ARGS = PARSER.parse_args()

        NCBI = NCBITaxa(taxdump_file=ARGS.taxdump)
        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)

        signature(
            ARGS.clp,
            ARGS.rank,
            ARGS.sig,
            ARGS.threads
        )
else:
    NCBI = NCBITaxa()
