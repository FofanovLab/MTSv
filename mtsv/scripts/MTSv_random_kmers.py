import argparse
import logging
import time
import mmap
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from collections import namedtuple
from mtsv.utils import config_logging
from mtsv.mtsv_prep.MTSv_prune import deserialization, get_tree
from mtsv.parsing import file_type, outfile_type, positive_int

DNA = set(b"ACTG")

POS_RANGE = namedtuple('Ranges', ['start', 'stop'])

def get_eof(file_name):
    with open(file_name, 'r') as infile:
        infile.seek(0, 2)     # go to the file end.
        return infile.tell()


def get_random_positions(ranges, n_samples):
    lens = [r.stop - r.start for r in ranges]
    total = sum(lens)
    if n_samples*2 >= total:
        pos = np.array([])
        for r in ranges:
            pos = np.append(pos, np.arange(r.start, r.stop))
        return pos
    p = [l/total for l in lens]
    ns = np.random.multinomial(n_samples*2, p, size=1)[0]
    pos = np.array([], dtype=int)
    for n, _range in zip(ns, ranges):
        pos = np.append(pos, np.random.choice(np.arange(_range.start, _range.stop), replace=True, size=n))
    return np.unique(pos)


def get_outstring(kmers, sample, total_samples):
    count_arr = np.array(np.zeros(total_samples, dtype=int), dtype=bytes)
    count_arr[sample] = b"1"
    count_arr = b"_".join(count_arr)
    s = b""
    for kmer in kmers:
        s += b">R" + next(ID) + b"_" + count_arr + b"\n" + kmer + b"\n"
    return s
    

def get_kmers_from_taxon(
    taxon, fasta, tax_ids, positions, kmer_size, n_kmers):
    taxons = set(get_tree(tax_ids, taxon, positions))
    tx_positions = []
    LOGGER.info("Running taxon: {}".format(taxon))
    for tx in taxons:
        if tx:
            try:
                tx_positions += positions[tx.encode()]
            except KeyError:
                continue
    tx_positions.sort()
    with open(fasta, 'rb') as fasta_handle:
        ranges = []
        kmers = []
        for position in tx_positions:
            fasta_handle.seek(position)
            fasta_handle.readline()  # header
            r1 = fasta_handle.tell()
            r2 = r1
            line = fasta_handle.readline()
            while line and chr(line[0]) != ">":
                r2 = fasta_handle.tell()
                line = fasta_handle.readline()
            r2 = r2 - kmer_size
            if r2 > r1 and (r2 - r1) > kmer_size :
                ranges.append(POS_RANGE(r1,r2))
        if not ranges:
            LOGGER.info("No kmers found for taxid {}".format(taxon))
            return ""
        rand_positions = get_random_positions(ranges, n_kmers)
        for pos in rand_positions:
            fasta_handle.seek(pos)
            kmer = fasta_handle.read(
                    kmer_size * 2).replace(
                        b"\n", b"")[:kmer_size]
            if valid_kmer(kmer, kmer_size):
                kmers.append(kmer)
        kmers = np.unique(kmers)
        np.random.shuffle(kmers)
        kmers = kmers[:n_kmers]
        LOGGER.info(
            "Found {0} kmers for taxon {1}".format(
                len(kmers), taxon))
        return kmers


def valid_kmer(kmer, kmer_size):
    '''check if kmer is correct size (in case where)
    end of file is reached and check whether it
    it contains a newline char which indicates that
    the kmer is bridging two sequences '''
    if not set(kmer).issubset(DNA):
        return False
    if b'\n' in kmer:
        return False
    if len(kmer) < kmer_size:
        return False
    return True
    

# def get_kmers_from_taxon(fasta, kmer_size, n_kmers):
#     len_of_file = get_eof(fasta)
#     if len_of_file < kmer_size:
#         LOGGER.warn("No kmers found for taxid {0}".format(fasta))
#         return set()
#     with open(fasta, 'r+b') as handle:
#         good_kmers = set()
#         mm = mmap.mmap(handle.fileno(), 0, mmap.PROT_READ)
#         n_tries = 0
#         while len(good_kmers) < n_kmers and n_tries < 3:
#             # over sample because some will be lost in QC
#             for i in range(n_kmers + n_kmers//2):
#                 position = np.random.randint(len_of_file)
#                 mm.seek(position)
#                 kmer = mm.read(kmer_size)
#                 if not valid_kmer(kmer, kmer_size):
#                     continue
#                 good_kmers.add(kmer)
#             n_tries += 1
#         mm.close()
#     good_kmers = list(good_kmers) if len(
#         good_kmers) <= n_kmers else list(good_kmers)[:n_kmers]
#     LOGGER.info(
#         "Sampled {0} kmers for taxid {1}".format(
#             len(good_kmers), fasta))
#     return good_kmers


def get_sample_kmers(
    taxa, kmer_size, n_kmers,
    outpath, fasta_path, pickle_path, threads):
    taxa = np.array(np.genfromtxt(taxa, dtype=str), ndmin=1)
    total_samples = len(taxa)
    LOGGER.info(
        "Generating random kmers for candidate taxa:\n {}".format(
            "\n".join(taxa)))
    LOGGER.info("deserializing {}".format(pickle_path))
    tx_ids, _, positions = deserialization(pickle_path)
    LOGGER.info("Finished deserialization")
    # write_seq_partial = partial(
    #     write_sequences_to_file,
    #         fasta=fasta_path,
    #         tax_ids=tx_ids,
    #         positions=positions,
    #         outpath=os.path.dirname(outpath)
    #     )
    # p = Pool(threads)    
    # seq_files = p.map(write_seq_partial, taxa)
    # p.close()
    # p.join()
    
    get_kmer_partial = partial(
        get_kmers_from_taxon,
        kmer_size=kmer_size,
        n_kmers=n_kmers,
        positions=positions,
        tax_ids=tx_ids,
        fasta=fasta_path
        )

    
    p = Pool(threads)
    LOGGER.info(
        "Starting {0} threads to generate {1} "
        "random kmers of size {2} from {3} candidate taxa".format(
            threads, n_kmers, kmer_size, len(taxa)
        ))
    
    LOGGER.info("Writing to file: {}".format(outpath))

    with open(outpath, 'wb') as handle:
        for sample, kmer in enumerate(p.imap(
            get_kmer_partial, taxa, chunksize=1)):
            LOGGER.info("Writing sample {} to file".format(sample))
            handle.write(
                get_outstring(
                    kmer, sample, total_samples))
            LOGGER.info("Done writing sample {} to file".format(sample))


    p.close()
    p.join()



def id_generator():
    i = 0
    while True:
        i += 1
        yield bytes(str(i), 'utf8')

ID = id_generator()



if __name__ == "__main__":
    try:
        config_logging(snakemake.log[0], "INFO")
        LOGGER = logging.getLogger(__name__)

        get_sample_kmers(
            snakemake.input[0],
            snakemake.params[0],
            snakemake.params[1],
            snakemake.output[0],
            snakemake.params[2],
            snakemake.params[3],
            snakemake.threads)
    except NameError:

        PARSER = argparse.ArgumentParser(
        prog="MTSv Random Kmers",
        description="Get random kmers from taxids.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )     

        PARSER.add_argument(
            "fasta", metavar="FASTAPATH", type=file_type,
            help="Path to database fasta file."
        )

        PARSER.add_argument(
            "taxa", metavar="TAXAFILE", type=file_type,
            help="File with list of taxa"
        )

        PARSER.add_argument(
            "pickle", metavar="PICKLEPATH", type=file_type,
            help="Path to database pickle file"
        )

        PARSER.add_argument(
            "outpath", metavar="OUTPATH", type=outfile_type,
            help="Output file"
        )

        PARSER.add_argument(
            "--threads", default=1, type=positive_int,
            help="Number of threads."
        )

        PARSER.add_argument(
            "--kmer_size", '-k', type=positive_int, default=50,
            help="Size of random kmers"
        )

        PARSER.add_argument(
            "--n_kmers", '-n', type=positive_int, default=100000,
            help="number of random kmers per taxid"
        )

        PARSER.add_argument(
            "--log", type=outfile_type, default="random_kmers.log",
            help="Path to log file."
        )

        ARGS = PARSER.parse_args()


        config_logging(ARGS.log, "INFO")
        LOGGER = logging.getLogger(__name__)

        get_sample_kmers(
            ARGS.taxa, ARGS.kmer_size, ARGS.n_kmers,
            ARGS.outpath, ARGS.fasta, ARGS.pickle, ARGS.threads)




        





