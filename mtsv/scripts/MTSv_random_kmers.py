import argparse
import logging
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from mtsv.utils import config_logging
from mtsv.mtsv_prep.MTSv_prune import deserialization, get_tree
from mtsv.parsing import file_type, outfile_type, positive_int

DNA = set(b"ACTG")

POSITIONS = None
SORTED_POSITIONS = None
TX_IDS = None
EOF = None

def get_eof(file_name):
    with open(file_name, 'r') as infile:
        infile.seek(0, 2)     # go to the file end.
        return infile.tell()


def get_counts(lens, total, n_samples):
    ns = np.random.multinomial(
        n_samples, pvals=[l/total for l in lens])
    remainder = 0
    for i, (n, l) in enumerate(zip(ns, lens)):
        if n > l:
            ns[i] = l
            remainder += (n - l)
        take_from_remainder = min(remainder, l - n)
        remainder -= take_from_remainder
        ns[i] += take_from_remainder
    return ns

def check_range(ranges, kmer_size):
    # adjust the end of the range
    # so that we don't go off the end
    good_ranges = []
    for start, stop in zip(ranges[0], ranges[1]):
        stop = stop - kmer_size
        if stop > start:
            good_ranges.append((start, stop))
    return np.array(good_ranges, dtype=np.uint64)


def get_random_positions(ranges, n_samples, kmer_size):
    ranges = check_range(ranges, kmer_size)
    n_samples = n_samples * 2 # over sample as some will be lost due to duplicates
    lengths = [r[1] - r[0] for r in ranges]
    total = sum(lengths)
    if n_samples >= total:
        positions = []
        for r in ranges:
            positions += range(r[0], r[1])
        return np.array(positions, dtype=int)
    positions = []
    # this may not add up perfectly to n_samples
    # but we over sampled so this shouldn't
    # be a problem.
    ns = get_counts(lengths, total, n_samples)
    for n, rng in zip(ns, ranges):
        positions += list(np.random.choice(
            np.arange(
                rng[0], rng[1], dtype=int),
                size=n, replace=False))
    positions = np.array(positions, dtype=int)
    positions.sort()
    return positions

def get_outstring(kmers, sample, total_samples):
    count_arr = np.array(
        np.zeros(total_samples, dtype=int), dtype=bytes)
    count_arr[sample] = b"1"
    count_arr = b"_".join(count_arr)
    s = b""
    for kmer in kmers:
        s += b">R" + next(ID) + b"_" + count_arr + b"\n" + kmer + b"\n"
    return s
    
def get_start_positions_for_taxon(
    taxon, fasta, tax_ids, positions, n_kmers):
    taxons = set(get_tree(tax_ids, taxon, positions))
    tx_positions = []
    for tx in taxons:
        if tx:
            try:
                tx_positions += positions[tx.encode()]
            except KeyError:
                continue
    if len(tx_positions) > 5 * n_kmers:
        tx_positions = np.random.choice(
            tx_positions, size=5 * n_kmers, replace=False)
    tx_positions.sort()
    with open(fasta, 'rb') as fasta_handle:
        start_positions = []
        seek = fasta_handle.seek
        readline = fasta_handle.readline
        tell = fasta_handle.tell
        append = start_positions.append
        for position in tx_positions:
            seek(position)
            readline()
            append(tell())
    return np.array(start_positions, dtype=np.uint64)        

def get_end_position_for_taxon(
    start_positions, sorted_positions, eof):
    end_positions = []
    append = end_positions.append
    length_positions = len(sorted_positions)
    for pos in start_positions:
        end = np.searchsorted(
            sorted_positions, np.uint64(pos), side='right')
        
        if end < length_positions:
            append(sorted_positions[end])
        else:
            append(eof)
    return np.array(end_positions, dtype=np.uint64)
        

def get_kmers(positions, fasta, kmer_size, n_kmers):
    with open(fasta, 'rb') as fasta_handle:
        seek = fasta_handle.seek
        read = fasta_handle.read
        replace = bytes.replace
        kmers = []
        for position in positions:
            seek(position)
            kmer = replace(
                read(kmer_size * 2),
                b"\n", b"")[:kmer_size]
            if valid_kmer(kmer, kmer_size):
                kmers.append(kmer)
        kmers = np.unique(kmers)
        np.random.shuffle(kmers)
        kmers = kmers[:n_kmers]
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


def get_sorted_positions():
    pos = []
    for p in POSITIONS.values():
        pos += p
    pos = np.array(pos, dtype=np.uint64)
    pos.sort()
    return pos


def get_kmers_from_taxon(
    taxon, fasta_path, n_kmers, kmer_size):
    LOGGER.info("Sampling kmers from taxid {}".format(taxon))
    start_positions = get_start_positions_for_taxon(
        taxon, fasta_path, TX_IDS, POSITIONS, n_kmers)

    return get_kmers(get_random_positions(
        (start_positions,
        get_end_position_for_taxon(
        start_positions, SORTED_POSITIONS, EOF)),
        n_kmers, kmer_size
    ), fasta_path, kmer_size, n_kmers)
    


def get_sample_kmers(
    taxa, kmer_size, n_kmers,
    outpath, fasta_path, pickle_path, threads):
    taxa = np.array(np.genfromtxt(taxa, dtype=str), ndmin=1)
    total_samples = len(taxa)
    LOGGER.info(
        "Generating random kmers for candidate taxa:\n{}".format(
            "\n".join(taxa)))
    LOGGER.info("Deserializing {}".format(pickle_path))
    global TX_IDS
    global POSITIONS
    global SORTED_POSITIONS
    global EOF
    TX_IDS, _, POSITIONS = deserialization(pickle_path)
    EOF = get_eof(fasta_path)
    SORTED_POSITIONS = get_sorted_positions()
    LOGGER.info("Finished deserialization")

    LOGGER.info(
        "Starting {0} threads to generate {1} "
        "random kmers of size {2} from {3} candidate taxa".format(
            threads, n_kmers, kmer_size, len(taxa)
        ))

    get_kmers_partial = partial(
        get_kmers_from_taxon, 
        fasta_path=fasta_path,
        n_kmers=n_kmers,
        kmer_size=kmer_size)

    p = Pool(threads)
    LOGGER.info("Writing to file {}".format(outpath))
    with open(outpath, 'wb') as handle:
        for sample, kmer in enumerate(p.imap(
            get_kmers_partial, taxa, chunksize=1)):
            taxid = taxa[sample]
            LOGGER.info("Found {0} kmers for taxid {1}.".format(
                len(kmer), taxid))
            LOGGER.info("Writing samples for taxid {} to file".format(
                taxid))
            handle.write(
                get_outstring(
                    kmer, sample, total_samples))
            LOGGER.info(
                "Done writing sample for taxid {} to file".format(
                    taxid))

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




        




