import logging
import time
import mmap
import os
import numpy as np
from multiprocessing import Pool
from functools import partial
from mtsv.utils import config_logging
from mtsv.mtsv_prep.MTSv_prune import deserialization, get_tree

DNA = set(b"ACTG")


def get_eof(file_name):
    with open(file_name, 'r') as infile:
        infile.seek(0, 2)     # go to the file end.
        return infile.tell()

def write_sequences_to_file(taxon, fasta, tax_ids, positions, outpath):
    outfile = os.path.join(outpath, str(taxon))
    taxons = set(get_tree(tax_ids, taxon, positions))
    tx_positions = []
    for tx in taxons:
        if tx:
            try:
                tx_positions += positions[tx.encode()]
            except KeyError:
                continue
    tx_positions.sort()
    with open(outfile, 'wb') as out: 
        with open(fasta, 'rb') as fasta_handle:
            for position in tx_positions:
                fasta_handle.seek(position)
                fasta_handle.readline()  # header
                line = fasta_handle.readline()
                while line and chr(line[0]) != ">":
                    out.write(line.replace(b"\n", b""))
                    line = fasta_handle.readline()
                out.write(b"\n")
    return outfile


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
    

def get_kmers_from_taxon(fasta, kmer_size, n_kmers):
    len_of_file = get_eof(fasta)
    print(len_of_file)
    if len_of_file < kmer_size:
        logger.warn("No kmers found for taxid {0}".format(fasta))
        return set()
    with open(fasta, 'r+b') as handle:
        good_kmers = set()
        mm = mmap.mmap(handle.fileno(), 0, mmap.PROT_READ)
        n_tries = 0
        while len(good_kmers) < n_kmers and n_tries < 3:
            # over sample because some will be lost in QC
            for i in range(n_kmers + n_kmers//2):
                position = np.random.randint(len_of_file)
                mm.seek(position)
                kmer = mm.read(kmer_size)
                if not valid_kmer(kmer, kmer_size):
                    continue
                good_kmers.add(kmer)
            print(len(good_kmers), n_tries)
            n_tries += 1
        mm.close()
    good_kmers = list(good_kmers) if len(
        good_kmers) <= n_kmers else list(good_kmers)[:n_kmers]
    logger.info(
        "Sampled {0} kmers for taxid {1}".format(
            len(good_kmers), fasta))
    return good_kmers


def get_sample_kmers(
    taxa, kmer_size, n_kmers,
    outpath, fasta_path, pickle_path, threads):
    taxa = np.array(np.genfromtxt(taxa, dtype=str), ndmin=1)
    logger.info(
        "Generating random kmers for candidate taxa:\n {}".format(
            "\n".join(taxa)))
    logger.info("deserializing {}".format(pickle_path))
    tx_ids, _, positions = deserialization(pickle_path)
    logger.info("Finished deserialization")
    write_seq_partial = partial(
        write_sequences_to_file,
            fasta=fasta_path,
            tax_ids=tx_ids,
            positions=positions,
            outpath=os.path.dirname(outpath)
        )
    p = Pool(threads)    
    seq_files = p.map(write_seq_partial, taxa)
    p.close()
    p.join()
    get_kmer_partial = partial(
        get_kmers_from_taxon,
        kmer_size=kmer_size,
        n_kmers=n_kmers
    )
    p = Pool(threads)
    logger.info(
        "Starting {0} threads to generate {1} "
        "random kmers of size {2} from {3} candidate taxa".format(
            threads, n_kmers, kmer_size, len(taxa)
        ))
    with open(outpath, 'w') as handle:
        for sample, kmers in enumerate(p.imap(get_kmer_partial, seq_files)):
            write_samples_to_file(handle, kmers, sample, len(taxa))
    logger.info("Writing to file: {}".format(outpath))
    p.close()
    p.join()



def id_generator():
    i = 0
    while True:
        i += 1
        yield i

ID = id_generator()

def write_samples_to_file(handle, kmers, sample, total_samples):
    count_arr = np.array(np.zeros(total_samples, dtype=int), dtype=str)
    count_arr[sample] = "1"
    count_arr = "_".join(count_arr)
    for kmer in kmers:
        handle.write(">R{0}_{1}\n{2}\n".format(
            next(ID), count_arr, kmer.decode('ascii')))


if __name__ == "__main__":
    config_logging(snakemake.log[0], "INFO")
    logger = logging.getLogger(__name__)

    get_sample_kmers(
        snakemake.input[0],
        snakemake.params[0],
        snakemake.params[1],
        snakemake.output[0],
        snakemake.params[2],
        snakemake.params[3],
        snakemake.threads)
