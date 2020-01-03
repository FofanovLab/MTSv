import logging
import os
import mmap
import click
import numpy as np
from mtsv.mtsv_prep.MTSv_prune import (deserialization, get_tree)

def random_kmers(taxid, fasta, pickle, kmer, n_kmers, outfile):
    logging.info("Deserializing {}".format(pickle))
    tx_ids, positions = deserialize_pickle(pickle)
    logging.info("Finished Deserializing")
    mmap = get_mmap(fasta)
    logging.info(
        "Sampling {0} {1}-mers from taxid {2}".format(
            n_kmers, kmer, taxid))

    kmers = set()
    while len(kmers) < n_kmers:
        kmers.update(get_random_kmers_from_position(
            mmap, taxid, tx_ids, positions, n_kmers, kmer))

    kmers = list(kmers)[:n_kmers]
    logging.info("Finished sampling.")
    logging.info("Writing to {0}".format(outfile))
    write_kmers_to_fasta(kmers, outfile)


def get_random_kmers_from_position(
        mmap, taxid, tx_ids, positions, n_kmers, kmer):
    kmers = []
    random_positions = get_random_start_positions_for_taxon(
        taxid, tx_ids, positions, n_kmers)
    for start, n in random_positions:
        kmers += random_kmers_from_sequence(
            start, n, kmer, mmap)
    return kmers


def write_kmers_to_fasta(kmers, outfile):
    fasta_record = "\n".join(
        [">R{0}_1\n{1}".format(i, kmer.decode())
         for i, kmer in enumerate(kmers)])
    with open(outfile, 'w') as out:
        out.write(fasta_record + "\n")


def deserialize_pickle(pickle):
    tx_ids, _, positions = deserialization(pickle)
    return tx_ids, positions


def get_random_start_positions_for_taxon(
        taxon, tax_ids, positions, n_kmers):
    """
    Returns a list of tuples representing the start position
    of random sequences for taxon and the number of kmers to
    sample from the sequence. The number of samples is 2X
    n_kmers to allow for loss due to duplication.  
    """
    taxons = set(get_tree(tax_ids, taxon, positions))
    tx_positions = []
    for tx in taxons:
        if tx:
            try:
                tx_positions += positions[tx.encode()]
            except KeyError:
                continue
    try:
        return list(zip(*np.unique(np.random.choice(
            tx_positions, size=2 * n_kmers, replace=True),
            return_counts=True)))
    except ValueError:
        if len(tx_positions) == 0:
            logging.warning(
                """
                Taxid {0} not in database, this should not
                if the same database is used for all steps""".format(taxon))
        return tx_positions


def get_mmap(fasta):
    mfd = os.open(fasta, os.O_RDONLY)
    return mmap.mmap(mfd, 0, prot=mmap.PROT_READ)


def random_kmers_from_sequence(position, n_sample, kmer_size, mmap):
    start = get_sequence_start(mmap, position)
    end = get_sequence_end(mmap, start)
    # subtract kmer_size from length to get only valid start positions
    seq_len = get_length_of_sequence(mmap, start, end) - kmer_size
    try:
        sample_positions = np.random.randint(
            0, seq_len, size=n_sample)
    except ValueError:  # if sequences is too short to have a range
        return []
    return get_kmers_at_positions(
        mmap, sample_positions, start, end, kmer_size)


def get_kmers_at_positions(
        mmap, sample_positions, start, end, kmer_size):
    """
    Return kmers from sequence and sampled positions
    """
    seq = mmap[start + 1: end].replace(b"\n", b"")
    return [
        seq[pos: pos + kmer_size] for
        pos in sample_positions]


def get_sequence_start(mmap, position):
    """
    Get start of sequence at position (after header)
    """
    return mmap.find(b"\n", position + 1)


def get_sequence_end(mmap, start):
    """
    Get end of sequence
    """
    return mmap.find(b">", start + 1)


def get_length_of_sequence(mmap, start, end):
    """
    Return length of sequence after removing
    newline characters
    """
    return len(
        mmap[start + 1: end].replace(b"\n", b""))


@click.command()
@click.option(
    '--fasta', '-f', type=click.Path(exists=True, dir_okay=False),
    required=True, help="Fasta database file.")
@click.option(
    '--serial', '-s', type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Serialized file from fasta.")
@click.option(
    '--outfile', '-o', type=click.Path(exists=True, writable=True, dir_okay=False),
    required=True, help="Output fasta file.")
@click.option(
    '--taxid', '-t', type=click.INT, required=True, help="Taxid to sample.")
@click.option('--kmer_size', default=50, help="Size of kmers to sample.")
@click.option('--n_kmers', default=100000, help="Number of kmers to sample.")
def main(fasta, serial, outfile, taxid, kmer_size, n_kmers):
    random_kmers(
        taxid, fasta, serial, kmer_size, n_kmers, outfile)

if __name__ == '__main__':

    main()
