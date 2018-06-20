import random
import numpy as np
from multiprocessing import Pool
from functools import partial
from mtsv_prep.MTSv_prune import deserialization, get_tree

strip = str.strip
DNA = set("ACTG")

def get_kmers_from_taxon(taxon, kmer_size, n_kmers, fasta_path, pickle_path):
    tx_ids, _, positions = deserialization(pickle_path)
    taxons = set(get_tree(tx_ids, taxon, positions))
    kmer_reservoir = set()
    n = 0
    with open(fasta_path, 'rb') as fasta:
        for tx in sorted(taxons):
            if tx:
                tx = tx.encode().strip()
                try:
                    pos = positions[tx]
                except KeyError:
                    continue
                pos.sort()
                for off in pos:
                    fasta.seek(off)
                    fasta.readline() # header
                    buff = strip(fasta.readline())
                    while True:
                        line = strip(fasta.readline())
                        if not line:
                            break
                        if chr(line[0]) == ">":
                            break
                        buff += line
                    # fill up reservoir
                    while len(kmer_reservoir) < n_kmers and len(buff) >= kmer_size:
                        k = buff[:kmer_size]
                        buff = buff[1:]
                        if set(k).issubset(DNA) and k not in kmer_reservoir:
                            kmer_reservoir.add(k)
                            n += 1
                    # sample reservoir
                    while len(buff) >= kmer_size:
                        k = buff[:kmer_size]
                        buff = buff[1:]
                        # don't consider invalid strings to have
                        # been seen (don't increment n) to ensure
                        # uniform sampling over collection of
                        # valid kmers
                        if set(k).issubset(DNA) and k not in kmer_reservoir:
                            n += 1
                            if np.random.binomial(1, n_kmers/(n+1)):
                                remove = random.sample(kmer_reservoir, 1)
                                kmer_reservoir.remove(remove)
                                kmer_reservoir.add(k)
                        
    return kmer_reservoir


def get_sample_kmers(
    taxa, kmer_size, n_kmers,
    outpath, fasta_path, pickle_path, threads):
    # fasta path
    # pickle path
    taxa = np.genfromtxt(taxa, dtype=int)
    get_kmer_partial = partial(
        get_kmers_from_taxon,
        kmer_size=kmer_size,
        n_kmers=n_kmers,
        fasta_path=fasta_path,
        pickle_path=pickle_path
    )
    p = Pool(threads)
    with open(outpath, 'w') as handle:
        for sample, kmer_set in enumerate(p.imap(get_kmer_partial, taxa)):
            write_samples_to_file(handle, kmer_set, sample, len(taxa))
    p.close()
    p.join()


def id_generator():
    i = 0
    while True:
        i += 1
        yield i

ID = id_generator()

def write_samples_to_file(handle, kmer_set, sample, total_samples):
    count_arr = np.array(np.zeros(total_samples, dtype=int), dtype=str)
    count_arr[sample] = "1"
    count_arr = "_".join(count_arr)
    for kmer in kmer_set:
        handle.write(">R{0}_{1}\n{2}\n".format(
            next(ID), count_arr, kmer.decode('ascii')))


if __name__ == "__main__":

    get_sample_kmers(
        snakemake.input[0],
        snakemake.params[0],
        snakemake.params[1],
        snakemake.output[0],
        snakemake.params[2],
        snakemake.params[3],
        snakemake.threads)
