from Bio import SeqIO
from ete3 import NCBITaxa
import logging
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from mtsv.utils import line_generator
from mtsv.parsing import parse_output_row


def get_decendants(taxids, des_flag):
    parents = []
    descendants = []
    for taxid in taxids:
        if des_flag:
            des = list(NCBI.get_descendant_taxa(
                taxid, collapse_subspecies=True))
            if taxid not in des:
                des.append(taxid)
        else:
            des = [taxid]
        descendants += des
        parents += [taxid for _ in range(len(des))]
    return np.array(parents, dtype=int),
            np.array(descendants, dtype=int)


def get_targets(clps_data_chunk, parents, descendants):
    target_dict = defaultdict(list)
    for line in clps_data_chunk:
        data = parse_output_row(line)
        hits = np.where(np.in1d(descendants, data.taxa))[0]
        for hit in hits:
            target_dict[parents[hit]].append(data.read_name)
    return target_dict

def reduce_targets(taxids, targets):
    # ensure each taxid is in final dict even if there
    # were no hits
    all_targets = {tax:[] for tax in taxids}
    for target in targets:
        for k, v in target:
            all_targets[k].extend(v)
    return all_targets

def mtsv_extract(
    taxids, clps_file, query_fasta,
    descendants, by_sample, outpath, threads, log):
    # for higher level taxids, include all descendants in output
    logger = logging.getLogger(__name__)
    parents, descendants = get_decendants(taxids, descentants)
    targets_partial = partial(
                        get_targets,
                        parents=parents,
                        descendants=descendants)
    get_lines = line_generator(clps_file, 5000)
    p = Pool(threads)
    target_dicts = p.imap(targets_partial, get_lines)
    p.close()
    p.join()
    targets = reduce_targets(taxids, target_dicts)
    query_fasta_dict = SeqIO.to_dict(
        SeqIO.parse(query_fasta, "fasta"))
    
    fun = write_sequences_by_sample if by_sample else write_sequences
    write_partial = partial(fun,
                            query_fasta=query_fasta_dict,
                            outpath=outpath)
    p = Pool(threads)
    p.imap(write_partial, targets.items())
    

def write_sequences_by_sample(
    targets_dict, query_fastas, outpath):
    if not len(targets_dict[1]):
        logger.info(
            "There were no sequences for taxid: {}".format(
                targets_dict[0]))
        return   
    sample_dict = defaultdict(list)
    for query in targets_dict[1]:
        samples = np.where(
            np.array(s.split("_")[1:], dtype=int))[0]
        for sample in samples:
            sample_dict[sample].append(query)
    for sample, queries in sample_dict.items():
        file_name = os.path.join(
            outpath,
            "{0}_{1}.fasta".format(targets_dict[0], sample))
        with open(file_name, 'w') as handle:
            records = [query_fastas[query] for query in queries]
            SeqIO.write(records, handle, "fasta")               

        

def write_sequences(
    targets_dict, query_fastas, outpath):
    file_name = os.path.join(
        outpath, targets_dict[0] + ".fasta")
    with open(file_name, 'w') as handle:
        if not len(targets_dict[1]):
            logger.info(
                "There were no sequences for taxid: {}".format(
                    targets_dict[0]))
        else:
            records = [query_fastas[target] for target in targets_dict[1]]
            SeqIO.write(records, handle, "fasta")


if __name__ == "__main__":

    NCBI = NCBITaxa(taxdump_file=snakemake.params[0])

    mtsv_extract(
        snakemake.params[1],
        snakemake.input[1],
        snakemake.input[0],
        snakemake.params[3],
        snakemake.params[4],
        snakemake.params[2],
        snakemake.threads,
        snakemake.log[0])


