import pandas as pd
import numpy as np
from ete3 import NCBITaxa

def get_sample_count(n_cols):
    n_cols -= 3
    return n_cols//4

def get_candidate_taxa(summary_file, outfile, signature_cutoff):
    df = pd.read_csv('summary_file')
    n_samples = get_sample_count(df.shape[1])
    taxa_set = set()
    for sample in range(n_samples):
        d = df[df[
            "Unique Signature Hits (S{})".format(
                sample + 1)] > signature_cutoff]
        taxa_set.update(list(d.TaxID))
    ranks = NCBI.get_rank(list(taxa_set))
    taxa = [k for k, v in ranks.items() if v == "species"]
    np.savetxt(outfile, taxa, fmt="%d")

if __name__ == "__main__":
    NCBI = NCBITaxa(snakemake.params[1])
    get_candidate_taxa(
        snakemake.input[0],
        snakemake.output[0],
        snakemake.params[0])
