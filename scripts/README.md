# Setup Environment on Monsoon
### Load the Python3 version of Anaconda
```
$ module load anaconda/3.latest
```

### Create Conda environment with Python3 and required packages.
The environment only needs to be created once.
```
$ conda create --name biopy3 python=3.5.2 biopython pandas
```
The name `biopy3` can be changed to whatever name you want.  

### Activate Conda Environment
```
$ source activate biopy3
```

### Install etetoolkit
```
$ conda install -c etetoolkit ete3 ete_toolchain
```

### Deactivate Conda Environment
```
$ source deactivate
```

If you already have a Python3.5.2 environment you can just use conda to install the required packages.
```
$ source activate ENV_NAME
$ conda install -c etetoolkit ete3 ete_toolchain
$ conda install pandas
$ conda install biopython
```
# Pre-Processing
## MTSv Prune
The 'MTSv_prune.py' is a ***work in progress*** module. The first step is to acquire the necessary GenBank flat files from NCBI
along with two data stores are needed to partition sequence data by NCBI taxonomy

### GenBank Flat file Download
This command will download and store needed files in a "./raw/" subfolder
```
python MTSv_prune.py -p -t <optional:threads>
```

### Fasta Database
```
python MTSv_prune.py -bdb -fl <file-list> \ 
-o <output name without extension> -t <threads:default 1>  
```
### Index Building
```
python MTSv_prune.py -biacc -a2t <list of *accession2taxi.seq.gz> \
-tp <path to taxdump.tar.gz> -fp <fasta database path> \ 
-o <name of serialization without extension>```
```
### Clipping
To obtain sequences associated with an NCBI taxonomic subtree use command
```
python MTSv_prune.py -c -txi <list of taxids to include> \
-txe <list of taxids to exclude> -fp <fasta database path> \
-sp <path to index> -rur <taxonomic rank to assign sequence in output> \
-o <output name with extension>
```
### Configuration JSON
Many of the parameters in pruning will not change often so  to save time a configuration file can be created
```
python MTSv_prune.py -gc -fp <fasta database path> -sp <path to index> \
-rur <taxonomic rank to assign sequence in output> -o <configuration path w/out ext>
```
This saved configuration file can then be accessed during clipping with the command
```
python MTSv_prune.py -c -cp <Config JSON path> -txi <list of taxids to include> \
-txe <list of taxids to exclude> -o <output name with extension>
```

### Monsoon PreBuilt Fasta Database and Configurations
```
/scratch/tes87/database/assembly_levels/nt_2016_seqs.json
/scratch/tes87/database/assembly_levels/Complete_Genome_assembly.json
/scratch/tes87/database/assembly_levels/Chromosome_assembly.json
/scratch/tes87/database/assembly_levels/Scaffolds_assembly.json
```


# MTSv Summary
The `MTSv_summary.py` script summarizes the number of hits per taxon per sample (sample colunns are in the same order as the FASTQ files that were passed to `MTSv-readprep`). The total number of reads mapped, the number of unique mapped reads, the number of signature hits, and the number of unique signature hits per taxon. Note: if the `--lca` option is modified to combine species up to the genus or family level in the `MTSv-inform` module, taxid for uncombined species in the `MTSv-collapse` output will be double counted in total hits and unique hits. 

### Input
**Required Positional Arguments**  
`project_name`: Provide a prefix that will be used to name all output files.  
`collapse_file`: Path to MTSv-collapse output file  
`signature_file`: Path to MTSv-inform output file  
**Optional Arguments**  
`--out_path`: Directory to write output (Default: ./)  
`--update`: Updates the ncbi taxonomy database by downloading and parsing the latest taxdump.tar.gz file from the NCBI FTP site (via HTTP).  
`--taxdump`: Path to alternative taxdump.tar.gz file.

### Output

A csv file <`out_path`/`program_name`_summary.csv> in the following format:  


| TaxID | Division | Sci. Name          | Total Hits (S1) | Unique Hits (S1) | Signature Hits (S1) | Unique Signature Hits (S1) | Total Hits (S2) | Unique Hits (S2) | Signature Hits (S2) | Unique Signature Hits (S2) |
|-------|----------|--------------------|-----------------|------------------|---------------------|----------------------------|-----------------|------------------|---------------------|----------------------------|
| 1392  | Bacteria | Bacillus anthracis | 12              | 2                | 2                   | 1                          | 13              | 1                | 0                   | 0                          |
| 1396  | Bacteria | Bacillus cereus    | 10              | 1                | 0                   | 0                          | 13              | 1                | 0                   | 0                          |

These results correspond to the following input files:

`sample_collapse.txt`
```
R1_10_13:1392,1396
R2_2_0:1392
```
`sample_signature.txt`
```
R2_2_0:1392
```

### Usage
```
$ python MTSv_summary.py --help
usage: MTSv Summary [-h] [-o OUT_PATH] [--update] [--taxdump TAXDUMP]
                    PROJECT_NAME COLLAPSE_FILE SIGNATURE_FILE

Summarize number of hits for each taxa, including signature hits.

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  COLLAPSE_FILE         Path to MTSv-collapse output file
  SIGNATURE_FILE        Path to MTSv-inform output file

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Output directory (default: ./)
  --update              Update taxdump (default: False)
  --taxdump TAXDUMP     Alternative path to taxdump. Default is home directory
                        where ete3 automatically downloads the file. (default:
                        None)
```
### Example Slurm Script
Change nauid to your nauid and modify `out_path` to test and run script  

```
#!/bin/bash
#SBATCH --job-name=MTSv-summary
#SBATCH --output=/scratch/nauid/output.txt 

module load anaconda/3.latest
source activate biopy3

python -u MTSv_summary.py test \
/scratch/tf362/vedro/merge/merged_results.txt \
/scratch/tf362/vedro/inform/informative.txt \
--out_path /scratch/nauid/path/to/output/ \
```

# MTSv Extract
The `MTSv_extract.py` script extracts all unique read sequences that aligned to a provided taxid or species name.

### Input
**Required Positional Arguments**  
`project_name`: provide prefix that will be used to name all output files  
`collapse_file`: Path to MTSv-collapse output file  
`signature_file`: Path to MTSv-inform output file  
`reads_fasta`: Path to FASTA file produced from MTSv-readprep  

**Required Mutually Exclusive Arguments**  
`taxid`: the taxid to extract  
`species`: the species name to extract

**Optional Arguments**  
`--out_path`: Directory to write output (Default: ./)  
`--update`: Updates the ncbi taxonomy database by downloading and parsing the latest taxdump.tar.gz file from the NCBI FTP site (via HTTP).  
`--taxdump`: Path to alternative taxdump.tar.gz file.


### Output
`PROJECT_NAME_TAXID_all.fasta`: Contains all sequence reads that aligned to taxid.  
`PROJECT_NAME_TAXID_signature.fasta`: Contains sequence reads that only aligned to this taxid and no other taxids.

### Usage
```
$ python MTSv_extract.py --help
usage: MTSv Extract [-h] [-o OUT_PATH] [--update] [--taxdump TAXDUMP]
                    (-t TAXID | -s SPECIES)
                    PROJECT_NAME COLLAPSE_FILE SIGNATURE_FILE READS_FASTA

Extracts all sequences, including signature hits, that aligned to a given
taxid or species name.

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  COLLAPSE_FILE         Path to MTSv-collapse output file
  SIGNATURE_FILE        Path to MTSv-inform output file
  READS_FASTA           Path to FASTA file from MTSv-readprep

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Output directory (default: ./)
  --update              Update taxdump (default: False)
  --taxdump TAXDUMP     Alternative path to taxdump. Default is home directory
                        where ete3 automatically downloads the file. (default:
                        None)
  -t TAXID, --taxid TAXID
                        Extract sequences by taxid (default: None)
  -s SPECIES, --species SPECIES
                        Extract sequences by species name (default: None)
```
### Example Slurm Script
Change nauid to your nauid and modify `out_path` to test and run script
```
#!/bin/bash
#SBATCH --job-name=MTSv-extract
#SBATCH --output=/scratch/nauid/output.txt      
module load anaconda/3.latest
source activate biopy3

python -u MTSv_extract.py test \
/scratch/tf362/vedro/merge/merged_results.txt \
/scratch/tf362/vedro/inform/informative.txt \
/scratch/tf362/vedro/readprep/seg50_minqual15-qualthresh-3.fasta \
--out_path /scratch/nauid/path/to/output/ \
--taxid 9606 # human
```
