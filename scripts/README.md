# Setup Environment on Monsoon
### Load the Python3 version of Anaconda
```
$ module load anaconda/3.latest
```

### Create Conda environment with Python3 and Biopython
Only needs to be created once and only if you don't already have an environment with python3 and biopython. 
```
$ conda create --name biopy3 python=3.5.2 biopython
```
The name `biopy3` can be changed to whatever name you want.

### Activate Conda Environment
```
$ source activate biopy3
```

### Deactivate Conda Environment
```
$ source deactivate
```

# MTSv Extract
The `MTSv_extract.py` script extracts all read sequences that aligned to a provided taxid.

### Input
**Required Positional Arguments**  
`project_name`: provide prefix that will be used to name all output files  
`taxid`: the taxid to extract  
`collapse_file`: Path to MTSv-collapse output file  
`signature_file`: Path to MTSv-inform output file  
`reads_fasta`: Path to FASTA file produced from MTSv-readprep  

**Optional Arguments**  
`out_path`: Directory to write output (Default: ./)

### Output
`PROJECT_NAME_all.fasta`: Contains all sequence reads that aligned to taxid.  
`PROJECT_NAME_signature.fasta`: Contains sequence reads that only aligned to this taxid and no other taxids.

### Usage
```
python MTSv_extract.py --help
usage: MTSv Extract [-h] [-o OUT_PATH]
                    PROJECT_NAME TAXID COLLAPSE_FILE SIGNATURE_FILE
                    READS_FASTA

Extracts all read sequences, including signature reads, that aligned to a
given taxid.

positional arguments:
  PROJECT_NAME          Project name and output file prefix
  TAXID                 Taxid to extract
  COLLAPSE_FILE         Path to MTSv-collapse output file
  SIGNATURE_FILE        Path to MTSv-inform output file
  READS_FASTA           Path to FASTA file from MTSv-readprep

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_PATH, --out_path OUT_PATH
                        Output directory (default: ./)
```
### Example Slurm Script
```
#!/bin/bash
#SBATCH --job-name=MTSv-extract
#SBATCH --output=/scratch/nauid/output.txt      
module load anaconda/3.latest
source activate biopy3

python -u MTSv_extract.py test_1392 1392 \ /scratch/tf362/vedro/merge/merged_results.txt \
/scratch/tf362/vedro/inform/informative.txt \
/scratch/tf362/vedro/readprep/seg50_minqual15-qualthresh-3.fasta \
--output /scratch/nauid/path/to/output/
```