# MTSv

MTSv is a suite of metagenomic binning and analysis tools. It attempts to accurately identify which species are present in a given DNA sample. It assumes that read fragments in samples will be in a "shotgun" or short read format, typically ~50-200 bases in length.

All commands listed in this document assume they're being executed from the repository's root directory. Adjust accordingly if you've installed the tools elsewhere.

## Building

MTSv is built in Rust, with a little bit of Python. You'll need:

* `rustc` and `cargo` >= 1.8.0 ([rustup.rs](https://rustup.rs) is the easiest installation method)
* a C compiler (tested with GCC and clang)

### Tests

To run tests:

~~~
$ cargo test
~~~

To generate a code coverage report, make sure [kcov >= 26](https://simonkagstrom.github.io/kcov/) is installed on your `PATH`, then install `cargo-kcov`:

~~~
$ cargo install cargo-kcov
~~~

To run coverage:

~~~
$ cargo kcov -- --exclude-pattern="/.cargo,vendor/,tests/,bench/,include/,bin/,ssw/"
~~~

This will place a code coverage report under `target/kcov/index.html`.

### Building

To build the MTSv binaries:

~~~
$ cargo build --release
~~~

They'll be available under `target/release/MTSv-*`.

## Documentation

To generate the internal documentation:

~~~
$ cargo doc [--open]
~~~

(pass the `--open` flag if you want to immediately open the docs in your browser)

## Usage

MTSv builds several binaries:

* `MTSv-chunk`
* `MTSv-binner`
* `MTSv-build`
* `MTSv-collapse`
* `MTSv-inform`
* `MTSv-readprep`
* `MTSv-tree-build`

All of these accept the `--help` flag to print a help message on their usage. See below for specific use instructions.

### Snakemake

There's an (in-progress) Snakefile in the repository root. This manages a MTSv-based metagenomic workflow. See the [Snakemake website](https://bitbucket.org/snakemake/snakemake/wiki/Home) for installation and usage instructions.

(TODO: include instructions on configuring snakefile for particular file sets)

### Index construction

MTSv uses several pre-constructed indices for running its queries.

#### Reference file format & taxdump.tar.gz

To construct the indices, you'll need two files:

1. A FASTA file of all reference sequences, with headers in the format `ACCESSION-TAXONOMICID`. So if a sequence has accession # 12345, and belongs to the NCBI taxonomic ID 987, the header for that sequence should read `12345-987`.
2. The `taxdump.tar.gz` file from NCBI which corresponds to the sequences in your FASTA file.

#### Chunking reference database

MTSv uses A LOT of memory for its indices. About 20x the space compared to the FASTA file its given. As a result, it's generally preferable to split the database into small chunks that can be processed iteratively. These chunks should, as much as possible, have all or most of a taxonomic ID in each of them, as MTSv achieves speedups by skipping queries once it's found a successful match in a taxonomic node. MTSv includes a utility for doing so. To split your reference database into 1GB chunks (resulting in 15-20GB needed for running queries):

~~~
$ target/release/MTSv-chunk -i PATH_TO_FASTA -o PATH_TO_CHUNK_FOLDER -g NUM_GBS_PER_CHUNK
~~~

This will write a series of chunk files into the directory specified. See the help message for further information.

#### Metagenomic index

Now that you have N chunks of your FASTA database, they need to be processed into indices which MTSv can use for querying.

~~~
$ target/release/MTSv-build --fasta /path/to/chunkN.fasta --index /path/to/write/chunkN.index
~~~

Using default settings, indices will use ~15-20x as much RAM as the reference file used for their creation (at a sampling interval of 512 bytes). This can be overridden by passing the `--sample-interval <FM_SAMPLE_INTERVAL>` flag. Lower than 512 will increase the size of the index and can provide a reduction in query time. Increasing the flag will decrease the size of the index up to a point (the size of the suffix array can't be reduced, this setting only changes the FM index size) while accepting a slower query time.

See the help message for other options.

#### Taxonomic tree index

To determine which reads are informative for particular taxonomic IDs, you'll need to construct an index from the NCBI `taxdump.tar.gz` which corresponds to your FASTA database.

~~~
$ target/release/MTSv-tree-build --dump /path/to/taxdump.tar.gz --index /path/to/write/tree.index
~~~

See the help message for other options.

### readprep

MTSv assumes that unidentified read fragments/sequences (referred to as "query reads") come in FASTA format and are of uniform length within a given query file. Often one needs to run some quality-control processes and combine several files. If you have a variety of FASTQ files to combine and QC, run `MTSv-readprep`. For the full list of configurations, see the help message:

~~~
$ target/release/MTSv-readprep --help
~~~

This will write reads with headers in the format `R_COUNT1_COUNT2_COUNT3` and so on, where each count corresponds to the number of times that read was found in each FASTQ file, in the order they were passed as arguments.

### Binning queries

Now that you have the indices MTSv needs, and have prepared the query reads, run `MTSv-binner` on each index chunk. In this example, 3 SNPs are tolerated in matches, and 8 threads are used for processing:

~~~
$ target/release/MTSv-binner --edits 3 --threads 8 \
    --index /path/to/chunkN.index \
    --fasta /path/to/prepared_reads.fasta \
    --results /path/to/write/chunkN_results.txt
~~~

See the help message for other options.

#### Output

`MTSv-binner` writes results for a single read per line. For example, if a read with the header `R1_0_1` maps to taxon IDs `562`, `9062`, and `100`:

~~~
R1_0_1:562,9062,100
~~~

### Collapsing chunks

Since each results file from the binner will only represent some of the species matches for a given query read, combine all of the chunked results into a single results file for further analysis:

~~~
$ target/release/MTSv-collapse /path/to/chunk1_results.txt /path/to/chunk2_results.txt ... \
    --output /path/to/collapsed_results.txt
~~~

Make sure to include all of the chunk files. While the collapser could be run in multiple phases, it's generally much faster to do them all at once.

See the help message for other options.

### Finding informative reads

At this point, you should have a single file which records all of the taxonomic IDs (species, usually) which your query reads have mapped to (within the specified number of edits, that is).

To determine which reads are "informative" for which taxonomic nodes, you'll run the `MTSv-inform` tool. The simplest way to do so (spawning 8 threads):

~~~
$ target/release/MTSv-inform \
    --index /path/to/tree.index \
    --input /path/to/collapsed_results.txt \
    --threads 8
    --lca 0
    --output /path/to/write/informatives.txt
~~~

The sensitivity of the analysis can be adjusted either by changing the `--lca` flag, by by specifying one of the `[--genus|--family]` flags. Changing the LCA (least common ancestor) value will affect how many jumps "up" the taxonomic tree to consider. Specifying a "logical" LCA flag will search for a common genus or family (depending on the flag) which covers all of the results found earlier.

# Analyzing output

## MTSv Summary
The `MTSv_summary.py` script summarizes the number of hits per taxon per sample (sample colunns are in the same order as the FASTQ files that were passed to `MTSv-readprep`). The total number of reads mapped, the number of unique mapped reads, the number of signature hits, and the number of unique signature hits per taxon. Note: if the `--lca` option is modified to combine species up to the genus or family level in the `MTSv-inform` module, taxid for uncombined species in the `MTSv-collapse` output will be double counted in total hits and unique hits. 

### Requirements
biopython=1.68
python=3.5.2  
pandas=0.22.0
ete3=3.1.1
ete_toolchain=3.0.0


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

