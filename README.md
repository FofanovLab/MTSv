[![Build Status](https://travis-ci.org/FofanovLab/MTSv.svg?branch=master)](https://travis-ci.org/FofanovLab/MTSv)
[![Anaconda-Server Badge](https://anaconda.org/fofanov/mtsv/badges/installer/conda.svg)](https://conda.anaconda.org/fofanov)
[![Anaconda-Server Badge](https://anaconda.org/tara_furstenau/mtsv/badges/license.svg)](https://anaconda.org/tara_furstenau/mtsv)
[![Anaconda-Server Badge](https://anaconda.org/tara_furstenau/mtsv/badges/platforms.svg)](https://anaconda.org/tara_furstenau/mtsv)
# MTSv Pipeline
[Change log](CHANGELOG.md)

MTSv is a suite of metagenomic binning and analysis tools. It attempts to accurately identify which taxa are present in a given DNA sample by identifying *signature* reads that are unique to a single taxa. It assumes that read fragments in samples will be in a "shotgun" or short read format, typically ~50-200 bases in length.

## Pipeline Overview
The pipeline is broken into two major sections: (1) the download and setup of sequences from the sequence database which only needs to be completed once per database and as needed to include updates and (2) the main binning and analysis pipeline that is completed for each new read set. 

0. [Installation](#installation)
1. [Sequence Download and Setup Quick Start Guide](https://github.com/FofanovLab/MTSv/wiki/Sequence-Download-and-Setup-Quick-Start-Guide)
2. [Binning and Analysis Quick Start Guide](https://github.com/FofanovLab/MTSv/wiki/Binning-and-Analysis-Quick-Start-Guide)


## Installation
MTSv is written in Python but calls compiled rust binaries for most of the core functionality. It is currently available as a conda package (on OSX and linux-64 platforms) and we recommend that it is installed into a Python 3.6 conda environment.

### Dependencies
1. [Anaconda or Miniconda](https://conda.io/docs/user-guide/install/index.html)
2. Greater than GCC/5.2.0 on the PATH for runtime compilation of the C++11 code required for initial setup of the sequence databases.


### Create Conda environment with Python3.6.
```
$ conda create -n [ENV_NAME] python=3.6
```

#### Activate Conda Environment
```
$ source activate [ENV_NAME]
```
#### Install MTSv
```
$ conda install mtsv -c fofanov -c bioconda -c etetoolkit -c conda-forge
```
#### Deactivate Conda Environment
```
$ source deactivate
```
---
## License & Copyright

© 2018 FofanovLab

Licensed under the [MIT License](LICENSE)
