# cfCOP

**!!!cfCOP is currently in its beta development phase. While we are actively improving and testing the tool, please be aware that it may still have some bugs or incomplete features. We encourage you to use it and share feedback, but do so with caution in critical environments!!!**

cfCOP (cell-free Community observational Pipeline) is a computational workflow for pathogen identification from high-thoughput cell-free DNA samples. It was specifically developed to analyze cell-free DNA samples sequenced on Illumina and Oxford Nanopore platforms.


# Installation

Hardware prerequisites:
* 250 GB RAM
* multiple CPUs are advantageous, but not obligatory 
* 310 GB free hard drive storage (for the database)
* ~ 4x of the fastq.gz size (for analysis) 

Software prerequisites:
* Linux OS (cfCOP was tested on Ubuntu distributions)
* miniconda (see miniconda quick-installation manual [here](https://docs.anaconda.com/miniconda/install/#quick-command-line-install))
* git (installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git))

cfCOP is a Nextflow-based pipeline. All dependencies including Nextflow are provided in conda environment cfCOPenv.yml

```bash
conda env create -f environment.yml 
``` 

# Connecting a database

To retrieve a cfCOP metagenomic database use the command:

```bash
wget "https://www.dropbox.com/scl/fi/x194ided1ivqakouzz5ud/2024_12_07_cfCOPdb_export.tar.gz?rlkey=wap4nv3o1n7q0o624a8zqfhtm&e=1&dl=1" -O /your/path/to/2024_12_07_cfCOPdb_export.tar.gz
``` 

Unpack the database with:

```bash
tar -zxvf /your/path/to/2024_12_07_cfCOPdb_export.tar.gz
``` 

If cfCOP database is not connected, the tool will communicate that to you. To initially connect a database:

```bash
/your/path/to/run_cfCOP.py --connect_db /your/path/to/2024_12_07_export/
``` 

 
Once connected all paths the database are stored in the nextflow.config file. 

Please pay attention, that only minimap2 index is provided. BWA index can be created on your machine from a provided fasta file using command:

```bash
bwa index /your/path/to/combined_db.fasta
``` 

# Basic usage

cfCOP is run via invoking run_cfCOP.py script. 

run_cfCOP.py script requires either: 

* an input fastq-file(s) of sequenced sample(s)
* a sample sheet (please see the sample sheet format here)

Do not supply both sample sheet and additional samples, instead add your additional samples to a samplesheet and supply the samplesheet.

## Running cfCOP with input fastq file(s)

Both fastq and fastq.gz are possible.

With a single fastq file:

```bash
run_cfCOP.py --inputs /path/to/file.fastq 
``` 

With multiple fastq files:

```bash
run_cfCOP.py --inputs /path/to/files/file1.fastq /path/to/files/file2.fastq.gz
``` 

or with bash wildcards

```bash
run_cfCOP.py --inputs /path/to/files/file*.fastq.gz
``` 

Without an explicitly given mapper cfCOP uses a default mapper from nextflow.config file supplied with the pipeline. With --mapper option you set a desired mapper for the analysis. Available mapper choices: minimap2, bwa.  

```bash
run_cfCOP.py --inputs /path/to/file.fastq --mapper minimap2
``` 

Without an explicitly given sequencing platform cfCOP assumes data was produced on a default sequencing platform from nextflow.config file supplied with the pipeline. With --platform option you set a desired sequencing platform for the analysis. Available platform choices: illumina, nanopore. Choice of a sequencing platform influences the data mapping mode. 

```bash
run_cfCOP.py --inputs /path/to/file.fastq --platform illumina
``` 


## Running cfCOP with a samplesheet

```bash
run_cfCOP.py --samplesheet /path/to/samplesheet.txt 
``` 

### cfCOP Samplesheet format
To be continued...

# Additional parameters

# Outputs
## Read statistics (.readstat)
To be continued...
## Taxonomy statistics (.taxstat)
To be continued...
## Report
To be continued...