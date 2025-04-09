# cfCOP

**!!!cfCOP is currently in its beta development phase. While we are actively improving and testing the tool, please be aware that it may still have some bugs or incomplete features. We encourage you to use it and share feedback, but do so with caution in critical environments!!!**

cfCOP (cell-free Community observational Pipeline) is a computational workflow for pathogen identification from high-thoughput cell-free DNA samples. It was specifically developed to analyze cell-free DNA samples sequenced on Illumina and Oxford Nanopore platforms.


# Installation

Hardware prerequisites:
* 250 GB RAM
* multiple CPUs are advantageous, but not obligatory 
* 310 GB free hard drive storage (for the database)
* temporary hard drive storage of ca. 4x of the fastq.gz size (for analysis) 

Software prerequisites:
* Linux OS (cfCOP was tested on Ubuntu distributions)
* miniconda (see miniconda quick-installation manual [here](https://docs.anaconda.com/miniconda/install/#quick-command-line-install))
* git (installation instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git))

cfCOP is a Nextflow-based pipeline. All dependencies including Nextflow are provided in conda environment cfCOPenv.yml

```bash
cd /your/path/to/cfCOP
conda env create -f cfCOPenv.yml 
``` 

Before using cfCOP activate its environment:

```bash
conda activate cfCOPenv.yml 
``` 

# Connecting a database

To retrieve a cfCOP metagenomic database use the command:

```bash
wget "https://www.dropbox.com/scl/fi/asu9s784zb2iy7cam7w4k/2024_12_07_cfCOPdb_export.tar.gz?rlkey=edlenq7i34mylmws7lfc5grcq&e=1&st=sn0nnbyf&dl=1" -O /your/path/to/2024_12_07_cfCOPdb_export.tar.gz
``` 

Unpack the database with:

```bash
tar -zxvf /your/path/to/2024_12_07_cfCOPdb_export.tar.gz
``` 

If cfCOP database is not connected, the tool will communicate that to you. To initially connect a database:

```bash
/your/path/to/run_cfCOP.py --connect_db /your/path/to/2024_12_07_cfCOPdb_export/
``` 

 
Once connected all paths the database are stored in the nextflow.config file. 

Please pay attention, that only minimap2 index is provided. If you wish to use BWA mapper for your analysis (preferably for Illumina, and not for Nanopore data), BWA index can be created on your machine from a provided fasta file using command:

```bash
bwa index /your/path/to/combined_db.fasta
``` 

Afterwards please manually change the default path to parameter **bwa_index** in nextflow.config file.

# Basic usage

cfCOP is run via invoking run_cfCOP.py script. 

run_cfCOP.py script requires either: 

* an input fastq-file(s) of sequenced sample(s)
* a sample sheet (please see the sample sheet format here)

Do not supply both sample sheet and additional samples, instead add your additional samples to a samplesheet and supply the samplesheet.
Without further specifications cfCOP will run with default parameters specified in nextflow.config and/or in the samplesheet.

## Running cfCOP with input fastq file(s)

Both fastq and fastq.gz are possible.

With a single fastq file:

```bash
/your/path/to/run_cfCOP.py --inputs /path/to/file.fastq 
``` 

With multiple fastq files:

```bash
/your/path/to/run_cfCOP.py --inputs /path/to/files/file1.fastq /path/to/files/file2.fastq.gz
``` 

or with bash wildcards

```bash
/your/path/to/run_cfCOP.py --inputs /path/to/files/file*.fastq.gz
``` 

Without an explicitly given mapper cfCOP uses a default mapper from nextflow.config file supplied with the pipeline. With --mapper option you set a desired mapper for the analysis. Available mapper choices: minimap2, bwa.  

```bash
/your/path/to/run_cfCOP.py --inputs /path/to/file.fastq --mapper minimap2
``` 

Without an explicitly given sequencing platform cfCOP assumes data was produced on a default sequencing platform from nextflow.config file supplied with the pipeline. With --platform option you set a desired sequencing platform for the analysis. Available platform choices: illumina, nanopore. Choice of a sequencing platform influences the data mapping mode. 

```bash
/your/path/to/run_cfCOP.py --inputs /path/to/file.fastq --platform illumina
``` 


## Running cfCOP with a samplesheet

```bash
/your/path/to/run_cfCOP.py --samplesheet /your/path/to/samplesheet.txt 
``` 

Samplesheet overwrites the default mapper and sequencing platform for each sample. Samplesheet enables analyss of illumina and nanopore samples with both bwa and minimap2 mapper in one run.

Samplesheet is a tab-delimited file with linux line endings. Column order is pre-defined in header.

Columns of the samplesheet are:
> 1. **id**	*(mandatory)* - sample id. It should be unique to each sample
> 2. **platform** *(mandatory)* - sequencing platform. Choices: illumina, nanopore
> 3. **mapper** *(mandatory)* - mapper that should be used for the mapping step. Choices: minimap2, bwa
> 4. **fastq1** *(mandatory)* - absolute path to the fastq or fastq.gz file
> 5. **fastq2** *(optional)* - reversed reads for Illumina paired end samples
> 6. **sorted_sam** *(optional)* - samfile sorted by read name. The mapping step will be skipped then.

As an example, see provided template_samplesheet.txt

## Comand-Line arguments

| **Argument**           | **Description**                                                                                     |
|-------------------------|-----------------------------------------------------------------------------------------------------|
| `-h, --help`           | Show this help message and exit.                                                                    |
| `--samplesheet SAMPLESHEET` | Path to the sample sheet.                                                   |
| `--inputs INPUTS [INPUTS ...]` | List of input files (for multiple files).                                                    |
| `--mapper {minimap2,bwa}` | Specify the aligner to use. Choices are: `minimap2` or `bwa`.                                      |
| `--platform {nanopore,illumina}` | Specify the sequencing platform of sample(s). Choices are: `nanopore` or `illumina`.       |
| `--kmer_cutoff KMER_CUTOFF` | Report generation: Kmer index cutoff should be a float number between 0 and 1.                  |
| `--count_cutoff COUNT_CUTOFF` | Report generation: Number of reads assigned to a taxon in order to be reported. Should be an integer. |
| `--fraction_cutoff FRACTION_CUTOFF` | Report generation: Fraction of reads (of total reads) assigned to a taxon in order to be reported. Should be a float number between 0 and 1. |
| `--connect_db CONNECT_DB` | Connect to a database by filling the paths of the database files into the `nextflow.config`.       |
| `--resume`             | Resume the pipeline with changed parameters or for debugging purposes. Cached unchanged processes will not be re-run. |
| `--ro_report`             | Per default, cfCOP creates an html report of the current run. To supress reporting, set this flag. |

# Advanced configuration

cfCOP can be further configured via nextflow.config file. Changing this file changes cfCOP default parameters. Below is a list of parameters that can be configured for the pipeline.


### General Parameters
- **`bwa_index`**  
  Path to the BWA index for alignment.

- **`minimap2_index`**  
  Path to the Minimap2 index for alignment.

- **`database_fasta`**  
  Path to the FASTA file for the database.

- **`taxonomy`**  
  Path to the taxonomy database folder.


### Filters
- **`complexity_filter`**  
  Path to the complexity filter file.

- **`homology_filter`**  
  Path to the homology filter file.

- **`plasmid_filter`**  
  Path to the plasmid filter file.


### Resource Parameters
- **`threads`**  
  Number of threads for parallel processing.

- **`threads_samtools`**  
  Number of threads for Samtools processing.

- **`kmer_size`**  
  Size of the k-mer used in analysis.


### Filtering and Reporting
- **`taxrank`**  
  Taxonomic ranks to consider in reporting. Options: `["species", "genus", "family"]`.

- **`platform`**  
  Sequencing platform used. Options: `"illumina"`, `"nanopore"`.

- **`mapper`**  
  Alignment tool to be used. Options: `"minimap2"`, `"bwa"`.

- **`kmer_cutoff`**  
  K-mer index cutoff for reporting. Must be a float between 0 and 1.

- **`count_cutoff`**  
  Minimum number of reads assigned to a taxon for reporting.

- **`fraction_cutoff`**  
  Minimum fraction of reads (relative to total reads) assigned to a taxon for reporting. Must be a float between 0 and 1.


### Input Parameters
- **`input_files`**  
  Path to input files. Can be null or a list of paths.

- **`samplesheet`**  
  Path to the sample sheet. Can be null or a valid file path.


### Executor Settings
- **`$sge.queueSize`**  
  Maximum number of queued jobs for SGE executor.

- **`$local.cpus`**  
  Number of CPUs for local execution.

- **`$local.memory`**  
  Memory allocation for local execution.


### Profiles
#### **`standard`**  
Uses the local executor.  
Configuration:  

process.executor = 'local'

#### **`cluster `**
Execution on a cluster (using PBSPro).
To enable the cluster profile, use the -profile cluster option when running the pipeline.
Adjust the cluster options (params.project_cluster_opts) as needed for your environment.

Executor settings:
- **`process.executor`**  
  Specifies the executor to use. For the cluster profile, this is set to `'pbspro'`.

- **`params.project_cluster_opts`**  
  Cluster-specific options passed to the job scheduler. 

For the cluster execution job are submitted depending on the memory and ime consumption.

**Big Memory Jobs**
- **`cpus`**: `params.threads`  
  Number of CPUs allocated for big memory jobs.
- **`memory`**: `300 GB`  
  Memory allocated for big memory jobs.
- **`time`**: `1h`  
  Maximum runtime for big memory jobs.

**Low Memory Jobs**
- **`cpus`**: `1`  
  Number of CPUs allocated for low memory jobs.
- **`memory`**: `10 GB`  
  Memory allocated for low memory jobs.
- **`time`**: `5h`  
  Maximum runtime for low memory jobs.
 

# Outputs
## Read statistics (.readstat)

Read statistics include mapping and filtering statistics for each read. A *readstat* file consists of one row per read. Read format includes 8 mandatory fields:

1. **`#flag`** - Status flag indicating the processing or filtering status of the read.  
2. **`name`** - Unique identifier or name of the read.  
3. **`taxid`** - Taxonomy ID assigned to the read based on mapping results.  
4. **`length`** - Length of the read in base pairs (bp).  
5. **`stats`** - General mapping statistics for the read, such as alignment quality or coverage.  
6. **`complexity`** - Indicates whether the read was filtered due to low sequence complexity *(e.g., 0 for pass, 1 for fail)*.  
7. **`homology`** - Indicates whether the read was filtered due to homology to undesired sequences *(e.g., host DNA)*.  
8. **`plasmid`** - Indicates whether the read was identified as originating from a plasmid.  

**Example:**

| **#flag** | **name**  | **taxid** | **length** | **stats** | **complexity** | **homology** | **plasmid** |
|-----------|-----------|-----------|------------|-----------|----------------|--------------|-------------|
| 0         | read1     | 562       | 150        | 95.2      | 0              | 0            | 0           |
| 1         | read2     | 9606      | 200        | 89.4      | 1              | 0            | 0           |
| 0         | read3     | 1234      | 120        | 92.8      | 0              | 1            | 1           |

- **`read1`**: Passed all filters and is assigned to `taxid` 562.  
- **`read2`**: Failed the complexity filter.  
- **`read3`**: Passed the complexity filter but failed the homology and plasmid filters.  

## Taxonomy statistics (.taxstat)

Taxonomy statistics summarized the read statistics on a taxon level.

## Report

The compact overview of the detected taxons on family, genus and species levels