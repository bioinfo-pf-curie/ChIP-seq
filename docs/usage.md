# Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`conda`](#conda)
        * [`toolsPath`](#toolsPath)
        * [`singularity`](#singularity)
        * [`cluster`](#cluster)
        * [`test`](#test)
    * [`--reads`](#--reads)
	* [`--samplePlan`](#--samplePlan)
	* [`--design`](#--design)
	* [`--singleEnd`](#--singleend)
* [Mapping](#mapping)
	* [`--aligner`](#--aligner)
	* [`--saveAlignedIntermediates`](#--saveAlignedIntermediates)
* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
	* [`--spike`](#--spike)
* [Filtering](#fitlering)
    * [`--mapQ`](#--mapq)
	* [`--keepDups`](#--keepDups)
	* [`--blacklist`](#--blacklist)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
	* [`--metadata`](#--metadta)
	* [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--multiqc_config`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile 'singularity'
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outdir/-w` options.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile singularity` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `toolsPath`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Use the conda images available on the cluster
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Use the singularity images available on the cluster
* `cluster`
    * Run the workflow on the computational cluster
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--samplePlan`

Use this to specify a sample plan file instead of a regular expression to find fastq files. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv
```

The sample plan is a csv file with the following information :

Sample ID | Sample Name | Path to R1 fastq file | Path to R2 fastq file

### `--design`

Specify a design for extended analysis which require metadata (peak calling, IDR, etc.).

The expected format is the following :

Sample ID | Control ID | Sample Name | Group | Peak Type

Note that the `Sample ID` and ` Control ID` must match the IDs of the samplePlan.

The `Group` information allows to consider as replicates all samples from the same group.

The peak type is used to infer the algorithm to use for the peak calling. Expected values are `sharp`, `broad` or `very-broad`.

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

## Mapping

### `--aligner`

The current version of the pipeline supports three different aligners;
- [`star`](https://github.com/alexdobin/STAR). Default value.
- [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [`bwa-mem`](http://bio-bwa.sourceforge.net/bwa.shtml)

By default, the `BWA-mem` mapper is run. You can specify the tool to use as follows:

```bash
--aligner 'bwa-mem'
```

### `--saveAlignedIntermediates`

By default, only the final bam files are saved in the results folder.
Activate this option, if you want to keep all intermediate bam files. Note that it will take a lot of disk space.

## Reference genomes

The pipeline config files come bundled with paths to the genomes reference files. 

### `--genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg38`
* Mouse
  * `--genome mm10`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the genomes resource. 
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'hg19' {
       fasta = path to genome fasta file
	   chrsize = path to chromosome size file
       bwaIndex = path to bwa-mem index
       bowtie2Index = path to bowtie index
       starIndex = path to star index
       geneBed = path to BED gene file
       gtf = path to GTF annotation file
       effGenomeSize = effective genome size value
       blacklist = path to blacklist regions (.bed) file
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

Note that these paths can be updated on command line using the following parameters:
- `--fasta` - Path to genome fasta file
- `--starIndex` - Path to STAR index
- `--bwaIndex` - Path to Bwa-mem index
- `--bowtie2Index` - Path to Bowtie2 index
- `--gtf` - Path to GTF file
- `--geneBed` - Path to gene bed file
- `--effGenomeSize` - Effective genome size. Mandatory for peak calling and deeptools usage.
- `--blacklist` - Path to black list genome file defined by ENCODE

### `--spike`

Same as `--genome`. Define the organism of spike-in which has to be used.

In this case, the reads are mapped on both the reference and the spike genomes.

### `--tssSize`

Define the size of the upstream/downstream region of the transcription start sites should be consider as promoter region.
Default: 2000 bp.

## Filtering

### `--mapQ`

Keep reads with a mapping quality value (mapQ) higher that `--mapQ`.
By default, no filtering is performed.

### `--keepDups`

Keep reads flagged as duplicates for downstream analysis. By default, these reads are discarded.

### `--blacklist`

Remove reads from blacklist regions at the peak calling level and in the deeptools analysis.
See the `conf/genomes.config` for details.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/base.conf` file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:
- `--skip_qc` - Skip all QC steps apart from MultiQC
- `--skip_rrna` - Skip rRNA mapping
- `--skip_fastqc` - Skip FastQC
- '--skip_genebody_coverage' - Skip genebody coverage step
- `--skip_saturation` - Skip Saturation qc
- `--skip_dupradar` - Skip dupRadar (and Picard MarkDups)
- `--skip_readdist` - Skip read distribution steps
- `--skip_expan` - Skip exploratory analysis
- `--skip_multiqc` - Skip MultiQC
				
### `--metadata`

Specify a two-columns (tab-delimited) metadata file to diplay in the final Multiqc report.

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
