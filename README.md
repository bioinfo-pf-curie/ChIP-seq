# Chip-seq

**Institut Curie - Nextflow Chip-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.8-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with containers making installation trivial and results highly reproducible.
The current workflow is based on the nf-core best practice. See the nf-core project from details on [guidelines](https://nf-co.re/).


### Pipeline Summary

1. Run quality control of raw sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Align reads on reference genome ([`BWA`](http://bio-bwa.sourceforge.net/) / [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [`STAR`](https://github.com/alexdobin/STAR))
    * If using spike-in normalization, ambiguous reads and unmapped reads will be removed from both BAM files generated ([`pysam`](https://pysam.readthedocs.io/en/latest/api.html), [`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html))
3. Sort aligned reads ([`SAMTools`](http://www.htslib.org/))
4. Mark duplicates ([`Picard`](https://broadinstitute.github.io/picard/))
5. Library complexity analysis ([`Preseq`](http://smithlabresearch.org/software/preseq/))
6. Filtering aligned BAM files ([`SAMTools`](http://www.htslib.org/) & [`BAMTools`](https://github.com/pezmaster31/bamtools))
7. Computing Normalized and Relative Strand Cross-correlation (NSC/RSC) ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
8. Diverse alignment QCs and BigWig file creation ([`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html))
    * If using spike-in normalization, a scaling factor will be computed for BigWig generation ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
9. Peak calling for sharp & broad peaks ([`MACS2`](https://github.com/taoliu/MACS)) and very broad peaks ([`epic2`](https://github.com/biocore-ntnu/epic2))
10. Feature counting for every sample ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
11. Calculation of Irreproducible Discovery Rate in case of multiple replicates ([`IDR`](https://github.com/nboley/idr))
12. Peak annotation ([`HOMER`](http://homer.ucsd.edu/homer/ngs/annotation.html))
13. Results summary ([`MultiQC`](https://multiqc.info/))

### Quick help

```bash

```


### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda --replicates

```
If the analysis does not contain replicates, the --replicates option is not useful, but the test sample requires it.

#### Run the pipeline from a sample plan

```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --genome 'hg19' --outdir MY_OUTPUT_DIR

```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using the global environment (build by conda)
-profile condaPath

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda

```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs.
The sample plan is expected to be created as below :

SAMPLE_ID | SAMPLE_NAME | FASTQ_R1 [Path to R1.fastq file] | FASTQ_R2 [For paired end, path to Read 2 fastq]

### Design control

A design control is a csv file that list all experimental samples, their IDs, the associated input control, the replicate number and the expected peak type.
The design control is expected to be created as below :

SAMPLE_ID | CONTROL_ID | SAMPLE_NAME [Without '-ReplicateNumber'] | REPLICATE [Only the number] | PEAK_TYPE

Both files will be checked by the pipeline and have to be rigorously defined in order to make the pipeline work.

### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

<!-- Add developer names -->
This pipeline has been written by the bioinformatics platform of the Institut Curie 

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

