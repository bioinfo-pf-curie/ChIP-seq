# ChIP-seq

**Institut Curie - Nextflow ChIP-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.11-blue.svg)](https://multiqc.info/)
[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

[![DOI](https://zenodo.org/badge/269415034.svg)](https://zenodo.org/badge/latestdoi/269415034)

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with containers making installation trivial and results highly reproducible.
The current workflow was initiated from the [nf-core ChIP-seq pipeline](https://github.com/nf-core/chipseq). See the nf-core project from details on [guidelines](https://nf-co.re/).

### Pipeline Summary

1. Trim adapters from sequencing reads ([`TrimGalore!`](https://github.com/FelixKrueger/TrimGalore)
2. Run quality control of raw sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align reads on reference genome ([`BWA`](http://bio-bwa.sourceforge.net/) / [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [`STAR`](https://github.com/alexdobin/STAR))
    * If spike-in are used, mapping on spike genome is run and ambiguous reads are removed from both BAM files ([`pysam`](https://pysam.readthedocs.io/en/latest/api.html))
4. Sort aligned reads ([`SAMTools`](http://www.htslib.org/))
5. Mark duplicates ([`Picard`](https://broadinstitute.github.io/picard/))
6. Library complexity analysis ([`Preseq`](http://smithlabresearch.org/software/preseq/))
7. Filtering aligned BAM files ([`SAMTools`](http://www.htslib.org/) & [`BAMTools`](https://github.com/pezmaster31/bamtools))
   - reads mapped to blacklisted regions
   - reads marked as duplicates
   - reads that arent marked as primary alignments
   - reads that are unmapped
   - reads mapped with a low mapping quality (multiple hits, secondary alignments, etc.)
8. Computing Normalized and Relative Strand Cross-correlation (NSC/RSC) ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
9. Diverse alignment QCs and bigWig file creation ([`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html))
    * If spike-in are used, a scaling factor is computed and additional bigWig are generated ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
10. Peak calling for sharp, broad peaks and very-broad peaks ([`MACS2`](https://github.com/taoliu/MACS)) and very broad peaks ([`epic2`](https://github.com/biocore-ntnu/epic2))
11. Feature counting for every sample at gene and transcription start sites loci ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
12. Calculation of Irreproducible Discovery Rate in case of multiple replicates ([`IDR`](https://github.com/nboley/idr))
13. Peak annotation and QC ([`HOMER`](http://homer.ucsd.edu/homer/ngs/annotation.html))
14. Results summary ([`MultiQC`](https://multiqc.info/))

### Quick help

```bash

N E X T F L O W  ~  version 21.10.6
Launching `main.nf` [tender_stallman] - revision: fcda6ad7de
------------------------------------------------------------------------

    _____ _    ___________ 
   /  __ \ |  |_   _| ___ \
   | /  \/ |__  | | | |_/ /_____ ___  ___  __ _ 
   | |   | '_ \ | | |  __/______/ __|/ _ \/ _` |
   | \__/\ | | || |_| |         \__ \  __/ (_| |
    \____/_| |_\___/\_|         |___/\___|\__, |
                                             | |
                                             |_|
																								   
                   v2.0.0
------------------------------------------------------------------------
																																  
Usage:
	
The typical command for running the pipeline is as follows:
		
nextflow run main.nf --reads PATH --samplePlan PATH --genome STRING -profile PROFILES
			

MANDATORY ARGUMENTS:
  --genome     STRING   Name of the reference genome.
  --reads      PATH     Path to input data (must be surrounded with quotes)
  --samplePlan PATH     Path to sample plan (csv format) with raw reads (if `--reads` is not specified)
			
INPUTS:
  --bam                    For aligned (BAM) input data
  --design       PATH      Path to design file (csv format)
  --fragmentSize INTEGER   Estimated fragment length used to extend single-end reads
  --singleEnd              For single-end input data
  --spike        INTEGER   Name of the genome used for spike-in analysis
  --tssSize      INTEGER   Distance (upstream/downstream) to transcription start point to consider
									
PREPROCESSING:
  --trimming           Trim adapters with TrimGalore
  
REFERENCES:
  --effectiveGenomeSize  INTEGER   Effective genome size
  --fasta                PATH      Path to genome fasta file
  --geneBed              PATH      Path to gene file (BED)
  --genomeAnnotationPath PATH      Path to genome annotations folder
  --gtf                  PATH      Path to GTF annotation file. Used in HOMER peak annotation
  --spikeFasta           PATH      Path to Fasta reference for spike-in
				
ALIGNMENT:
  --spikeBwaIndex     PATH                             Spike-in Index for Bwa-mem aligner
  --aligner           STRING [bwa-mem, star, bowtie2]  Tool for reads alignment
  --bowtie2Index      PATH                             Indexes for Bowtie2 aligner
  --bwaIndex          PATH                             Path to Bwa-mem indexes
  --spikeBowtie2Index PATH                             Spike-in indexes for Bowtie2 aligner
  --spikeStarIndex    PATH                             Path to STAR indexes of spike-in reference
  --starIndex         PATH                             Indexes for STAR aligner
  
OTHER OPTIONS:
  --metadata          PATH      Specify a custom metadata file for MultiQC
  --multiqcConfig     PATH      Specify a custom config file for MultiQC
  --name              STRING    Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
  --outDir            PATH      The output directory where the results will be saved
  --saveIntermediates           Save intermediates files
					
ANALYSIS:
  --noReadExtension           Do not extend reads to fragment length
					
FILTERING:
  --blacklist          PATH      Path to black list regions (.bed). See the genome.config for details
  --keepDups                     Do not remove duplicates afer marking
  --keepSingleton                Keep unpaired reads
  --mapq               INTEGER   Minimum mapping quality to consider
  --spikePercentFilter INTEGER   Minimum percent of reads aligned to spike-in genome
  
SKIP OPTIONS:
  --skipCounts                Disable counts analysis
  --skipDeeptools             Disable deeptools QC
  --skipFastqc                Disable Fastqc
  --skipIDR                   Disable IDR analysis
  --skipMultiqc               Disable MultiQC
  --skipPPQT                  Disable phantompeakqualtools QC
  --skipPeakCalling           Disable peak calling analysis
  --skipPeakanno              Disable peaks annotation
  --skipSaturation            Disable saturation analysis with Preseq
									
=======================================================
Available Profiles
  -profile test                        Run the test dataset
  -profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
  -profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
  -profile docker                      Use the Docker images for each process
  -profile singularity                 Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
  -profile cluster                     Run the workflow on the cluster, instead of locally
															
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda

```

#### Run the pipeline from a `sample plan` and a `design` file
```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --design MY_DESIGN --genome 'hg19' --genomeAnnotationPath ANNOTATION_PATH --outDir MY_OUTPUT_DIR

```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
-profile path --globalPath INSTALLATION_PATH

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityPath SINGULARITY_PATH

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE

```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs.
The sample plan is expected to be created as below :

SAMPLE_ID,SAMPLE_NAME,FASTQ_R1 [Path to R1.fastq file],FASTQ_R2 [For paired end, path to Read 2 fastq]

Note that if you already have aligned data, you can specified BAM file in the sample plan and add the `--bam` to specify that inputs are bam files.

SAMPLE_ID,SAMPLE_NAME,BAM [Path to bam file]

### Design control

A design control is a csv file that list all experimental samples, their IDs, the associated input control (or IgG), the replicate number and the expected peak type.
The design control is expected to be created as below :

SAMPLE_ID,CONTROL_ID,GROUP,PEAKTYPE

Both files will be checked by the pipeline and have to be rigorously defined in order to make the pipeline work.  
Note that the control is optional if not available but is highly recommanded.  
If the `design` file is not specified, the pipeline will run until the alignment, QCs and track generation. The peak calling and the annotation will be skipped.


### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (Valentin Laroche, Nicolas Servant)

#### Citation

If you use this pipeline for your project, please cite it using the following doi: [10.5281/zenodo.7443721](https://doi.org/10.5281/zenodo.7443721)
Do not hesitate to use the Zenodo doi corresponding to the version you used !

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

