#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
            ChIP-seq
========================================================================================
ChIP-seq Analysis Pipeline.
#### Homepage / Documentation
https://gitlab.curie.fr/data-analysis/chip-seq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/devMessage.txt")
  log.info devMess.text
  }

  log.info"""

  ChIP-seq v${workflow.manifest.version}
  ======================================================================

  Usage:

  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19' 
  nextflow run main.nf --samplePlan 'sample_plan.csv' --design 'design.csv' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19'

  Mandatory arguments:
  --reads [file]                     Path to input data (must be surrounded with quotes)
  --samplePlan [file]                Path to sample plan file if '--reads' is not specified
  --genome [str]                     Name of genome reference. See the `--genomeAnnotationPath` to defined the annotations path.
  -profile [str]                     Configuration profile to use. Can use multiple (comma separated)

  Inputs:
  --design [file]                    Path to design file for downstream analysis
  --singleEnd [bool]                 Specifies that the input is single end reads
  --fragmentSize [int]               Estimated fragment length used to extend single-end reads. Default: 200
  --spike [str]                      Name of the genome used for spike-in analysis. Default: false

  References: If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --genomeAnnotationPath [file]      Path  to genome annotation folder
  --fasta [file]                     Path to Fasta reference
  --spikeFasta [file]                Path to Fasta reference for spike-in

  Alignment:
  --aligner [str]                    Alignment tool to use ['bwa-mem', 'star', 'bowtie2']. Default: 'bwa-mem'
  --saveAlignedIntermediates [bool]  Save all intermediates mapping files. Default: false  
  --starIndex [dir]                  Index for STAR aligner
  --spikeStarIndex [dir]             Spike-in Index for STAR aligner
  --bwaIndex [file]                  Index for Bwa-mem aligner
  --spikeBwaIndex [file]             Spike-in Index for Bwa-mem aligner
  --bowtie2Index [file]              Index for Bowtie2 aligner
  --spikeBowtie2Index [file]         Spike-in Index for Bowtie2 aligner

  Filtering:
  --mapq [int]                       Minimum mapping quality to consider. Default: 10
  --keepDups [bool]                  Do not remove duplicates afer marking. Default: false
  --keepSingleton [bool]             Keep unpaired reads. Default: false
  --blacklist [file]                 Path to black list regions (.bed).
  --spikePercentFilter [float]       Minimum percent of reads aligned to spike-in genome. Default: 0.2

  Analysis:
  --noReadExtension [bool]           Do not extend reads to fragment length. Default: false

  Annotation:          If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --genomeAnnotationPath [dir]       Path to genome annotations.
  --geneBed [file]                   BED annotation file with gene coordinate.
  --gtf [file]                       GTF annotation file. Used in HOMER peak annotation
  --effGenomeSize [int]              Effective Genome size
  --tssSize [int]                    Distance (upstream/downstream) to transcription start point to consider. Default: 2000

  Skip options:        All are false by default
  --skipFastqc [bool]                Skips fastQC
  --skipPreseq [bool]                Skips preseq QC
  --skipPPQT [bool]                  Skips phantompeakqualtools QC
  --skipDeepTools [bool]             Skips deeptools QC
  --skipPeakcalling [bool]           Skips peak calling
  --skipPeakanno [bool]              Skips peak annotation
  --skipIDR [bool]                   Skips IDR QC
  --skipFeatCounts [bool]            Skips feature count
  --skipMultiQC [bool]               Skips MultiQC step

  Other options:
  --metadata [file]                  Path to metadata file for MultiQC report
  --outDir [dir]                     The output directory where the results will be saved
  -w/--work-dir [dir]                The temporary directory where intermediate data will be saved
  -name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

  =======================================================
  Available Profiles
    -profile test                    Run the test dataset
    -profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
    -profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
    -profile path                    Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
    -profile multipath               Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path 
    -profile docker                  Use the Docker images for each process
    -profile singularity             Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
    -profile cluster                 Run the workflow on the cluster, instead of locally

  """.stripIndent()
}

// Show help messsage
if (params.help){
  helpMessage()
  exit 0
}

if (params.aligner != 'bwa-mem' && params.aligner != 'star' && params.aligner != 'bowtie2' ) {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'bowtie2' or 'bwa-mem'"
}

/*********************
 * Fasta file
 */

// Configurable reference genomes
genomeRef = params.genome

// Genome Fasta file
params.fasta = genomeRef ? params.genomes[ genomeRef ].fasta ?: false : false
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set{chFastaHomer}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

// Chromosome size file
params.chrsize = genomeRef ? params.genomes[ genomeRef ].chrsize ?: false : false
if ( params.chrsize ){
  Channel
    .fromPath(params.chrsize, checkIfExists: true)
    .set{chChromSize}
}
else{
  exit 1, "Chromosome size file not found: ${params.chrsize}"
}

// spike
if (params.spike || (params.spikeFasta && (params.spikeBwaIndex || params.spikeBt2Index || params.spikeStarIndex))){
  useSpike = true
}else{
  useSpike = false
}

params.spikeFasta = params.spike ? params.genomes[ params.spike ].fasta ?: false : false
if ( params.spikeFasta ){
  Channel
    .fromPath(params.spikeFasta, checkIfExists: true)
    .set{chFastaSpike}
}else{
  chFastaSpike = Channel.empty()
}

/********************
 * Bwa-mem Index
 */

params.bwaIndex = genomeRef ? params.genomes[ genomeRef ].bwaIndex ?: false : false
if (params.bwaIndex){
  lastPath = params.bwaIndex.lastIndexOf(File.separator)
  bwaDir =  params.bwaIndex.substring(0,lastPath+1)
  bwaBase = params.bwaIndex.substring(lastPath+1)
  Channel
    .fromPath(bwaDir, checkIfExists: true)
    .ifEmpty {exit 1, "BWA index file not found: ${params.bwaIndex}"}
    .combine( [ bwaBase ] )
    .set { chBwaIndex }
} else {
  exit 1, "BWA index file not found: ${params.bwaIndex}"
}

params.spikeBwaIndex = params.spike ? params.genomes[ params.spike ].bwaIndex ?: false : false
if (params.spikeBwaIndex){
  lastPath = params.spikeBwaIndex.lastIndexOf(File.separator)
  bwaDirSpike =  params.spikeBwaIndex.substring(0,lastPath+1)
  spikeBwaBase = params.spikeBwaIndex.substring(lastPath+1)
  Channel
    .fromPath(bwaDirSpike, checkIfExists: true)
    .ifEmpty {exit 1, "Spike BWA index file not found: ${params.spikeBwaIndex}"}
    .combine( [ spikeBwaBase ] ) 
    .set { chSpikeBwaIndex }
  
  chBwaIndex = chBwaIndex.concat(chSpikeBwaIndex)
}

/*********************
 * Bowtie2 indexes
 */

params.bt2Index = genomeRef ? params.genomes[ genomeRef ].bowtie2Index ?: false : false
if (params.bt2Index){
  lastPath = params.bt2Index.lastIndexOf(File.separator)
  bt2Dir =  params.bt2Index.substring(0,lastPath+1)
  bt2Base = params.bt2Index.substring(lastPath+1)
  Channel
    .fromPath(bt2Dir, checkIfExists: true)
    .ifEmpty {exit 1, "Bowtie2 index file not found: ${params.bt2Index}"}
    .combine( [ bt2Base ] ) 
    .set { chBt2Index }
} else {
  exit 1, "Bowtie2 index file not found: ${params.bt2Index}"
}

params.spikeBt2Index = params.spike ? params.genomes[ params.spike ].bowtie2Index ?: false : false
if (params.spikeBt2Index){
  lastPath = params.spikeBt2Index.lastIndexOf(File.separator)
  bt2DirSpike =  params.spikeBt2Index.substring(0,lastPath+1)
  spikeBt2Base = params.spikeBt2Index.substring(lastPath+1)
  Channel
    .fromPath(bt2DirSpike, checkIfExists: true)
    .ifEmpty {exit 1, "Spike Bowtie2 index file not found: ${params.spikeBt2Index}"}
    .combine( [ spikeBt2Base ] ) 
    .set { chSpikeBt2Index }

  chBt2Index = chBt2Index.concat(chSpikeBt2Index)
}

/********************
 * STAR indexes
 */

params.starIndex = genomeRef ? params.genomes[ genomeRef ].starIndex ?: false : false
if (params.starIndex){
  Channel
    .fromPath(params.starIndex, checkIfExists: true)
    .ifEmpty {exit 1, "STAR index file not found: ${params.starIndex}"}
    .combine( [ genomeRef ] ) 
    .set { chStarIndex }
} else {
  exit 1, "STAR index file not found: ${params.starIndex}"
}

params.spikeStarIndex = params.spike ? params.genomes[ params.spike ].starIndex ?: false : false
if (params.spikeStarIndex){
  if (params.spike){
    genomeSpike = params.spike
  }else{
    genomeSpike = 'spikeGenome'
  }

  Channel
    .fromPath(params.spikeStarIndex, checkIfExists: true)
    .ifEmpty {exit 1, "Spike STAR index file not found: ${params.spikeStarIndex}"}
    .combine( [ genomeSpike ] )
    .set { chSpikeStarIndex }

  chStarIndex = chStarIndex.concat(chSpikeStarIndex)
}

/*********************
 * Annotations
 */

params.gtf = genomeRef ? params.genomes[ genomeRef ].gtf ?: false : false
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .into{chGtfHomer; chGtfFeatCounts}
}
else {
  exit 1, "GTF annotation file not specified!"
}

params.geneBed = genomeRef ? params.genomes[ genomeRef ].geneBed ?: false : false
if (params.geneBed) {
  Channel
    .fromPath(params.geneBed, checkIfExists: true)
    .ifEmpty {exit 1, "BED file ${geneBed} not found"}
    .into{chGeneBedDeeptools; chGenePrepareAnnot; chGeneFeatCounts}
}else{
  chGeneBedDeeptools = Channel.empty()
  chGenePrepareAnnot = Channel.empty()
  chGeneFeatCounts = Channel.empty()
}

params.blacklist = genomeRef ? params.genomes[ genomeRef ].blacklist ?: false : false
if (params.blacklist) { 
  Channel
    .fromPath(params.blacklist, checkIfExists: true)
    .into {chBlacklistBigWig; chBlacklistBigWigSpike; chBlacklistCorrelation} 
}else{
  chBlacklistBigWig = Channel.empty()
  chBlacklistBigWigSpike = Channel.empty()
  chBlacklistCorrelation = Channel.empty()
}


/***********************
 * Header and conf
 */

//PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

//Peak Calling
params.effGenomeSize = genomeRef ? params.genomes[ genomeRef ].effGenomeSize ?: false : false
if (!params.effGenomeSize) {
  log.warn "=================================================================\n" +
            "  WARNING! Effective Genome Size is not defined.\n" +
            "  Peak calling and annotation will be skipped.\n" +
            "  Please specify value for '--effGenomeSize' to run these steps.\n" +
            "================================================================"
}

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .into { chPeakCountHeaderSharp; chPeakCountHeaderBroad; chPeakCountHeaderVeryBroad }

Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .into{ chFripScoreHeaderSharp; chFripScoreHeaderBroad; chFripScoreHeaderVeryBroad }

Channel
  .fromPath("$baseDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }

//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
chOutputDocs = file("$baseDir/docs/output.md", checkIfExists: true)
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

//Has the run name been specified by the user?
//This has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

//Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/devMessage.txt")
  log.info devMess.text
}

log.info """=======================================================

ChIP-seq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'ChIP-seq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = customRunName ?: workflow.runName
if (params.samplePlan) {
  summary['SamplePlan'] = params.samplePlan
} else{
  summary['Reads']      = params.reads
}
summary['Design']       = params.design ?: "None"
summary['Annotation']   = params.genomeAnnotationPath
summary['Fasta Ref']    = params.fasta
summary['Spikes'] = params.spike ? "${params.spike}" : useSpike ? "Yes" : "False"
if (params.spikeFasta)  summary["Fasta spike"] = params.spikeFasta
summary['GTF']          = params.gtf
summary['Genes']        = params.geneBed
if (params.blacklist)  summary['Blacklist '] = params.blacklist
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
if (params.singleEnd)  summary['Fragment Size '] = params.fragmentSize
summary['Aligner'] = params.aligner
if (params.keepDups)  summary['Keep Duplicates'] = 'Yes'
if (params.mapq)  summary['Min MapQ'] = params.mapq
summary['Max Memory']   = params.maxMemory
summary['Max CPUs']     = params.maxCpus
summary['Max Time']     = params.maxTime
summary['Output dir']   = params.outDir
summary['Working dir']  = workflow.workDir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * CHANNELS
 */

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}else{
  chMetadata=Channel.empty()
}                                                                                                                                                                                           
 

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
  if(params.singleEnd && !params.inputBam){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
  }else if (!params.singleEnd && !params.inputBam){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]]}
      .set { chAlignReads }
   params.reads=false
  }
} else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
   } else {
     Channel
       .from(params.readPaths)
       .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
       .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
       .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
   }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
}


/**************************
 * Make sample plan if not available
 */

if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .into {chSplan; chSplanCheck}
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .into{ chSplan; chSplanCheck }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .into{ chSplan; chSplanCheck }
  }
} else if(params.bamPaths){
  Channel
    .from(params.bamPaths)
    .collectFile() {
      item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
     }
    .into{ chSplan; chSplanCheck }
  params.aligner = false
} else {
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
      }     
      .into { chSplan; chSplanCheck }
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }     
      .into { chSplan; chSplanCheck }
   }
}
/******************************
 * Design file
 */

if (!params.design) {
  log.info "=================================================================\n" +
            "  INFO: No design file detected.\n" +
            "  Peak calling and annotation will be skipped.\n" +
            "  Please set up a design file '--design' to run these steps.\n" +
            "================================================================"
}

if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { chDesignCheck; chDesignControl; chDesignMqc }

  chDesignControl
    .splitCsv(header:true)
    .map { row ->
      if(row.CONTROLID==""){row.CONTROLID='NO_INPUT'}
      return [ row.SAMPLEID, row.CONTROLID, row.SAMPLENAME, row.GROUP, row.PEAKTYPE ]
     }
    .set { chDesignControl }

  // Create special channel to deal with no input cases
  Channel
    .from( ["NO_INPUT", ["NO_FILE","NO_FILE"]] )
    .toList()
    .set{ chNoInput }
}else{
  chDesignControl = Channel.empty()
  chDesignCheck = Channel.empty()
  chDesignMqc = Channel.empty()
}

/*
 * Check design and sampleplan files
 */

process checkDesign{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.summaryDir}/", mode: 'copy'

  when:
  params.design

  input:
  file design from chDesignCheck
  file samplePlan from chSplanCheck

  script:
  optSE = params.singleEnd ? "--singleEnd" : ""
  """
  checkDesign.py -d $design -s $samplePlan ${optSE}
  """
}

/*
 * FastQC
 */

process fastQC{
  tag "${prefix}"
  label 'fastqc'
  label 'lowCpu'
  label 'minMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy'

  when:
  !params.skipFastqc && !params.inputBam

  input:
  set val(prefix), file(reads) from rawReadsFastqc

  output:
  file "*_fastqc.{zip,html}" into chFastqcMqc
  file("v_fastqc.txt") into chFastqcVersion

  script:
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads -t ${task.cpus}
  """
}

/*
 * Alignment on reference genome
 */

/* BWA-MEM */
process bwaMem{
  tag "${sample} on ${genomeBase}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename -> 
	     if (filename.indexOf(".log") > 0) "logs/$filename" 
	     else if (params.saveAlignedIntermediates) filename}

  when:
  params.aligner == "bwa-mem" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase) from rawReadsBWA.combine(chBwaIndex)

  output:
  set val(sample), file("*.bam") into chAlignReadsBwa
  file("*.log") into chBwaMqc
  file("v_bwa.txt") into chBwaVersion

  script:
  prefix = genomeBase == genomeRef ? sample : sample + '_spike'
  opts = params.bwaOpts
  """
  echo \$(bwa 2>&1) &> v_bwa.txt
  bwa mem -t ${task.cpus} \
           ${index}/${genomeBase} \
          ${opts} \
          $reads | samtools view -bS - > ${prefix}.bam
  getBWAstats.sh -i ${prefix}.bam -p ${task.cpus} > ${prefix}_bwa.log
  """
}

/* BOWTIE2 */
process bowtie2{
  tag "${sample} on ${genomeBase}"
  label 'bowtie2'
  label 'highCpu' 
  label 'highMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
              saveAs: {filename ->
	      if (filename.indexOf(".log") > 0) "logs/$filename"  
	      else if (params.saveAlignedIntermediates) filename}
  when:
  params.aligner == "bowtie2" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase) from rawReadsBt2.combine(chBt2Index)

  output:
  set val(sample), file("*.bam") into chAlignReadsBowtie2
  file("*.log") into chBowtie2Mqc
  file("v_bowtie2.txt") into chBowtie2Version

  script:
  prefix = genomeBase == genomeRef ? sample : sample + '_spike'
  readCommand = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  opts = params.bowtie2Opts
  """
  bowtie2 --version &> v_bowtie2.txt
  bowtie2 -p ${task.cpus} \
          ${opts} \
           -x ${index}/${genomeBase} \
          $readCommand > ${prefix}.bam 2> ${prefix}_bowtie2.log
  """
}

/* STAR */
process star{
  tag "${sample} on ${genomeBase}"
  label 'star'
  label 'highCpu'
  label 'extraMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename ->
	     if (filename.indexOf(".log") > 0) "logs/$filename"  
 	     else if (params.saveAlignedIntermediates) filename}
  when:
  params.aligner == "star" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase) from rawReadsSTAR.combine(chStarIndex)

  output:
  set val(sample), file('*.bam') into chAlignReadsStar
  file ("*Log.final.out") into chStarMqc
  file("v_star.txt") into chStarVersion

  script:
  prefix = genomeBase == genomeRef ? sample : sample + '_spike'
  opts = params.starOpts
  """
  STAR --version &> v_star.txt
  STAR --genomeDir $index \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --readFilesCommand zcat \
       --runDirPerm All_RWX \
       --outSAMunmapped Within \
       --outTmpDir /local/scratch/rnaseq_\$(date +%d%s%S%N) \
       --outFileNamePrefix $prefix  \
       --outSAMattrRGline ID:$prefix SM:$prefix LB:Illumina PL:Illumina \
       ${opts}
  """
}

if (params.aligner == "bowtie2"){
  chAlignReads = chAlignReadsBowtie2
  chMappingMqc = chBowtie2Mqc
} else if (params.aligner == "bwa-mem"){
  chAlignReads = chAlignReadsBwa
  chMappingMqc = chBwaMqc
} else if (params.aligner == "star"){
  chAlignReads = chAlignReadsStar
  chMappingMqc = chStarMqc
}


if (params.inputBam){
  chFastqcMqc = Channel.empty()
  chMappingMqc = Channel.empty()
}

/****************
 * Spike-in
 */

spikes_poor_alignment = []
def checkMappingLog(logs, t='1') {
  def nb_ref = 0;
  def nb_spike = 0;
  def percent_spike = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /Reads on ref\t([\d\.]+)/)) {
      nb_ref = matcher[0][1]                                                                                                                                                                   
    }else if ((matcher = line =~ /Reads on spike\t([\d\.]+)/)) {
      nb_spike = matcher[0][1]                                                                                                                                                                 
    }
  }
  logname = logs.getBaseName() - '_ref_bamcomp.mqc'
  percent_spike = nb_spike.toFloat() / (nb_spike.toFloat() + nb_ref.toFloat()) * 100
  percent_spike = percent_spike.round(5)
  if(percent_spike.toFloat() <= t.toFloat() ){
      log.info "#################### VERY POOR SPIKE ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_spike}% <<"
      spikes_poor_alignment << logname
      return false
  }else {
      log.info "          Passed alignment ($logname)   >> ${percent_spike}% <<"
      return true
  }
}

if (useSpike){

   /* Split and rebuild Channel to be sure of order between bams */
   chAlignRef = Channel.create()
   chAlignSpike = Channel.create()
   chAlignReads.choice( chAlignSpike, chAlignRef ){ it -> it[1] =~ 'spike' ? 1 : 0 }

   chCompAln = chAlignRef
      .join(chAlignSpike)

   // Merging, if necessary reference aligned reads and spike aligned reads
   process compareRefSpike{
     tag "${sample}"
     label 'compbam'
     label 'minCpu'
     label 'medMem'
     publishDir "${params.outDir}/spike", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename }

     input:
     set val(sample), file(unsortedBamSpike), file(unsortedBamRef) from chCompAln

     output:
     set val(sample), file('*_clean.bam') into chRefBams
     set val(sampleSpike), file('*.mqc'), file('*_clean_spike.bam') into chSpikeBams
     file ('*.mqc') into chMappingSpikeMqc
     file ('*.log') into chMappingSpikeLog

     script:
     sampleSpike = sample + '_spike'
     """
     samtools sort $unsortedBamRef -n -@ ${task.cpus} -T ${sample}_ref -o ${sample}_ref_sorted.bam
     samtools sort $unsortedBamSpike -n -@ ${task.cpus} -T ${sample}_spike -o ${sample}_spike_sorted.bam
     compareAlignments.py -a ${sample}_ref_sorted.bam -b ${sample}_spike_sorted.bam -oa ${sample}_clean.bam -ob ${sample}_clean_spike.bam 2> ${sample}_compareAlignments.log
     """
   }

  // Filter removes all 'aligned' channels that fail the check
  chSpikeBams
        .filter { sample, logs, bams -> checkMappingLog(logs, t="$params.spikePercentFilter") }
        .map { row -> [row[0], row[2]]}
        .set { chSpikeCheckBams }


  // concat spike and ref bams
  chRefBams
    .concat(chSpikeCheckBams)
    .set {chAllBams}
}else{
  chAllBams = chAlignReads
  chMappingSpikeMqc = Channel.empty()
}


/*
 * Sorting BAM files
 */

process bamSort{
  tag "${prefix}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if ( filename.endsWith("stats") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null
            }

  input:
  set val(prefix), file(unsortedBam) from chAllBams

  output:
  set val(prefix), file('*sorted.{bam,bam.bai}') into chSortBams
  file("*stats") into chStats
  file("*mqc") into chStatsMqc
  file("v_samtools.txt") into chSamtoolsVersionBamSort

  script:
  """
  samtools --version &> v_samtools.txt
  samtools sort $unsortedBam -@ ${task.cpus} -T ${prefix} -o ${prefix}_sorted.bam
  samtools index ${prefix}_sorted.bam
  samtools flagstat ${prefix}_sorted.bam > ${prefix}_sorted.flagstats
  samtools idxstats ${prefix}_sorted.bam > ${prefix}_sorted.idxstats
  samtools stats ${prefix}_sorted.bam > ${prefix}_sorted.stats

  aligned="\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${prefix}_sorted.bam)"
  hqbam="\$(samtools view -@ $task.cpus -F 0x100 -F 0x800 -F 0x4 -q 10 -c ${prefix}_sorted.bam)"
  lqbam="\$((\$aligned - \$hqbam))"
  echo -e "Mapped,\${aligned}" > ${prefix}_mappingstats.mqc
  echo -e "HighQual,\${hqbam}" >> ${prefix}_mappingstats.mqc
  echo -e "LowQual,\${lqbam}" >> ${prefix}_mappingstats.mqc
  """
}

/*
 * Marking duplicates
 */

process markDuplicates{
  tag "${prefix}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && !filename.endsWith(".bam.bai") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null
            }

  input:
  set val(prefix), file(sortedBams) from chSortBams

  output:
  set val(prefix), file("*marked.{bam,bam.bai}") into chMarkedBams, chMarkedBamsFilt, chMarkedPreseq
  set val(prefix), file("*marked.flagstat") into chMarkedFlagstat
  file "*marked.{idxstats,stats}" into chMarkedStats
  file "*metrics.txt" into chMarkedPicstats
  file("v_picard.txt") into chPicardVersion

  script:
  """
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  picard -Xmx4g MarkDuplicates \\
    INPUT=${sortedBams[0]} \\
    OUTPUT=${prefix}_marked.bam \\
    ASSUME_SORTED=true \\
    REMOVE_DUPLICATES=false \\
    METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR=tmpdir
  samtools index ${prefix}_marked.bam
  samtools idxstats ${prefix}_marked.bam > ${prefix}_marked.idxstats
  samtools flagstat ${prefix}_marked.bam > ${prefix}_marked.flagstat
  samtools stats ${prefix}_marked.bam > ${prefix}_marked.stats
  """
}

/*
 * Preseq (before alignment filtering and only on ref mapped reads)
 */

process preseq {
  tag "${prefix}"
  label 'preseq'
  label 'lowCpu'
  label 'minMem'
  publishDir "${params.outDir}/preseq", mode: 'copy'

  when:
  !params.skipPreseq

  input:
  set val(prefix), file(bam) from chMarkedPreseq.filter{ it[0][-5..-1] != "spike" }

  output:
  file "*.ccurve.txt" into chPreseqStats
  file("v_preseq.txt") into chPreseqVersion

  script:
  defectMode = params.preseqDefect ? '-D' : ''
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -v $defectMode -output ${prefix}.ccurve.txt -bam ${bam[0]} -e 200e+06
  """
}

/*
 * BAM Filtering
 */

process bamFiltering {
  tag "${prefix}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "stats/$filename"
             else if (filename.endsWith("_filtered.bam") || (filename.endsWith("_filtered.bam.bai"))) filename
             else null}

  input:
  set val(prefix), file(markedBam) from chMarkedBamsFilt

  output:
  set val(prefix), file("*filtered.{bam,bam.bai}") into chFilteredBams
  set val(prefix), file("*filtered.flagstat") into chFilteredFlagstat
  file "*filtered.{idxstats,stats}" into chFilteredStats
  file("v_samtools.txt") into chSamtoolsVersionBamFiltering

  script:
  filterParams = params.singleEnd ? "-F 0x004" : params.keepSingleton ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
  dupParams = params.keepDups ? "" : "-F 0x0400"
  mapqParams = params.mapq > 0 ? "-q ${params.mapq}" : ""
  nameSortBam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
  """
  samtools --version &> v_samtools.txt
  samtools view \\
    $filterParams \\
    $dupParams \\
    $mapqParams \\
    -b ${markedBam[0]} > ${prefix}_filtered.bam
  samtools index ${prefix}_filtered.bam
  samtools flagstat ${prefix}_filtered.bam > ${prefix}_filtered.flagstat
  samtools idxstats ${prefix}_filtered.bam > ${prefix}_filtered.idxstats
  samtools stats ${prefix}_filtered.bam > ${prefix}_filtered.stats
  $nameSortBam
  """
}

chFilteredBams
  .into{ chBams; chGroupBamNameFeatCounts }
chFilteredFlagstat
  .set{ chFlagstat }

/*
 * Separate sample BAMs and spike BAMs
 */

chFlagstatChip=Channel.create()
chFlagstatSpikes=Channel.create()
chFlagstat.choice( chFlagstatSpikes, chFlagstatChip ) { it -> it[0] =~ 'spike' ? 0 : 1 }

chBamsChip=Channel.create()
chBamsSpikes=Channel.create()
chBams.choice( chBamsSpikes, chBamsChip ) { it -> it[0] =~ 'spike' ? 0 : 1 }

// Preparing all ChIP data for further analysis
chBamsChip
  .dump (tag:'cbams')
  .into { chBamsMacs1; chBamsMacs2; chBamsPPQT;
          chBamsBigWig; chBamsBigWigSF; 
          chBamDTCor ; chBaiDTCor ; chSampleDTCor ;
          chBamDTFingerprint ; chBaiDTFingerprint ; chSampleDTFingerprint ;
          chBamsCounts }

chFlagstatChip
  .into { chFlagstatMacs; chFlagstatMqc }

/*
 * PhantomPeakQualTools QC
 */

process PPQT{
  tag "${prefix}"
  label 'ppqt'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/ppqt", mode: "copy"

  when:
  !params.skipPPQT

  input:
  set val(prefix), file(bams) from chBamsPPQT
  file sppCorrelationHeader from chPpqtCorHeader
  file sppNSCHeader from chPpqtNSCHeader
  file sppRSCHeader from chPpqtRSCHeader

  output:
  file '*.pdf' into chPpqtPlot
  file '*.spp.out' into chPpqtOutMqc
  file '*_mqc.tsv' into chPpqtCsvMqc
  file("v_R.txt") into chPPQTVersion

  script:
  """
  R --version &> v_R.txt
  RUN_SPP=`which run_spp.R`
  Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bams[0]}" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
  cp $sppCorrelationHeader ${prefix}_spp_correlation_mqc.tsv
  Rscript -e "load('${prefix}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
  awk -v OFS='\t' '{print "${prefix}", \$9}' ${prefix}.spp.out | cat $sppNSCHeader - > ${prefix}_spp_nsc_mqc.tsv
  awk -v OFS='\t' '{print "${prefix}", \$10}' ${prefix}.spp.out | cat $sppRSCHeader - > ${prefix}_spp_rsc_mqc.tsv
  """
}

/* 
 * BigWig Tracks
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/bigWig", mode: "copy",
    saveAs: {filename ->
    	     if ( filename.endsWith(".bigwig") ) "$filename"
             else null}

  input:
  set val(prefix), file(filteredBams) from chBamsBigWig
  file(BLbed) from chBlacklistBigWig.collect()

  output:
  set val(prefix), file('*.bigwig') into chBigWig
  file("v_deeptools.txt") into chDeeptoolsVersion

  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  }
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  bamCoverage --version &> v_deeptools.txt
  nbreads=\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${filteredBams[0]})
  sf=\$(echo "10000000 \$nbreads" | awk '{printf "%.2f", \$1/\$2}')

  bamCoverage -b ${filteredBams[0]} \\
              -o ${prefix}_norm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
              ${extend} \\
              --scaleFactor \$sf
  """
}

/*
 * Spike-in normalization
 */

if (useSpike){
  chBamsSpikes
    .into{chBamsSpikesBam; chBamsSpikesBai}

  process getSpikeCountPerBin {
    label 'deeptools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy"

    input:
    file(allBams) from chBamsSpikesBam.map{it[1][0]}.collect()
    file(allBai) from chBamsSpikesBai.map{it[1][1]}.collect()

    output:
    file "readCounts_spike_10kbins.tab" into chTabCounts

    script:
    """
    multiBamSummary bins \\
                   -b $allBams \\
                   --binSize 10000 \\
                   -o results.npz \\
                   --outRawCounts readCounts_spike_10kbins.tab
    """
  }

 process getSpikeScalingFactor {
    label 'r'
    label 'minCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy"

    input:
    file(tab) from chTabCounts

    output:
    file "*.sf" into chTabSF

    script:
    """
    getDESeqSF.r ${tab}
    """
  }

  chTabSF
    .splitCsv(header:false, sep:',')
    .map { row -> [row[0], row[1]]}
    .set{chScaleFactor}

  chBamsBigWigSF
    .combine(chScaleFactor)
    .filter{it[0] == it[2]}
    .map { it -> it[0,1,3]}
    .set{chBigWigScaleFactor}

  process bigWigSpikeNorm{
    tag "${prefix}"
    label 'deeptools'
    label 'medCpu'
    label 'medMem'
    publishDir "${params.outDir}/bigWigSpike", mode: "copy",
      saveAs: {filename ->
        if ( filename.endsWith(".bigwig") ) "$filename"
        else null}
 
    input:
    set val(prefix), file(filteredBams), val(normFactor) from chBigWigScaleFactor
    file(BLbed) from chBlacklistBigWigSpike.collect()

    output:
    set val(prefix), file('*.bigwig') into chBigWigSF

    script:
    if (params.singleEnd){
      extend = params.fragmentSize > 0 && !params.noReadExtension ? "--extendReads ${params.fragmentSize}" : ""
    }else{
      extend = params.noReadExtension ? "" : "--extendReads"
    }
    blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : ""
    effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
    """
    bamCoverage -b ${filteredBams[0]} \\
                -o ${prefix}_spikenorm.bigwig \\
                -p ${task.cpus} \\
                 ${blacklistParams} \\
                 ${effGsize} \\
                 ${extend} \\
                --scaleFactor ${normFactor}
    """
  }
}

/*
 * DeepTools QC
 */

process deepToolsComputeMatrix{
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/computeMatrix", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  set val(prefix), file(bigwig) from chBigWig
  file geneBed from chGeneBedDeeptools.collect()

  output:
  file("*.{mat,gz,pdf}") into chDeeptoolsSingle
  file("*mqc.tab") into chDeeptoolsSingleMqc

  script:
  """
  computeMatrix scale-regions \\
                -R ${geneBed} \\
                -S ${bigwig} \\
                -o ${prefix}_matrix.mat.gz \\
                --outFileNameMatrix ${prefix}.computeMatrix.vals.mat \\
                --downstream 2000 --upstream 2000 --skipZeros --binSize 100\\
                -p ${task.cpus}

  plotProfile -m ${prefix}_matrix.mat.gz \\
              -o ${prefix}_bams_profile.pdf \\
              --outFileNameData ${prefix}.plotProfile.tab
  
  sed -e 's/.0\t/\t/g' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}_plotProfile_mqc.tab
  """
}

process deepToolsCorrelationQC{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/correlationQC", mode: "copy"

  when:
  allPrefix.size() >= 2 && !params.skipDeepTools

  input:
  file(allBams) from chBamDTCor.map{it[1][0]}.collect()
  file(allBai) from chBaiDTCor.map{it[1][1]}.collect()
  val (allPrefix) from chSampleDTCor.map{it[0]}.collect()
  file(BLbed) from chBlacklistCorrelation.ifEmpty([])

  output:
  file "bams_correlation.pdf" into chDeeptoolsCorrel
  file "bams_correlation.tab" into chDeeptoolsCorrelMqc
  
  script:
  blacklistParams = params.blacklist ? "--blackListFileName ${BLbed}" : "" 
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ")
  allPrefix = allPrefix.replace("]","")
  """
  multiBamSummary bins -b $allBams \\
                       --binSize=50000 \\
                        -o bams_summary.npz \\
                        -p ${task.cpus} \\
                        ${blacklistParams}

  plotCorrelation -in bams_summary.npz \\
                  -o bams_correlation.pdf \\
                  -c spearman -p heatmap -l $allPrefix \\
                  --outFileCorMatrix bams_correlation.tab
  """
}

process deepToolsFingerprint{
  label 'deeptools'
  label 'highCpu'
  label 'lowMem'
  publishDir "${params.outDir}/deepTools/fingerprintQC", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  file(allBams) from chBamDTFingerprint.map{it[1][0]}.collect()
  file(allBai) from chBaiDTFingerprint.map{it[1][1]}.collect()
  val (allPrefix) from chSampleDTFingerprint.map{it[0]}.collect()
 
  output:
  file "bams_fingerprint.pdf" into chDeeptoolsFingerprint
  file "plotFingerprint*" into chDeeptoolsFingerprintMqc 
 
  script:
  if (params.singleEnd){
    extend = params.fragmentSize > 0 ? "--extendReads ${params.fragmentSize}" : ""
  }else{
    extend = params.noReadExtension ? "" : "--extendReads"
  }
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ") 
  allPrefix = allPrefix.replace("]","")
  """
  plotFingerprint -b $allBams \\
                  -plot bams_fingerprint.pdf \\
                  -p ${task.cpus} \\
                  -l $allPrefix \\
                  ${extend} \\
                  --skipZeros \\
                  --outRawCounts plotFingerprint.raw.txt \\
                  --outQualityMetrics plotFingerprint.qmetrics.txt
  """
}

/*#########################################################################
  /!\ From this point, 'design' is mandatory /!\
###########################################################################*/

/*
 * Prepare channels
 */

if (params.design){
  chBamsMacs1
    .join(chFlagstatMacs)
    .combine(chNoInput.concat(chBamsMacs2))
    .set { chBamsMacs1 }

  chDesignControl
    .combine(chBamsMacs1)
    .filter { it[0] == it[5] && it[1] == it[8] }
    .map { it ->  it[2..-1] }
    .into { chGroupBamMacsSharp; chGroupBamMacsBroad; chGroupBamMacsVeryBroad}

  chGroupBamMacsSharp
    .filter { it[2] == 'sharp' }
    .dump(tag:'peakCall')
    .set { chGroupBamMacsSharp }

  chGroupBamMacsBroad
    .filter { it[2] == 'broad' }
    .dump(tag:'peakCall')
    .set { chGroupBamMacsBroad }

  chGroupBamMacsVeryBroad
    .filter { it[2] == 'very-broad' }
    .dump(tag:'peakCall')
    .set { chGroupBamMacsVeryBroad }
}else{
  chGroupBamMacsSharp=Channel.empty()
  chGroupBamMacsBroad=Channel.empty()
  chGroupBamMacsVeryBroad=Channel.empty()
}

/***********************
 * Peak calling 
 */

/*
 * MACS2 - sharp mode
 */

process sharpMACS2{
  tag "${sampleID} - ${controlID}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/sharp", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else filename
            }
 
  when:
  !params.skipPeakCalling && params.effGenomeSize

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsSharp
  file peakCountHeader from chPeakCountHeaderSharp.collect()
  file fripScoreHeader from chFripScoreHeaderSharp.collect()

  output:
  file("*.xls") into chMacsOutputSharp
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.narrowPeak") into chPeaksMacsSharp
  file "*_mqc.tsv" into chMacsCountsSharp
  file("v_macs2.txt") into chMacs2VersionMacs2Sharp

  script:
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  """
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    -f $format \\
    -g $params.effGenomeSize \\
    -n $sampleID \\
    --SPMR --trackline --bdg \\
    --keep-dup all
  cat ${sampleID}_peaks.narrowPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.narrowPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
 }

/*
 * MACS2  - Broad
 */

process broadMACS2{
  tag "${sampleID} - ${controlID}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/broad", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else filename
            }
   
  when:
  !params.skipPeakCalling && params.effGenomeSize

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsBroad
  file peakCountHeader from chPeakCountHeaderBroad.collect()
  file fripScoreHeader from chFripScoreHeaderBroad.collect()

  output:
  file("*.xls") into chMacsOutputBroad
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chPeaksMacsBroad
  file "*_mqc.tsv" into chMacsCountsBroad
  file("v_macs2.txt") into chMacs2VersionMacs2Broad

  script:
  broad = "--broad --broad-cutoff ${params.broadCutoff}"
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  """
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    ${broad} \\
    -f $format \\
    -g $params.effGenomeSize \\
    -n $sampleID \\
    --SPMR --trackline --bdg \\
    --keep-dup all
  cat ${sampleID}_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.broadPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
}

/*
 * EPIC2 - very broad
 */

process veryBroadEpic2{
  tag "${sampleID} - ${controlID}"
  label 'epic2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/very-broad", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else filename
            }
 
  when:
  !params.skipPeakCalling && params.effGenomeSize

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsVeryBroad
  file peakCountHeader from chPeakCountHeaderVeryBroad.collect()
  file fripScoreHeader from chFripScoreHeaderVeryBroad.collect()
  file chromsize from chChromSize.collect()

  output:
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chPeaksEpic
  file ("*.out") into chEpicOutput
  file "*_mqc.tsv" into chMacsCountsVbroad
  file("v_epic2.txt") into chEpic2Version

  script:
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  opts = (params.singleEnd && params.fragmentSize > 0) ? "--fragment-size ${params.fragmentSize}" : "" 
  """
  epic2 --version &> v_epic2.txt
  epic2 -t ${sampleBam[0]} \\
    ${ctrl} \\
    --chromsizes ${chromsize} \\
    --effective-genome-fraction ${params.effGenomeSize} \\
    -a \\
    --bin-size 200 \\
    --gaps-allowed 3 \\
    --false-discovery-rate-cutoff 0.05 \\
    -o ${sampleID}_epic.out

  echo "track type=broadPeak name=\"${sampleID}\" description=\"${sampleID}\" nextItemButton=on" > ${sampleID}_peaks.broadPeak
  awk -v id="${sampleID}" 'NF<10{fc=0;pval=-1;qval=-1}NF==10{fc=\$10;pval=\$4;qval=\$9}NR>1{OFS="\t"; print \$1,\$2,\$3,id"_peak_"NR,\$5,".",fc,pval,qval}' ${sampleID}_epic.out >> ${sampleID}_peaks.broadPeak
  cat ${sampleID}_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.broadPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
}

// Join the results of all peaks callers
chPeaksMacsSharp
  .mix(chPeaksMacsBroad, chPeaksEpic)
  .into{ chPeaksHomer; chIDRpeaks; chIDRid; chPeakQC }

/************************************
 * Peaks Annotation
 */

process peakAnnoHomer{
  tag "${sampleID}"
  label 'homer'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/annotation/", mode: 'copy'

  when:
  !params.skipPeakAnno

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file (peakfile) from chPeaksHomer
  file gtfFile from chGtfHomer.collect()
  file fastaFile from chFastaHomer.collect()

  output:
  file "*.txt" into chHomerMqc

  script:
  """
  annotatePeaks.pl $peakfile \\
        $fastaFile \\
        -gtf $gtfFile \\
        -cpu ${task.cpus} \\
        > ${sampleID}_annotated_peaks.txt
  """
}


/*
 * Peak calling & annotation QC
 */

process peakQC{
  label 'r'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/peakCalling/QC/", mode: 'copy'

  when:
  !params.skipPeakQC && params.design

  input:
  file peaks from chPeakQC.collect{ it[-1] }
  file annotations from chHomerMqc.collect()
  file peakHeader from chPeakAnnotationHeader

  output:
  file "*.{txt,pdf}" into chMacsQcOutput
  file "*.tsv" into chPeakMqc

  script:
  """
  plot_macs_qc.r \\
    -i ${peaks.join(',')} \\
    -s ${peaks.join(',').replaceAll("_peaks.narrowPeak","").replaceAll("_peaks.broadPeak","")} \\
    -o ./ \\
    -p peak
  plot_homer_annotatepeaks.r \\
    -i ${annotations.join(',')} \\
    -s ${annotations.join(',').replaceAll("_annotated_peaks.txt","")} \\
    -o ./ \\
    -p annotatePeaks
  cat $peakHeader annotatePeaks.summary.txt > annotatedPeaks.summary_mqc.tsv
  """
}


/**********************************
 * Per Group Analysis
 */

/*
 * Irreproducible Discovery Rate
 */

chIDRpeaks
  .map { it -> [it[0],it[4]] }
  .groupTuple()
  .dump (tag:'rep')
  .set{ chPeaksPerGroup }

process IDR{
  tag "${group}"
  label 'idr'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/IDR", mode: 'copy'

  when:
  allPeaks.toList().size > 1 && !params.skipIDR

  input:
  set val(group), file(allPeaks) from chPeaksPerGroup

  output:
  file "*idrValues.txt" into chIdr
  file "*log.txt" into chMqcIdr
  file("v_idr.txt") into chIdrVersion

  script:
  peaktype = allPeaks[0].toString()
  peaktype = peaktype.substring(peaktype.lastIndexOf(".") + 1)
  """
  idr --version &> v_idr.txt
  idr --samples ${allPeaks} \\
      --input-file-type  ${peaktype} \\
      -o ${group}_idrValues.txt \\
      -l ${group}_log.txt \\
      --plot
  """
}


/**************************************
 * Feature counts
 */

process prepareAnnotation{
  label 'unix'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.outDir}/featCounts/", mode: "copy"

  when:
  !params.skipFeatCounts

  input:
  file(bed) from chGenePrepareAnnot.collect()

  output:
  file("*tss.bed") into chTSSFeatCounts

  script:
  prefix = bed.toString() - ~/(.bed)?$/
  """
  awk -F"\t" -v win=${params.tssSize} 'BEGIN{OFS="\t"} \$6=="+"{s=\$2-win;e=\$2+win;if(s<0){s=0}; print \$1,s,e,\$4,\$5,\$6} \$6=="-"{print \$1,\$3-win,\$3+win,\$4,\$5,\$6}' ${bed} > ${prefix}_tss.bed
  """
}
    
process featureCounts{
  tag "${bed}"
  label 'featureCounts'
  label 'medCpu'
  label 'lowMem'
  publishDir "${params.outDir}/featCounts/", mode: "copy"

  when:
  !params.skipFeatCounts

  input:
  file(bams) from chBamsCounts.map{items->items[1][0]}.collect()
  each file(annot) from chGeneFeatCounts.concat(chTSSFeatCounts)

  output:
  file("*csv") into chFeatCounts
  file("*summary") into chFeatCountsMqc
  file("v_featurecounts.txt") into chFeaturecountsVersion

  script:
  prefix = annot.toString() - ~/(\.bed)?$/
  paramsPairedEnd = params.singleEnd ? '' : '-p -C -B'
  """
  featureCounts -v &> v_featurecounts.txt
  awk '{OFS="\t";print \$4,\$1,\$2,\$3,\$6}' ${annot} > ${prefix}.saf
  featureCounts -a ${prefix}.saf -F SAF \\
                -o allchip_counts_${prefix}.csv \\
                -T ${task.cpus} \\
                -s 0 ${paramsPairedEnd} \\
                -O ${bams} 2> featureCounts_${prefix}.log
  """
}


/*
 * MultiQC
 */
process getSoftwareVersions{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file 'v_fastqc.txt' from chFastqcVersion.first().ifEmpty([])
  file 'v_bwa.txt' from chBwaVersion.first().ifEmpty([])
  file 'v_bowtie2.txt' from chBowtie2Version.first().ifEmpty([])
  file 'v_star.txt' from chStarVersion.first().ifEmpty([])
  file 'v_samtools.txt' from chSamtoolsVersionBamSort.concat(chSamtoolsVersionBamFiltering).first().ifEmpty([])
  file 'v_picard.txt' from chPicardVersion.first().ifEmpty([])
  file 'v_macs2.txt' from chMacs2VersionMacs2Sharp.concat(chMacs2VersionMacs2Broad).first().ifEmpty([])
  file 'v_epic2.txt' from chEpic2Version.first().ifEmpty([])
  file 'v_preseq.txt' from chPreseqVersion.first().ifEmpty([])
  file 'v_idr.txt' from chIdrVersion.first().ifEmpty([])
  file 'v_R.txt' from chPPQTVersion.first().ifEmpty([])
  file 'v_deeptools.txt' from chDeeptoolsVersion.first().ifEmpty([])
  file 'v_featurecounts.txt' from chFeaturecountsVersion.first().ifEmpty([])

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}


process workflowSummaryMqc {
  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/data-analysis/chip-seq'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}


process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file multiqcConfig from chMultiqcConfig
  file design from chDesignMqc.collect().ifEmpty([])
  file metadata from chMetadata.ifEmpty([])
  file ('software_versions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflow_summary/*') from workflowSummaryYaml.collect()
  file ('fastqc/*') from chFastqcMqc.collect().ifEmpty([])
  file ('mapping/*') from chMappingMqc.collect().ifEmpty([])
  file ('mapping/*') from chMappingSpikeMqc.collect().ifEmpty([])
  file ('mapping/*') from chMarkedPicstats.collect().ifEmpty([])
  file ('mapping/*') from chStatsMqc.collect().ifEmpty([])
  file ('preseq/*') from chPreseqStats.collect().ifEmpty([])
  file ('ppqt/*') from chPpqtOutMqc.collect().ifEmpty([])
  file ('ppqt/*') from chPpqtCsvMqc.collect().ifEmpty([])
  file ('deepTools/*') from chDeeptoolsSingleMqc.collect().ifEmpty([])
  file ("deepTools/*") from chDeeptoolsCorrelMqc.collect().ifEmpty([])
  file ("deepTools/*") from chDeeptoolsFingerprintMqc.collect().ifEmpty([])
  file ('peakCalling/sharp/*') from chMacsOutputSharp.collect().ifEmpty([])
  file ('peakCalling/broad/*') from chMacsOutputBroad.collect().ifEmpty([])
  file ('peakCalling/sharp/*') from chMacsCountsSharp.collect().ifEmpty([])
  file ('peakCalling/broad/*') from chMacsCountsBroad.collect().ifEmpty([])
  file ('peakCalling/very-broad/*') from chMacsCountsVbroad.collect().ifEmpty([])
  file ('peakQC/*') from chPeakMqc.collect().ifEmpty([])
  
  output:
  file splan
  file "*_report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_chipseq_report" : "--filename chipseq_report" 
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  isPE = params.singleEnd ? "" : "-p"
  designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content -m fastqc -m bowtie2 -m star -m preseq -m picard -m phantompeakqualtools -m deeptools -m macs2 -m homer"
  """
  stats2multiqc.sh -s ${splan} ${designOpts} -a ${params.aligner} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')" 
  mqc_header.py --splan ${splan} --name "ChIP-seq" --version ${workflow.manifest.version} ${metadataOpts} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

/*
 * Sub-routine
 */
process outputDocumentation {
    label 'python'
    label 'minCpu'
    label 'minMem'
    publishDir "${params.summaryDir}/", mode: 'copy'

    input:
    file output_docs from chOutputDocs
    file images from chOutputDocsImages

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}


/* Creates a file at the end of workflow execution */
workflow.onComplete {

  /*pipeline_report.html*/

  def report_fields = [:]
  report_fields['pipeline'] = workflow.manifest.name
  report_fields['version'] = workflow.manifest.version
  report_fields['runName'] = customRunName ?: workflow.runName
  report_fields['success'] = workflow.success
  report_fields['dateComplete'] = workflow.complete
  report_fields['duration'] = workflow.duration
  report_fields['exitStatus'] = workflow.exitStatus
  report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  report_fields['errorReport'] = (workflow.errorReport ?: 'None')
  report_fields['commandLine'] = workflow.commandLine
  report_fields['projectDir'] = workflow.projectDir
  report_fields['summary'] = summary
  report_fields['summary']['Date Started'] = workflow.start
  report_fields['summary']['Date Completed'] = workflow.complete
  report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
  report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision

  report_fields['spikes_poor_alignment'] = spikes_poor_alignment

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
  def txt_template = engine.createTemplate(tf).make(report_fields)
  def report_txt = txt_template.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/onCompleteTemplate.html")
  def html_template = engine.createTemplate(hf).make(report_fields)
  def report_html = html_template.toString()

  // Write summary e-mail HTML to a file
  def output_d = new File( "${params.summaryDir}/" )
  if( !output_d.exists() ) {
    output_d.mkdirs()
  }
  def output_hf = new File( output_d, "pipelineReport.html" )
  output_hf.withWriter { w -> w << report_html }
  def output_tf = new File( output_d, "pipelineReport.txt" )
  output_tf.withWriter { w -> w << report_txt }

  /*oncomplete file*/
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'

  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  if(spikes_poor_alignment.size() > 0){
    log.info "[chIP-seq] WARNING - ${spikes_poor_alignment.size()} samples skipped due to poor alignment scores!"
  }
                                                                                                                                                                                         
  if(workflow.success){
    log.info "[ChIP-seq] Pipeline Complete"
  }else{
    log.info "[ChIP-seq] FAILED: $workflow.runName"
  }
}
