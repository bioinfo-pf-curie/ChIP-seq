#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
            Chip-seq
========================================================================================
Chip-seq Analysis Pipeline.
#### Homepage / Documentation
https://gitlab.curie.fr/chipseq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/dev_message.txt")
  log.info devMess.text
  }

  log.info"""

  Chip-seq v${workflow.manifest.version}
  ======================================================================

  Usage:

  nextflow run main.nf -profile test,toolsPath --genome 'hg19' --singleEnd

  Mandatory arguments:
  --samplePlan             Path to sample plan file if '--reads' is not specified
  --genome                 Name of iGenomes reference
  -profile                 Configuration profile to use. Can use multiple (comma separated)
                           Available: conda, docker, singularity, awsbatch, test, toolsPath and more.

  Inputs:
  --design                 Path to design file for downstream analysis
  --singleEnd              Specifies that the input is single end reads
  --spike                  Indicates if the experiment includes a spike-in normalization.
                           Default : false. Available : 'spike' to use metagenome with reference genome
                           '[spike genome]' to use a specific second genome

  Tools:
  --aligner                Alignment tool to use ['bwa-mem', 'star', 'bowtie2']. Default: 'bwa-mem'
  
  Filtering:
  --mapQ                   Minimum mapping quality to consider
  --keepDups               Do not remove duplicates afer marking
  --blacklist              Path to black list regions (.bed)

  Peak calling
  --macsGzise              Reference genome size for MACS2

  References           If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --fasta                  Path to Fasta reference

  Indexes              Path to the indexes for aligners
  --starIndex              Index for STAR aligner
  --bwaIndex               Index for BWA MEM aligner
  --bowtie2Index           Index for Bowtie2 aligner

  Annotation
  --bed                    BED annotation file. Used with samtools to filter the reads, and in Deeptools ComputeMatrix function
  --gtf                    GTF annotation file. Used in HOMER peak annotation

  Skip options:        All are false by default
  --skipMultiqc            Skips final report writing
  --skipFastqc             Skips fastQC
  --skipAlignment          Skips alignments
  --skipPreseq             Skips preseq QC
  --skipFiltering          Skips filtering
  --skipPPQT               Skips phantompeakqualtools QC
  --skipDeepTools          Skips deeptools QC
  --skipPeakcalling        Skips peak calling
  --skipPeakanno           Skips peak annotation
  --skipIDR                Skips IDR QC
  --skipFeatCounts         Skips feature count
  --skipMultiQC            Skips MultiQC step

  Other options:
  --outdir                 The output directory where the results will be saved
  --email                  Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
  -name                    Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

  """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help messsage
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes
if (params.spike == 'spike'){
  genomeRef = params.genome + '_spike'
}
else {
  genomeRef = params.genome
}

params.fasta = genomeRef ? params.genomes[ genomeRef ].fasta ?: false : false
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .into{chFastaHomer; chFastaBwa; chFastaBt2; chFastaStar}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}


if (params.spike && params.spike != 'spike'){
  params.spikeFasta = params.spike ? params.genomes[ params.spike ].fasta ?: false : false
  if ( params.spikeFasta ){
    Channel
      .fromPath(params.spikeFasta, checkIfExists: true)      
      .into{chSpikeFaBwa; chSpikeFaBt2; chSpikeFaStar}
  } else {
    exit 1, "Spike-in fasta file not found: ${params.spikeFasta}"
  }
  println "Fasta ok"
}

//Bwa-mem Index
params.bwaIndex = genomeRef ? params.genomes[ genomeRef ].bwaIndex ?: false : false
if (params.bwaIndex){
  lastPath = params.fasta.lastIndexOf(File.separator)
  bwaDir =  params.bwaIndex.substring(0,lastPath+1)
  bwaBase = params.bwaIndex.substring(lastPath+1)
  Channel
    .fromPath(bwaDir, checkIfExists: true)
    .ifEmpty {exit 1, "BWA index file not found: ${params.bwaIndex}"}
    .set { chBwaIndex }
} else {
  exit 1, "BWA index file not found: ${params.bwaIndex}"
}

if (params.spike && params.spike != 'spike'){
  params.spikeBwaIndex = params.spike ? params.genomes[ params.spike ].bwaIndex ?: false : false
  if (params.spikeBwaIndex){
    lastPath = params.spikeFasta.lastIndexOf(File.separator)
    bwaDirSpike =  params.spikeBwaIndex.substring(0,lastPath+1)
    spikeBwaBase = params.spikeBwaIndex.substring(lastPath+1)
    Channel
      .fromPath(bwaDirSpike, checkIfExists: true)
      .ifEmpty {exit 1, "Spike BWA index file not found: ${params.spikeBwaIndex}"}
      .set { chSpikeBwaIndex }
  } else {
    chSpikeBwaIndex=Channel.empty()
    //exit 1, "Spike BWA index file not found: ${params.spikeBwaIndex}"
  }
  println "Bwa OK"
}

//Bowtie2 indexes
params.bt2Index = genomeRef ? params.genomes[ genomeRef ].bowtie2Index ?: false : false
if (params.bt2Index){
  lastPath = params.fasta.lastIndexOf(File.separator)
  bt2Dir =  params.bt2Index.substring(0,lastPath+1)
  bt2Base = params.bt2Index.substring(lastPath+1)
  Channel
    .fromPath(bt2Dir, checkIfExists: true)
    .ifEmpty {exit 1, "Bowtie2 index file not found: ${params.bt2Index}"}
    .set { chBt2Index }
} else {
  exit 1, "Bowtie2 index file not found: ${params.bt2Index}"
}

if (params.spike && params.spike != 'spike'){
  params.spikeBt2Index = params.spike ? params.genomes[ params.spike ].bowtie2Index ?: false : false
  if (params.spikeBt2Index){
    lastPath = params.spikeFasta.lastIndexOf(File.separator)
    bt2DirSpike =  params.spikeBt2Index.substring(0,lastPath+1)
    spikeBt2Base = params.spikeBt2Index.substring(lastPath+1)
    Channel
      .fromPath(bt2DirSpike, checkIfExists: true)
      .ifEmpty {exit 1, "Spike Bowtie2 index file not found: ${params.spikeBt2Index}"}
      .set { chSpikeBt2Index }
  } else {
    chSpikeBt2Index=Channel.empty()
    //exit 1, "Spike bowtie2 index file not found: ${params.spikeBt2Index}"
  }
  println "Bt2 OK"
}

//Star indexes
params.starIndex = genomeRef ? params.genomes[ genomeRef ].starIndex ?: false : false
if (params.starIndex){
  lastPath = params.fasta.lastIndexOf(File.separator)
  starDir =  params.starIndex.substring(0,lastPath+1)
  starBase = params.starIndex.substring(lastPath+1)
  Channel
    .fromPath(starDir, checkIfExists: true)
    .ifEmpty {exit 1, "STAR index file not found: ${params.starIndex}"}
    .set { chStarIndex }
} else {
  exit 1, "STAR index file not found: ${params.starIndex}"
}

if (params.spike && params.spike != 'spike'){
  params.spikeStarIndex = params.spike ? params.genomes[ params.spike ].starIndex ?: false : false
  if (params.spikeStarIndex){
    lastPath = params.spikeFasta.lastIndexOf(File.separator)
    starDirSpike =  params.spikeStarIndex.substring(0,lastPath+1)
    spikeStarBase = params.spikeStarIndex.substring(lastPath+1)
    Channel
      .fromPath(starDirSpike, checkIfExists: true)
      .ifEmpty {exit 1, "Spike STAR index file not found: ${params.spikeStarIndex}"}
      .set { chSpikeStarIndex }
  } else {
    chSpikeStarIndex = Channel.empty()
    //exit 1, "Spike STAR index file not found: ${params.spikeStarIndex}"
  }
  println "starOk"
}


// Other inputs
params.gtf = genomeRef ? params.genomes[ genomeRef ].gtf ?: false : false
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .into{chGtfHomer; chGtfFeatCounts}
}
else {
  exit 1, "GTF annotation file not specified!"
}

params.geneBed = genomeRef ? params.genomes[ genomeRef ].bed12 ?: false : false
if (params.geneBed) {
  Channel
    .fromPath(params.geneBed, checkIfExists: true)
    .into{chGeneBed; chGeneBedDeeptools}
}

if (params.spike && params.spike != 'spike'){
  params.spikeGeneBed = params.spike ? params.genomes[ params.spike ].bed12 ?: false : false
  if (params.spikeGeneBed) {
    Channel
      .fromPath(params.spikeGeneBed, checkIfExists: true)
      .set{chSpikeGeneBed}
  }
  else{
  exit 1, "BED file not specified!"
  }
  println "bed ok"
}

//PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

//Peak Calling
params.macsGsize = genomeRef ? params.genomes[ genomeRef ].macsGsize ?: false : false
params.blacklist = genomeRef ? params.genomes[ genomeRef ].blacklist ?: false : false

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .into { chPeakCountHeaderSharp; chPeakCountHeaderBroad; chPeakCountHeaderVeryBroad }

Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .into{ chFripScoreHeaderSharp; chFripScoreHeaderBroad; chFripScoreHeaderVeryBroad }

//Channel
//  .fromPath("$baseDir/assets/peak_annotation_header.txt")
//  .set{ chPeakAnnotationHeader }


//Filtering
if (params.singleEnd) {
  chBamtoolsFilterConfig = Channel.fromPath(params.bamtoolsFilterSEConfig, checkIfExists: true)
} else {
  chBamtoolsFilterConfig = Channel.fromPath(params.bamtoolsFilterPEConfig, checkIfExists: true)
}


//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
Channel
  .fromPath(params.outputDoc, checkIfExists: true)
  .set{chOutputDocs}


//Has the run name been specified by the user?
//This has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

//Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/dev_message.txt")
  log.info devMess.text
}

log.info """=======================================================

Chip-seq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'Chip-seq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = customRunName ?: workflow.runName
if (params.samplePlan) {
  summary['SamplePlan']   = params.samplePlan
}
else{
  summary['Reads']        = params.reads
}
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
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
}

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .into { rawReadsFastqc; rawReadsBWA; rawReadsBt2; rawReadsSTAR; rawSpikeReadsBWA; rawSpikeReadsBt2; rawSpikeReadsSTAR }
   }
   params.reads=false
}
else if(params.readPaths){
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


/*
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
       .set{ chSplan; chSplanCheck }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan; chSplanCheck }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { chSplan; chSplanCheck }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { chSplan; chSplanCheck }
   }
}

/*
 * Design file
 */

if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { chDesignCheck; chDesignControl }

  chDesignControl
    .splitCsv(header:false, sep:',')
    .map { row ->
      if(row[1]==""){row[1]='NO_INPUT'}
      return [ row[0], row[1], row[2], row[3], row[4] ]
     }
    .dump(tag:'design')
    .into { chDesignControl; chDesignReplicate }

  // Create special channel to deal with no input cases
  Channel
    .from( ["NO_INPUT", ["NO_FILE","NO_FILE"]] )
    .toList()
    .set{ chNoInput }
}else{
  chDesignControl=Channel.empty()
  chDesignReplicate=Channel.empty()
  chDesignCheck=Channel.empty()
}

/*
 * Check design and sampleplan files
 */

process checkDesign{
  publishDir "${params.outdir}/design", mode: 'copy'

  when:
  params.design

  input:
  file design from chDesignCheck
  file samplePlan from chSplanCheck

  script:
  """
  check_designs.py $design $samplePlan ${params.singleEnd} $baseDir
  """
}


/*
 * FastQC
 */

process fastQC{
  tag "${prefix}"
  publishDir "${params.outdir}/fastQC"

  when:
  !params.skipFastqc

  input:
  set val(prefix), file(reads) from rawReadsFastqc

  output:
  file "*_fastqc.{zip,html}" into chFastqcResults

  script:
  """
  fastqc -q $reads -t 6
  """
}

/*
 * Alignment on reference genome
 */

// BWA-MEM
process bwaMem{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/bwa_alignment"

  when:
  !params.skipAlignment && params.aligner == "bwa-mem"

  input:
  set val(prefix), file(reads) from rawReadsBWA
  file bwaIndex from chBwaIndex.collect()

  output:
  set val(prefix), file("*.bam") into chAlignReadsBwa

  script:
  alnMult = params.spike == 'spike' ? '-a' : ''
  """
  bwa mem -t ${task.cpus} $alnMult $bwaIndex${bwaBase} $reads \\
  | samtools view -b -h -F 256 -o ${prefix}.bam
  """
}

// BOWTIE2
process bowtie2{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/bowtie2_alignment"

  when:
  !params.skipAlignment && params.aligner == "bowtie2"

  input:
  set val(prefix), file(reads) from rawReadsBt2
  file bt2Index from chBt2Index.collect()

  output:
  set val(prefix), file("*.bam") into chAlignReadsBowtie2

  script:
  readCommand = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  alnMult = params.spike ? alnMult = '-k 3' : ''
  """
  bowtie2 -p ${task.cpus} $alnMult -x $bt2Index${bt2Base} $readCommand \\
  | samtools view -b -h -F 256 -F 2048 -o ${prefix}.bam
  """
}

// STAR
process star{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/star_alignment"

  when:
  !params.skipAlignment && params.aligner == "star"

  input:
  set val(prefix), file(reads) from rawReadsSTAR
  file starIndex from chStarIndex.collect()

  output:
  set val(prefix), file('*.bam') into chAlignReadsStar

  script:
  """
  STAR --genomeDir $starIndex/${starBase} --runThreadN ${task.cpus} --readFilesIn $reads --outSAMtype BAM Unsorted --readFilesCommand zcat \\
  | samtools view -b -h -o ${prefix}.bam
  """
}

if (params.aligner == "bowtie2"){
  chAlignReads = chAlignReadsBowtie2
}else if (params.aligner == "bwa-mem"){
  chAlignReads = chAlignReadsBwa
}else if (params.aligner == "star"){
  chAlignReads = chAlignReadsStar
}

/*
 * SPIKES GENOME
 */
if(params.spike && params.spike != 'spike'){
// BWA-MEM
process spikeBwaMem{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/spike/bwa_alignment"

  when:
  params.spike && params.spike != 'spike' && !params.skipAlignment && params.aligner == "bwa-mem"

  input:
  set val(prefix), file(reads) from rawSpikeReadsBWA
  file spikeBwaIndex from chSpikeBwaIndex.collect()

  output:
  set val(spikeprefix), file("*_spike.bam") into chSpikeAlignReadsBwa

  script:
  spikeprefix = "${prefix}_spike"
  """
  bwa mem -t ${task.cpus} $spikeBwaIndex${spikeBwaBase} $reads \\
  | samtools view -b -h -F 256 -F 2048 -o ${spikeprefix}.bam
  """
}

// BOWTIE2
process spikeBowtie2{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/spike/bowtie2_alignment"

  when:
  params.spike && params.spike != 'spike' && !params.skipAlignment && params.aligner == "bowtie2"

  input:
  set val(prefix), file(reads) from rawSpikeReadsBt2
  file spikeBt2Index from chSpikeBt2Index.collect()

  output:
  set val(spikeprefix), file("*.bam") into chSpikeAlignReadsBowtie2

  script:
  spikeprefix = "${prefix}_spike"
  readCommand = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  bowtie2 -p ${task.cpus} -x $spikeBt2Index/${spikeBt2Base} $readCommand \\
  | samtools view -b -h -o ${spikeprefix}.bam
  """
 }
 
// STAR
process spikeStar{
  tag "${prefix}"
  publishDir "${params.outdir}/alignment/spike/star_alignment"

  when:
  params.spike && params.spike != 'spike' && !params.skipAlignment && params.aligner == "star"

  input:
  set val(prefix), file(reads) from rawSpikeReadsSTAR
  file spikeStarIndex from chSpikeStarIndex.collect()

  output:
  set val(spikeprefix), file('*.bam') into chSpikeAlignReadsStar

  script:
  spikeprefix = "${prefix}_spike"
  """
  STAR --genomeDir $spikeStarIndex/${spikeStarBase} --runThreadN ${task.cpus} --readFilesIn $reads --outSAMtype BAM Unsorted --readFilesCommand zcat \\
  | samtools view -b -h -o ${spikeprefix}.bam
  """
}
}

// Split Reference and Spike genomes
if(params.spike && params.spike != 'spike'){

  if (params.aligner == "bowtie2"){
    chSpikeAlignReads = chSpikeAlignReadsBowtie2
  }else if (params.aligner == "bwa-mem"){
    chSpikeAlignReads = chSpikeAlignReadsBwa
  }else if (params.aligner == "star"){
    chSpikeAlignReads = chSpikeAlignReadsStar
  }

  chAlignReads
    .combine(chSpikeAlignReads)
    .filter{ it[0] == it[2][0..-7] }
    .map{ it[0,1,3] }    
    .set{chAllReads}

  // Merging, if necessary reference aligned reads and spike aligned reads
  process compareRefSpike{
    tag "${prefixRef}"
    publishDir "${params.outdir}/alignment/compareRefSpike"

    input:
      set val(prefixRef), file(unsortedBamRef), file(unsortedBamSpike) from chAllReads

    output:
      set val(prefixRef), file('*_ref.bam') into chAlignedReads
      set val(prefixSpike), file('*Spikein.bam') into chSpikeAlignedReads
      file ('*.txt') into chLogCompare

    script:
    prefixSpike = prefixRef + 'Spike'
    """
    samtools sort $unsortedBamRef -n -T ${prefixRef} -o ${prefixRef}_sorted.bam
    samtools sort $unsortedBamSpike -n -T ${prefixSpike} -o ${prefixSpike}_sorted.bam
    compareAlignments.py -IR ${prefixRef}_sorted.bam -IS ${prefixSpike}_sorted.bam -OR ${prefixRef}_ref.bam -OS ${prefixSpike}in.bam -SE ${params.singleEnd}
    """
  }
  
  chAlignedReads
    .mix(chSpikeAlignedReads)
    .dump (tag:'bams')
    .set{chAlignedReads}

} else if(params.spike && params.spike == 'spike'){

  process sepMetagenome{
    tag "${prefix}"
    publishDir "${params.outdir}/alignment/sepMetagenome"

    input:
      set val(prefix), file(unsortedBam) from chAlignReads

    output:
      set val(prefix), file('*_ref.bam') into chAlignedReads
      set val(prefix), file('*Spike.bam') into chSpikeAlignedReads
      file ('*.txt') into chLogSep

    script:
    """
    samtools sort $unsortedBam -n -T ${prefix} -o ${prefix}_sorted.bam
    sepMetagenome.py -I ${prefix}_sorted.bam -OR ${prefix}_ref.bam -OS ${prefix}Spike.bam -SE ${params.singleEnd}
    """
  }

  chAlignedReads
    .mix(chSpikeAlignedReads)
    .dump (tag:'bams')
    .set{chAlignedReads}
} else {
  chAlignReads
    .dump (tag:'bams')
    .set{chAlignedReads}
}

/*
 * Sorting BAM files
 */

process bamSort{
  tag "${name}"
  publishDir path: "${params.outdir}/filtering/sortedBams", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
             else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
             else null
            }
  when:
  !params.skipAlignment

  input:
  set val(prefix), file(unsortedBam) from chAlignedReads

  output:
  set val(name), file('*.{bam,bam.bai}') into chSortBam
  file "*.{flagstat,idxstats,stats}" into chSortBamStats

  script:
  String unsortedBamsName = unsortedBam.toString()
  name = unsortedBamsName.contains('Spike') ? "${prefix}Spike" : "${prefix}"
  """
  samtools sort $unsortedBam -@ ${task.cpus} -T ${name} -o ${name}_sorted.bam
  samtools index ${name}_sorted.bam
  samtools flagstat ${name}_sorted.bam > ${name}_sorted.flagstat
  samtools idxstats ${name}_sorted.bam > ${name}_sorted.idxstats
  samtools stats ${name}_sorted.bam > ${name}_sorted.stats
  """
}

/*
 * Marking duplicates
 */

process markDuplicates{
  tag "${prefix}"
  publishDir path: "${params.outdir}/filtering/markedBams", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
             else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
             else null
            }

  input:
  set val(prefix), file(sortedBams) from chSortBam

  output:
  set val(prefix), file("*.{bam,bam.bai}") into chMarkedBams, chMarkedBamsFilt, chMarkedPreseq
  set val(prefix), file("*.flagstat") into chMarkedFlagstat
  file "*.{idxstats,stats}" into chMarkedStats
  file "*.txt" into chMarkedPicstats

  script:
  """
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
 * Preseq (before alignment filtering)
 */

process preseq {
  tag "${prefix}"
  publishDir "${params.outdir}/preseq"

  when:
  !params.skipPreseq

  input:
  set val(prefix), file(bam) from chMarkedPreseq.filter{ it[0][-6..-1] != '_spike'}

  output:
  file "*.ccurve.txt" into chPreseqStats

  script:
  defectMode = params.preseq_defect ? '-D' : ''
  """
  preseq lc_extrap -v $defectMode -output ${prefix}.ccurve.txt -bam ${bam[0]}
  """
}

/*
 * BAM Filtering
 */

process bamFiltering {
  tag "${prefix}"
  publishDir path: "${params.outdir}/filtering/filteredBams", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
             else if (filename.endsWith("_filtered.bam") || (filename.endsWith("_filtered.bam.bai"))) filename
             else null
            }

  when:
  !params.skipFiltering

  input:
  set val(prefix), file(markedBam) from chMarkedBamsFilt
  file bed from chGeneBed.collect()
  file bamtoolsFilterConfig from chBamtoolsFilterConfig.collect()

  output:
  set val(prefix), file("*.{bam,bam.bai}") into chFilteredBams
  set val(prefix), file("*.flagstat") into chFilteredFlagstat
  file "*.{idxstats,stats}" into chFilteredStats

  script:
  filterParams = params.singleEnd ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
  dupParams = params.keepDups ? "" : "-F 0x0400"
  mapQParams = params.mapQ ? "" : "-q ${params.mapQ}"
  blacklistParams = params.blacklist ? "-L $bed" : ""
  nameSortBam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
  """
  samtools view \\
    $filterParams \\
    $dupParams \\
    -F 0x08 \\
    $mapQParams \\
    $blacklistParams \\
    -b ${markedBam[0]} \\
    | bamtools filter \\
      -out ${prefix}_filtered.bam \\
      -script $bamtoolsFilterConfig
  samtools index ${prefix}_filtered.bam
  samtools flagstat ${prefix}_filtered.bam > ${prefix}_filtered.flagstat
  samtools idxstats ${prefix}_filtered.bam > ${prefix}_filtered.idxstats
  samtools stats ${prefix}_filtered.bam > ${prefix}_filtered.stats
  $nameSortBam
  """
}


// Preparing channels for all subsequent processes using filtered bams
//if(!params.singleEnd){
//  process mergeBamRemoveOrphan {
//    tag "$prefix"
//    publishDir path: "${params.outdir}/filtering/filteredBams/mergedLibrary", mode: 'copy',
//      saveAs: { filename ->
//              if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
//              else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
//              else if (filename.endsWith(".stats")) "samtools_stats/$filename"
//              else if (filename.endsWith("_sorted.bam")) filename
//              else if (filename.endsWith("_sorted.bam.bai")) filename
//              else null
//          }
//
//    input:
//    set val(prefix), file(bam) from chFilteredBams
//
//    output:
//    set val(prefix), file("*_sorted.{bam,bam.bai}") into chFilteredBamsFinal, chGroupBamNameFeatCounts
//    set val(prefix), file("*.flagstat") into chFilteredFlagstatFinal, chFilteredFlagstatSpikes, chFilteredFlagstatMqc
//    file "*.{idxstats,stats}" into chFilteredStats
//
//    script:
//    """
//    bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
//    samtools sort -@ $task.cpus -o ${prefix}_sorted.bam -T $prefix ${prefix}.bam
//    samtools index ${prefix}_sorted.bam
//    samtools flagstat ${prefix}_sorted.bam > ${prefix}_sorted.flagstat
//    samtools idxstats ${prefix}_sorted.bam > ${prefix}_sorted.idxstats
//    samtools stats ${prefix}_sorted.bam > ${prefix}_sorted.stats
//    """
//    }
//} else {
//}

if (!params.skipFiltering){
  chFilteredBams
    .set{ chBams }
  chFilteredFlagstat
    .set{ chFlagstat }
  chFilteredStats
    .set{ chStats }
}else{
  chMarkedBams
    .set{ chBams }
  chMarkedFlagstat
    .set{ chFlagstat }
  chFilteredStats
    .set{ chStats }
}


/*
 * Separate sample BAMs and spike BAMs
 */

chFlagstatChip=Channel.create()
chFlagstatSpikes=Channel.create()
chFlagstat.choice( chFlagstatSpikes, chFlagstatChip ) { it -> it[0][-5..-1] == 'Spike' ? 0 : 1 }

chBamsChip=Channel.create()
chBamsSpikes=Channel.create()
chBams.choice( chBamsSpikes, chBamsChip ) { it -> it[0][-5..-1] == 'Spike' ? 0 : 1 }

// Preparing all ChIP data for further analysis
chBamsChip
  .dump (tag:'cbams')
  .into { chBamsMetrics; chBamsMacs1; chBamsMacs2; chBamsPPQT;
          chBamsBigWig; chBamsBigWigSF; 
          chBamsDeeptoolsCorBam; chBamsDeeptoolsCorBai; chBamsDeeptoolsCorSample;
          chBamsForCount; chBamsCounts }

chFlagstatChip
  .into { chFlagstatMacs; chFlagstatMqc }

chStats
  .set { chStatsMqc }


/*
 * PhantomPeakQualTools QC
 */

process PPQT{
  tag "${prefix}"
  publishDir "${params.outdir}/ppqt", mode: "copy"

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

  script:
  """
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
  publishDir "${params.outdir}/bigWig", mode: "copy"

  input:
  set val(prefix), file(filteredBams) from chBamsBigWig

  output:
  set val(prefix), file('*.bigwig') into chBigWig

  script:
  """
  bamCoverage -b ${filteredBams[0]} -o ${prefix}.bigwig -p ${task.cpus}
  """
}

// With Spikes
if (params.spike){
  chBamsSpikes
    .into{chBamsSpikesBam;chBamsSpikesBai}
  process binSpikes {
    publishDir "${params.outdir}/bigwig/10kbins"

    input:
    file(allBams) from chBamsSpikesBam.map{it[1][0]}.collect()
    file(allBai) from chBamsSpikesBai.map{it[1][1]}.collect()
 
    output:
    file "*.tab" into chDeseq2tab

    script:
    """
    multiBamSummary bins -b $allBams --binSize 10000 -o results.npz --outRawCounts readCounts_10kbins.tab
    """
  }

  process getDeseqScaleFactor{
    publishDir "${params.outdir}/bigwig/scalefactors"

    input:
    file binsTab from chDeseq2tab

    output:
    file "*.sf" into chTabSF

    script:
    """
    getDESeqSF.R '${binsTab}'
    """
  }

  chTabSF
    .splitCsv(header:false, sep:',')
    .view()
    .map { row -> [row[0], row[1]]}
    .set{chScaleFactor}

  chBamsBigWigSF
    .combine(chScaleFactor)
    .view()
    .filter{it[0] == it[2]}
    .map { it -> it[0,1,3]}
    .set{chBigWigScaleFactor}

  process bigWigGenerationScaled{
    tag "${prefix}"
    publishDir "${params.outdir}/bigWig", mode: "copy"

    input:
    set val(prefix), file(filteredBams), val(normFactor) from chBigWigScaleFactor

    output:
    set val(prefix), file('*.bigwig') into chBigWigSF

    script:
    """
    bamCoverage -b ${filteredBams[0]} -o ${prefix}_spikes.bigwig -p ${task.cpus} --scaleFactor ${normFactor}
    """
  }
}


/*
 * DeepTools QC
 */

process deepToolsComputeMatrix{
  tag "${prefix}"
  publishDir "${params.outdir}/deepTools/singleBam", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  set val(prefix), file(bigwig) from chBigWig
  file geneBed from chGeneBedDeeptools.collect()

  output:
  set val(prefix), file("*.{gz,pdf}") into chDeeptoolsSingle
  set val(prefix), file("*_corrected.tab") into chDeeptoolsSingleMqc
  script:
  """
  computeMatrix scale-regions -R $geneBed -S ${bigwig} \\
                -o ${prefix}_matrix.mat.gz \\
                --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.gz \\
                --downstream 1000 --upstream 1000 --skipZeros\\
                -p ${task.cpus}

  plotProfile -m ${prefix}_matrix.mat.gz -o ${prefix}_bams_profile.pdf \\
        --outFileNameData ${prefix}.plotProfile.tab
  sed -e 's/.0\t/\t/g' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}_plotProfile_corrected.tab
  """
}

process deepToolsCorrelationQC{
  publishDir "${params.outdir}/deepTools/multipleBams"

  input:
  file(allBams) from chBamsDeeptoolsCorBam.map{it[1][0]}.collect()
  file(allBai) from chBamsDeeptoolsCorBai.map{it[1][1]}.collect()
  val (allPrefix) from chBamsDeeptoolsCorSample.map{it[0]}.collect()

  output:
  file "bams_correlation.pdf" into chDeeptoolsCorrel
  file "bams_coverage.pdf" into chDeeptoolsCoverage
  file "bams_fingerprint.pdf" into chDeeptoolsFingerprint
  file "bams_correlation.tab" into chDeeptoolsCorrelMqc
  file "bams_coverage_raw.txt" into chDeeptoolsCovMqc
  file "bams_fingerprint*" into chDeeptoolsFingerprintMqc

  when:
  allPrefix.size() > 2 && !params.skipDeepTools

  script:
  """
  multiBamSummary bins -b $allBams -o bams_summary.npz  -p ${task.cpus}
  plotCorrelation -in bams_summary.npz -o bams_correlation.pdf \\
          -c spearman -p heatmap -l $allPrefix \\
          --outFileCorMatrix bams_correlation.tab

  plotCoverage -b $allBams -o bams_coverage.pdf -p ${task.cpus} -l $allPrefix --outRawCounts bams_coverage_raw.txt
  plotFingerprint -b $allBams -plot bams_fingerprint.pdf -p ${task.cpus} -l $allPrefix --outRawCounts bams_fingerprint_raw.txt --outQualityMetrics bams_fingerprint_qmetrics.tab
  """
}


/*#########################################################################
  /!\ From this point, 'design' is mandatory /!\
###########################################################################*/

/*
 * Peak calling index build
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
    .into { chGroupBamMacsSharp; chGroupBamMacsBroad; chGroupBamMacsVeryBroad; chGroupBamFeatCounts}

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

/*
 * Peak calling 
 */

//MACS2 SHARP
process sharpMACS2{
  tag "${sampleID} - ${controlID}"
  publishDir path: "${params.outdir}/peakCalling/sharp", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
            }
 
  when:
  !params.skipPeakCalling

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsSharp
  file peakCountHeader from chPeakCountHeaderSharp.collect()
  file fripScoreHeader from chFripScoreHeaderSharp.collect()

  output:
  set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputSharp
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.narrowPeak") into chPeaksMacsSharp
  //set val(group), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.narrowPeak") into chMacsHomerSharp, chMacsQcSharp, chMacsIdrSharp
  //file "*igv.txt" into chMacs_igv_sharp
  file "*_mqc.tsv" into chMacsCountsSharp

  script:
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  peaktypeMacs = "narrowPeak"
  """
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    -f $format \\
    -g $params.macsGsize \\
    -n $sampleID \\
    --keep-dup all
  cat ${sampleID}_peaks.$peaktypeMacs | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.$peaktypeMacs -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  find * -type f -name "*.$peaktypeMacs" -exec echo -e "peakCalling/sharp/"{}"\\t0,0,178" \\; > ${sampleID}_peaks.igv.txt
  """
 }

//MACS2 BROAD
process broadMACS2{
  tag "${sampleID} - ${controlID}"
  publishDir path: "${params.outdir}/peakCalling/broad", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
            }
   
  when:
  !params.skipPeakCalling

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsBroad
  file peakCountHeader from chPeakCountHeaderBroad.collect()
  file fripScoreHeader from chFripScoreHeaderBroad.collect()

  output:
  set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputBroad
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chPeaksMacsBroad
  //set val(group), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.broadPeak") into chMacsHomerBroad, chMacsQcBroad, chMacsIdrBroad
  //file "*igv.txt" into chMacs_igv_broad
  file "*_mqc.tsv" into chMacsCountsBroad

  script:
  broad = "--broad --broad-cutoff ${params.broad_cutoff}"
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  peaktypeMacs = "broadPeak"
  """
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    ${broad} \\
    -f $format \\
    -g $params.macsGsize \\
    -n $sampleID \\
    --keep-dup all
    cat ${sampleID}_peaks.$peaktypeMacs | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
    READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.$peaktypeMacs -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
    find * -type f -name "*.$peaktypeMacs" -exec echo -e "peakCalling/broad/"{}"\\t0,0,178" \\; > ${sampleID}_peaks.igv.txt
  """
}

//EPIC2 VERY BROAD
process veryBroadEpic2{
  tag "${sampleID} - ${controlID}"
  publishDir path: "${params.outdir}/peakCalling/very_broad", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
            }
 
  when:
  !params.skipPeakCalling

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsVeryBroad
  file peakCountHeader from chPeakCountHeaderVeryBroad.collect()
  file fripScoreHeader from chFripScoreHeaderVeryBroad.collect()

  output:
  // set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacs_output_vbroad ===>> Mandatory data for MQC in Sharp&Broad MACS2
  set val(group), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chPeaksEpic
  //set val(group), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.broadPeak") into chMacsHomerVbroad, chMacsQcVbroad, chMacsIdrVbroad
  //file "*igv.txt" into chMacs_igv_vbroad
  file "*_mqc.tsv" into chMacsCountsVbroad

  script:
  peaktypeEpic = "broadPeak"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  """
  epic2 -t ${sampleBam[0]} \\
    ${ctrl} \\
    -gn ${params.genome} \\
    -kd -a \\
    -o ${sampleID}_peaks.$peaktypeEpic
  cat ${sampleID}_peaks.$peaktypeEpic | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.$peaktypeEpic -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  find * -type f -name "*.$peaktypeEpic" -exec echo -e "peakCalling/very_broad/"{}"\\t0,0,178" \\; > ${sampleID}_peaks.igv.txt
  """
}

// Join the results of all peaks callers
chPeaksMacsSharp
  .mix(chPeaksMacsBroad, chPeaksEpic)
  .into{ chPeaksHomer; chIDRpeaks; chIDRid }


/*
 * Peaks Annotation
 */

process peakAnnoHomer{
  tag "${sampleID}"
  publishDir path: "${params.outdir}/peak_annotation", mode: 'copy'

  when:
  !params.skipPeakAnno

  input:
  set val(group), val(replicate), val(peaktype), val(sampleID), file (peakfile) from chPeaksHomer
  file gtfFile from chGtfHomer.collect()
  file fastaFile from chFastaHomer.collect()

  output:
  file "*.txt" into chHomerPeakAnnotated

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
  publishDir path: "${params.outdir}/IDR", mode: 'copy'

  when:
  allPeaks.toList().size > 1 && !params.skipIDR

  input:
  set val(group), file(allPeaks) from chPeaksPerGroup

  output:
  file "*idrValues.txt" into chIdr
  file "*log.txt" into chMqcIdr

  script:
  peaktype = allPeaks[0].toString()
  peaktype = peaktype.substring(peaktype.lastIndexOf(".") + 1)
  """
  idr --samples ${allPeaks} \\
      --input-file-type  ${peaktype} \\
      -o ${group}_idrValues.txt \\
      -l ${group}_log.txt \\
      --plot
  """
}


//
// Peak calling & annotation QC
//
//if (!params.skipPeakQC && !params.noInput && params.design){
//  if (params.design){
//    chMacsQcSharp
//      .mix(chMacsQcBroad)
//      .set{chMacsQc}
//  }
//  process peakQC{
//    publishDir "${params.outdir}/peak_QC/", mode: 'copy'
//
//    input:
//    file peaks from chMacsQc.collect{ it[-1] }
//    file annotations from chHomerPeakAnnotated.collect()
//    file peakHeader from chPeakAnnotationHeader
//
//    output:
//    file "*.{txt,pdf}" into chMacsQcOutput
//    file "*.tsv" into chPeakMqc
//
//    script:
//    """
//    plot_macs_qc.r \\
//      -i ${peaks.join(',')} \\
//      -s ${peaks.join(',').replaceAll("_peaks.narrowPeak","").replaceAll("_peaks.broadPeak","")} \\
//      -o ./ \\
//      -p peak
//    plot_homer_annotatepeaks.r \\
//      -i ${annotations.join(',')} \\
//      -s ${annotations.join(',').replaceAll("_annotated_peaks.txt","")} \\
//      -o ./ \\
//      -p annotatePeaks
//    cat $peakHeader annotatePeaks.summary.txt > annotatedPeaks.summary_mqc.tsv
//    """
//  }
//} else if (!params.skipPeakQC && params.noInput && params.design){
//  process peakQCNoInput{
//    publishDir "${params.outdir}/peak_QC/", mode: 'copy'
//
//   input:
//    file peaks from chMacsQc.collect{ it[-1] }
//    file annotations from chHomerPeakAnnotated.collect()
//    file peakHeader from chPeakAnnotationHeader
//
//    output:
//    file "*.{txt,pdf}" into chMacsQcOutput
//    file "*.tsv" into chPeakMqc
//
//    script:
//    """
//    plot_macs_qc.r \\
//      -i ${peaks.join(',')} \\
//      -s ${peaks.join(',').replaceAll("_peaks.narrowPeak","").replaceAll("_peaks.broadPeak","")} \\
//      -o ./ \\
//      -p peak
//    plot_homer_annotatepeaks.r \\
//      -i ${annotations.join(',')} \\
//      -s ${annotations.join(',').replaceAll("_annotated_peaks.txt","")} \\
//      -o ./ \\
//      -p annotatePeaks
//    cat $peakHeader annotatePeaks.summary.txt > annotatedPeaks.summary_mqc.tsv
//    """
//  }
//}


/*
 * Feature counts


if (!params.skipFeatCounts){
  chGroupBamFeatCounts
    .map { it -> [ it[3], [ it[0], it[1], it[2] ] ] }
    .join(chGroupBamNameFeatCounts)
    .map { it -> [ it[1][0], it[1][1], it[1][2], it[2] ] }
    .groupTuple()
    .map { it -> [ it[0], it[3].flatten().sort() ] }
    .set { chGroupBamFeatCounts }

  process featureCounts{
    tag "${sampleName}"
    publishDir "${params.outdir}/featCounts/${sampleName}"

    input:
    set val(sampleName), file(bams) from chGroupBamFeatCounts
    file gtf from chGtfFeatCounts

    output:
    file "*_featCounts.txt" into chFeatCounts
    file "*_featCounts.txt.summary" into chFeatCountsMqc

    script:
    bamFiles = bams.findAll { it.toString().endsWith('.bam') }.sort()
    paramsPairedEnd = params.singleEnd ? '' : "-p"
    """
    featureCounts -T ${task.cpus} -a $gtf -o ${sampleName}_featCounts.txt $paramsPairedEnd ${bamFiles.join(' ')}
    """
  }
}
*/






/*
 * MultiQC
 */

process getSoftwareVersions{
  publishDir path: "${params.outdir}/software_versions", mode: "copy"

  when:
  !params.skipMultiQC

  output:
    file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  fastqc --version &> v_fastqc.txt
  multiqc --version &> v_multiqc.txt
  echo \$(bwa 2>&1) &> v_bwa.txt
  bowtie2 --version &> v_bowtie2.txt
  STAR --version &> v_star.txt
  samtools --version &> v_samtools.txt
  bedtools --version &> v_bedtools.txt
  echo \$(bamtools --version 2>&1) > v_bamtools.txt
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  preseq &> v_preseq.txt
  echo \$(plotFingerprint --version 2>&1) > v_deeptools.txt || true
  R --version &> v_R.txt
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  epic2 --version &> v_epic2.txt
  idr --version &> v_idr.txt
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
  section_href: 'https://gitlab.curie.fr/rnaseq'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}


process multiqc {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file metadata from chMetadata.ifEmpty([])
  file multiqcConfig from chMultiqcConfig
  file ('software_versions/*') from softwareVersionsYaml.collect()
  file ('workflow_summary/*') from workflowSummaryYaml.collect()

  file ('fastQC/*') from chFastqcResults.collect().ifEmpty([])

  file ('/filtering/stats/*') from chMarkedStats.collect().ifEmpty([])
  file ('/filtering/stats/*') from chMarkedPicstats.collect().ifEmpty([])
  file ('/filtering/stats/*') from chFlagstatMqc.collect().ifEmpty([])
  file ('/filtering/stats/*') from chStatsMqc.collect().ifEmpty([])

  file ('preseq/*') from chPreseqStats.collect().ifEmpty([])

  file ('ppqt/*') from chPpqtOutMqc.collect().ifEmpty([])
  file ('ppqt/*') from chPpqtCsvMqc.collect().ifEmpty([])

  file ('deepTools/*') from chDeeptoolsSingle.collect().ifEmpty([])
  file ('deepTools/*_corrected.tab') from chDeeptoolsSingleMqc.collect().ifEmpty([])
  file ("deepTools/bams_correlation.tab") from chDeeptoolsCorrelMqc.collect().ifEmpty([])
  file ("deepTools/bams_coverage_raw.txt") from chDeeptoolsCovMqc.collect().ifEmpty([])
  file ("deepTools/bams_fingerprint_*") from chDeeptoolsFingerprintMqc.collect().ifEmpty([])

  file ('peakCalling/sharp/*') from chMacsOutputSharp.collect().ifEmpty([])
  file ('peakCalling/broad/*') from chMacsOutputBroad.collect().ifEmpty([])
  file ('peakCalling/sharp/*') from chMacsCountsSharp.collect().ifEmpty([])
  file ('peakCalling/broad/*') from chMacsCountsBroad.collect().ifEmpty([])
  //file ('peakCalling/very_broad/*') from chMacsCountsVbroad.collect().ifEmpty([])

  //file('peak_QC/*') from chPeakMqc.collect().ifEmpty([])
  //file('featCounts/*') from chFeatCountsMqc.collect().ifEmpty([])

  output:
  file splan
  file "*multiqc_report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  metadata_opts = params.metadata ? "--metadata ${metadata}" : ""
  splan_opts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  modules_list = "-m custom_content -m fastqc -m preseq -m phantompeakqualtools -m deeptools -m macs2 -m homer -m featureCounts -m samtools -m picard"
  """
  mqc_header.py --name "ChIP-seq" --version ${workflow.manifest.version} ${metadata_opts} ${splan_opts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml $modules_list -c $multiqcConfig -s
  """
}


/* Creates a file at the end of workflow execution */
workflow.onComplete {

  /*pipeline_report.html*/

  def report_fields = [:]
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

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/oncomplete_template.txt")
  def txt_template = engine.createTemplate(tf).make(report_fields)
  def report_txt = txt_template.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/oncomplete_template.html")
  def html_template = engine.createTemplate(hf).make(report_fields)
  def report_html = html_template.toString()

  // Write summary e-mail HTML to a file
  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
  output_d.mkdirs()
  }
  def output_hf = new File( output_d, "pipeline_report.html" )
  output_hf.withWriter { w -> w << report_html }
  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << report_txt }

  /*oncomplete file*/
  File woc = new File("${params.outdir}/workflow.oncomplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'

  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")    println endWfSummary
  String execInfo = "Execution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
  woc.write(execInfo)

  /*]      = workflow.success     endSummary['exit status']  = workflow.exitStatus     endSummary['Error report'] = workflow.errorReport ?: '-' final logs*/
  if(workflow.success){
    log.info "[Chip-seq] Pipeline Complete"
  }else{
    log.info "[Chip-seq] FAILED: $workflow.runName"
  }
}
