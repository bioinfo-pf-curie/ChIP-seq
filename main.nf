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
  // TODO: Add to this help message with new command line parameters

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
  --aligner                Alignment tool to use:
                              Available : bwa, star, bowtie2

  Options:
  --singleEnd              Specifies that the input is single end reads

  SpikeIn analysis
  --spike                  Indicates if the experiment includes a spike-in normalization.
                           Default : false. Available : 'spike' to use metagenome with reference genome
                           '[spike genome]' to use a specific second genome

  References           If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --fasta                  Path to Fasta reference
  Indexes              Path to the indexes for aligners
  --starIndex              Index for STAR aligner
  --bwaIndex               Index for BWA MEM aligner
  --bowtie2Index           Index for Bowtie2 aligner
  Annotation
  --bed                    BED annotation file. Used with samtools to filter the reads, and in Deeptools ComputeMatrix function
  --gtf                    GTF annotation file. Used in HOMER peak annotation
  Peak calling
  --macsGzise              Reference genome size for MACS2
  --noInput                Default : false. Use --noInput to indicate no input controls in peak calling tools.

  Skip options:        All are false by default
  --skipMultiqc            Skips final report writing
  --skipFastqc             Skips fastQC
  --skipAlignment          Skips alignments
  --skipPreseq             Skips preseq QC
  --skipFiltering          Skips filtering
  --skipPpqt               Skips phantompeakqualtools QC
  --skipDeepTools          Skips deeptools QC
  --skipPeakcalling        Skips peak calling
  --skipPeakanno           Skips peak annotation
  --skipPeakQC             Skips peak QC
  --skipIdr                Skips IDR QC
  --skipFeatCounts         Skips feature count

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

// TODO - Add any reference files that are needed - see igenome.conf
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
    .into{chFastaHomer;
          chFastaBwa;
          chFastaBt2;
          chFastaStar}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

if (params.spike && params.spike != 'spike'){
  params.spikeFasta = params.spike ? params.genomes[ params.spike ].fasta ?: false : false
  if ( params.spikeFasta ){
    Channel
      .fromPath(params.spikeFasta, checkIfExists: true)
      .into{chSpikeFaBwa;
            chSpikeFaBt2;
            chSpikeFaStar}
  } else {
    exit 1, "Spike-in fasta file not found: ${params.spikeFasta}"
  }
}

if (params.aligner == "bwa-mem"){
  params.bwaIndex = genomeRef ? params.genomes[ genomeRef ].bwaIndex ?: false : false
  if (params.bwaIndex){
    lastPath = params.fasta.lastIndexOf(File.separator)
    bwaDir =  params.bwaIndex.substring(0,lastPath+1)
    bwaBase = params.bwaIndex.substring(lastPath+1)
    Channel
      .fromPath(bwaDir, checkIfExists: true)
      .set { chBwaIndex }
  } else {
    exit 1, "BWA index file not found: ${params.bwaIndex}"
  }

  if (params.spike && params.spike != 'spike'){
    params.spikeBwaIndex = params.spike ? params.genomes[ params.spike ].bwaIndex ?: false : false
    if (params.spikeBwaIndex){
      lastPath = params.spikeFasta.lastIndexOf(File.separator)
      bwaDir =  params.spikeBwaIndex.substring(0,lastPath+1)
      spikeBwaBase = params.spikeBwaIndex.substring(lastPath+1)
      Channel
        .fromPath(bwaDir, checkIfExists: true)
        .set { chSpikeBwaIndex }
    } else {
      exit 1, "Spike BWA index file not found: ${params.spikeBwaIndex}"
    }
  }
}

if (params.aligner == "bowtie2"){
  params.bt2Index = genomeRef ? params.genomes[ genomeRef ].bowtie2Index ?: false : false
  if (params.bt2Index){
    lastPath = params.fasta.lastIndexOf(File.separator)
    bt2Dir =  params.bt2Index.substring(0,lastPath+1)
    bt2Base = params.bt2Index.substring(lastPath+1)
    Channel
      .fromPath(bt2Dir, checkIfExists: true)
      .set { chBt2Index }
  } else {
    exit 1, "Bowtie2 index file not found: ${params.bt2Index}"
  }

  if (params.spike && params.spike != 'spike'){
        params.spikeBt2Index = params.spike ? params.genomes[ params.spike ].bowtie2Index ?: false : false
    if (params.spikeBt2Index){
      lastPath = params.spikeFasta.lastIndexOf(File.separator)
      bt2Dir =  params.spikeBt2Index.substring(0,lastPath+1)
      spikeBt2Base = params.spikeBt2Index.substring(lastPath+1)
      Channel
        .fromPath(bt2Dir, checkIfExists: true)
        .set { chSpikeBt2Index }
    } else {
      exit 1, "Spike bowtie2 index file not found: ${params.spikeBt2Index}"
    }
  }
}

if (params.aligner == "star"){
  params.starIndex = genomeRef ? params.genomes[ genomeRef ].starIndex ?: false : false
  if (params.starIndex){
    lastPath = params.fasta.lastIndexOf(File.separator)
    starDir =  params.starIndex.substring(0,lastPath+1)
    starBase = params.starIndex.substring(lastPath+1)
    Channel
      .fromPath(starDir, checkIfExists: true)
      .set { chStarIndex }
  } else {
    exit 1, "STAR index file not found: ${params.starIndex}"
  }

  if (params.spike && params.spike != 'spike'){
    params.spikeStarIndex = params.spike ? params.genomes[ params.spike ].starIndex ?: false : false
    if (params.spikeStarIndex){
      lastPath = params.spike.lastIndexOf(File.separator)
      starDir =  params.spikeStarIndex.substring(0,lastPath+1)
      spikeStarBase = params.spikeStarIndex.substring(lastPath+1)
      Channel
        .fromPath(starDir, checkIfExists: true)
        .set { chSpikeStarIndex }
    } else {
      exit 1, "Spike STAR index file not found: ${params.spikeStarIndex}"
    }
  }
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
    .into{chGeneBed;
          chGeneBedDeeptools}
}

if (params.spike && params.spike != 'spike'){
  params.spikeGeneBed = params.spike ? params.genomes[ params.spike ].bed12 ?: false : false
  if (params.spikeGeneBed) {
    Channel
      .fromPath(params.spikeGeneBed, checkIfExists: true)
      .set{chSpikeGeneBed}
  }
}

params.macsGsize = genomeRef ? params.genomes[ genomeRef ].macsGsize ?: false : false
params.blacklist = genomeRef ? params.genomes[ genomeRef ].blacklist ?: false : false

//PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

//Peak Calling headers
chPeakCountHeader = file("$baseDir/assets/peak_count_header.txt", checkIfExists: true)
chFripScoreHeader = file("$baseDir/assets/frip_score_header.txt", checkIfExists: true)
chPeakAnnotationHeader = file("$baseDir/assets/peak_annotation_header.txt", checkIfExists: true)

// Stage config files
  Channel
    .fromPath(params.multiqcConfig, checkIfExists: true)
    .set{chMultiqcConfig}
  Channel
    .fromPath(params.outputDoc, checkIfExists: true)
    .set{chOutputDocs}
if (params.singleEnd) {
  chBamtoolsFilterConfig = Channel.fromPath(params.bamtoolsFilterSEConfig, checkIfExists: true)
} else {
  chBamtoolsFilterConfig = Channel.fromPath(params.bamtoolsFilterPEConfig, checkIfExists: true)
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
customRunName = workflow.runName
}

// Header log info
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
// TODO : Report custom parameters here
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
 * Retrieve sample plan & design control files
 */

if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .ifEmpty { exit 1, "Sample Plan file not found: ${params.samplePlan}" }
    .into { chSplan;chSplanCheck;chReadsToMap }
}


if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { chDesignCheck;chDesignControl }

  chDesignControl
    .splitCsv(header:false, sep:',')
    .map { row -> [ row[0], row[1], row[2], row[3], row[4] ] }
    .set { chDesignControl }
} else {
  Channel
    .from(params.design)
    .set {chDesignCheck}
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
 * Create a channel for input read files
 */

if(params.samplePlan){
  if(params.singleEnd){
    chReadsToMap
      .splitCsv(header: false, sep:',')
      .map{ row -> [ row[0], [file(row[2], checkIfExists: true)]] }
      .into { rawReadsFastqc;rawReadsBWA;rawReadsBt2;rawReadsSTAR;rawSpikeReadsBWA;rawSpikeReadsBt2;rawSpikeReadsSTAR }
  }
  else{
    chReadsToMap
      .splitCsv(header: false, sep:',')
      .map{ row -> [ row[0], [file(row[2], checkIfExists: true), file(row[3], checkIfExists: true)]] }
      .into { rawReadsFastqc;rawReadsBWA;rawReadsBt2;rawReadsSTAR;rawSpikeReadsBWA;rawSpikeReadsBt2;rawSpikeReadsSTAR }
  }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { rawReadsFastqc;rawReadsBWA;rawReadsBt2;rawReadsSTAR }
  }
  else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
      .into { rawReadsFastqc;rawReadsBWA;rawReadsBt2;rawReadsSTAR }
  }
}
else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { rawReadsFastqc;rawReadsBWA;rawReadsBt2;rawReadsSTAR }
}

/*
 * FastQC
 */

if (!params.skipFastqc){
  process fastQC{
    tag "${prefix}"
    publishDir "${params.outdir}/fastQC"

    input:
      set val(prefix), file(reads) from rawReadsFastqc

    output:
      file "*_fastqc.{zip,html}" into chFastqcResults

    script:
      """
      fastqc -q $reads -t 6
      """
  }
}

/*
 * Alignment on reference genome
 */

// BWA-MEM
if (!params.skipAlignment && params.aligner == "bwa-mem"){
  process BWA{
    tag "${prefix}"
    publishDir "${params.outdir}/alignment/bwa_alignment"

    input:
      set val(prefix), file(reads) from rawReadsBWA
      file bwaIndex from chBwaIndex.collect()

    output:
      set val(prefix), file("*.bam") into chAlignReads

    script:
    """
    bwa mem -t ${task.cpus} $bwaIndex${bwaBase} $reads \\
    | samtools view -b -h -o ${prefix}.bam
    """
  }
}

// BOWTIE2
if (!params.skipAlignment && params.aligner == "bowtie2"){
  process bowtie2{
    tag "${prefix}"
    publishDir "${params.outdir}/alignment/bowtie2_alignment"

    input:
      set val(prefix), file(reads) from rawReadsBt2
      file bt2Index from chBt2Index.collect()

    output:
      set val(prefix), file("*.bam") into chAlignReads

    script:
    if (params.singleEnd){
      readcommand = "-1 ${reads[0]}"
    }
    else{
      readcommand = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    if (params.spike == 'spike'){
      alnMult = '-k 3'
    }
    else{
      alnMult = ''
    }
    """
    bowtie2 -p ${task.cpus} $alnMult -x $bt2Index/${bt2Base} $readcommand \\
    | samtools view -b -h -o ${prefix}.bam
    """
  }
}

// STAR
if (!params.skipAlignment && params.aligner == "star"){
  process star{
    tag "${prefix}"
    publishDir "${params.outdir}/alignment/star_alignment"

    input:
      set val(prefix), file(reads) from rawReadsSTAR
      file starIndex from chStarIndex.collect()

    output:
      set val(prefix), file('*.bam') into chAlignReads

    script:
    """
    STAR --genomeDir $starIndex/${starBase} --runThreadN ${task.cpus} --readFilesIn $reads --outSAMtype BAM Unsorted --readFilesCommand zcat \\
    | samtools view -b -h -o ${prefix}.bam
    """
  }
}

/*
 * Alignment on spikein genome
 */
if (params.spike && params.spike != 'spike'){
  // BWA-MEM
  if (!params.skipAlignment && params.aligner == "bwa-mem"){
    process spikeBWA{
      tag "${prefix}"
      publishDir "${params.outdir}/alignment/spike/bwa_alignment"

      input:
        set val(prefix), file(reads) from rawSpikeReadsBWA
        file spikeBwaIndex from chSpikeBwaIndex.collect()

      output:
        set val(spikeprefix), file("*_spike.bam") into chSpikeAlignReads

      script:
      spikeprefix = "${prefix}_spike"
      """
      bwa mem -t ${task.cpus} $spikeBwaIndex${spikeBwaBase} $reads \\
      | samtools view -b -h -o ${spikeprefix}.bam
      """
    }
  }

  // BOWTIE2
  if (!params.skipAlignment && params.aligner == "bowtie2"){
    process spikeBowtie2{
      tag "${prefix}"
      publishDir "${params.outdir}/alignment/spike/bowtie2_alignment"

      input:
        set val(prefix), file(reads) from rawSpikeReadsBt2
        file spikeBt2Index from chSpikeBt2Index.collect()

      output:
        set val(spikeprefix), file("*.bam") into chSpikeAlignReads

      script:
      spikeprefix = "${prefix}_spike"
      if (params.singleEnd){
        readcommand = "-1 ${reads[0]}"
      }
      else{
        readcommand = "-1 ${reads[0]} -2 ${reads[1]}"
      }
      """
      bowtie2 -p ${task.cpus} -x $spikeBt2Index/${spikeBt2Base} $readcommand \\
      | samtools view -b -h -o ${spikeprefix}.bam
      """
    }
  }

  // STAR
  if (!params.skipAlignment && params.aligner == "star"){
    process spikeStar{
      tag "${prefix}"
      publishDir "${params.outdir}/alignment/spike/star_alignment"

      input:
        set val(prefix), file(reads) from rawSpikeReadsSTAR
        file spikeStarIndex from chSpikeStarIndex.collect()

      output:
        set val(spikeprefix), file('*.bam') into chSpikeAlignReads

      script:
      spikeprefix = "${prefix}_spike"
      """
      STAR --genomeDir $spikeStarIndex/${spikeStarBase} --runThreadN ${task.cpus} --readFilesIn $reads --outSAMtype BAM Unsorted --readFilesCommand zcat \\
      | samtools view -b -h -o ${spikeprefix}.bam
      """
    }
  }
}

// Merging, if necessary reference aligned reads and spike aligned reads
if(params.spike && params.spike != 'spike'){
  process compareRefSpike{
    tag "${prefixRef}"
    publishDir "${params.outdir}/alignment/compareRefSpike"

    input:
      set val(prefixRef), file(unsortedBamRef) from chAlignReads
      set val(prefixSpike), file(unsortedBamSpike) from chSpikeAlignReads

    output:
      set val(prefixRef), file('*_ref.bam') into chAlignedReads
      set val(prefixSpike), file('*_spike.bam') into chSpikeAlignedReads
      file ('*.txt') into chLogCompare

    script:
    """
    samtools sort $unsortedBamRef -n -T ${prefixRef} -o ${prefixRef}_sorted.bam
    samtools sort $unsortedBamSpike -n -T ${prefixSpike} -o ${prefixSpike}_sorted.bam
    compareAlignments.py -IR $unsortedBamRef -IS $unsortedBamSpike -OR ${prefixRef}.bam -OS ${prefixSpike}.bam -SE ${params.singleEnd}
    """
  }
  
  chAlignedReads
    .mix(chSpikeAlignedReads)
    .set{chAlignedReads}
} else if(params.spike && params.spike == 'spike'){
  process sepMetagenome{
    tag "${prefix}"
    publishDir "${params.outdir}/alignment/sepMetagenome"

    input:
      set val(prefix), file(unsortedBam) from chAlignReads

    output:
      set val(prefix), file('*_ref.bam') into chAlignedReads
      set val(prefix), file('*_spike.bam') into chSpikeAlignedReads
      file ('*.txt') into chLogSep

    script:
    """
    samtools sort $unsortedBam -n -T ${prefix} -o ${prefix}_sorted.bam
    sepMetagenome.py -I ${prefix}_sorted.bam -OR ${prefix}_ref.bam -OS ${prefix}_spike.bam -SE ${params.singleEnd}
    """
  }

  chAlignedReads
    .mix(chSpikeAlignedReads)
    .set{chAlignedReads}
} else {
  chAlignReads
    .set{chAlignedReads}
}

/*
 * Sorting BAM files
 */

if (!params.skipAlignment){
  process bamSort{
    tag "${name}"
    publishDir path: "${params.outdir}/filtering/sortedBams", mode: 'copy',
          saveAs: {filename ->
                if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
                else if (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) filename
                else null
              }

    input:
      set val(prefix), file(unsortedBam) from chAlignedReads

    output:
      set val(name), file('*.{bam,bam.bai}') into chSortBam
      file "*.{flagstat,idxstats,stats}" into chSortBamStats

    script:
    String unsortedBamsName = unsortedBam.toString()
    if (unsortedBamsName.contains('spike')){
    name = "${prefix}_spike"
    } else {
    name = "${prefix}"
    }
    """
    samtools sort $unsortedBam -T ${name} -o ${name}_sorted.bam
    samtools index ${name}_sorted.bam
    samtools flagstat ${name}_sorted.bam > ${name}_sorted.flagstat
    samtools idxstats ${name}_sorted.bam > ${name}_sorted.idxstats
    samtools stats ${name}_sorted.bam > ${name}_sorted.stats
    """
  }
}

/*
 * Marking duplicates & filtering
 */

// Marking
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
    set val(prefix), file("*.{bam,bam.bai}") into chMarkedBam, chMarkedPreseq
    file "*.{flagstat,idxstats,stats}" into chMarkedBamSamstats
    file "*.txt" into chMarkedBamPicstats

  script:
  bamFiles = sortedBams.findAll { it.toString().endsWith('.bam') }.sort()
  """
  picard -Xmx4g MarkDuplicates \\
    INPUT=${bamFiles[0]} \\
    OUTPUT=${prefix}_marked.bam \\
    ASSUME_SORTED=true \\
    REMOVE_DUPLICATES=false \\
    METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR=tmpstararked.bam
  samtools index ${prefix}_marked.bam
  samtools idxstats ${prefix}_marked.bam > ${prefix}_marked.idxstats
  samtools flagstat ${prefix}_marked.bam > ${prefix}_marked.flagstat
  samtools stats ${prefix}_marked.bam > ${prefix}_marked.stats
  """
}

// Preseq complexity analysis before filtering
if(!params.skipPreseq) {
  chMarkedPreseq
    .filter{ it[0][-6..-1] != '_spike'}
    .set{chMarkedPreseq}

  process preseqAnalysis{
    tag "${prefix}"
    publishDir "${params.outdir}/preseq"

    input:
    set val(prefix), file(bam) from chMarkedPreseq

    output:
    file "*.ccurve.txt" into chPreseqStats

    script:
    if (params.preseq_defect){
      defectMode = '-D'
    }
    else{
      defectMode = ''
    }
    """
    preseq lc_extrap -v $defectMode -output ${prefix}.ccurve.txt -bam ${bam[0]}
    """
  }
}

// Filtering
if (!params.skipFiltering){
  process bamFiltering{
    tag "${prefix}"
    publishDir path: "${params.outdir}/filtering/filteredBams", mode: 'copy',
          saveAs: {filename ->
              if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "samtools_stats/$filename"
              else if (filename.endsWith("_filtered.bam") || (filename.endsWith("_filtered.bam.bai"))) filename
              else null
            }

    input:
    set val(prefix), file(markedBam) from chMarkedBam
    file bed from chGeneBed.collect()
    file bamtoolsFilterConfig from chBamtoolsFilterConfig.collect()

    output:
    set val(prefix), file("*.{bam,bam.bai}") into chFilteredBams
    set val(prefix), file("*.flagstat") into chFilteredFlagstat
    file "*.{idxstats,stats}" into chFilteredStats

    script:
    filterParams = params.singleEnd ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
    dupParams = params.keep_dups ? "" : "-F 0x0400"
    multimapParams = params.keep_multi_map ? "" : "-q 1"
    blacklistParams = params.blacklist ? "-L $bed" : ""
    nameSortBam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
    """
    samtools view \\
      $filterParams \\
      $dupParams \\
      -F 0x08 \\
      $multimapParams \\
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
}

// Preparing channels for all subsequent processes using filtered bams
if(!params.singleEnd){
  process mergeBamRemoveOrphan {
    tag "$prefix"
    publishDir path: "${params.outdir}/filtering/filteredBams/mergedLibrary", mode: 'copy',
      saveAs: { filename ->
              if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
              else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
              else if (filename.endsWith(".stats")) "samtools_stats/$filename"
              else if (filename.endsWith("_sorted.bam")) filename
              else if (filename.endsWith("_sorted.bam.bai")) filename
              else null
          }

    input:
    set val(prefix), file(bam) from chFilteredBams

    output:
    set val(prefix), file("*_sorted.{bam,bam.bai}") into chFilteredBamsFinal, chGroupBamNameFeatCounts
    set val(prefix), file("*.flagstat") into chFilteredFlagstatFinal, chFilteredFlagstatSpikes, chFilteredFlagstatMqc
    file "*.{idxstats,stats}" into chFilteredStatsFinal

    script:
    """
    bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
    samtools sort -@ $task.cpus -o ${prefix}_sorted.bam -T $prefix ${prefix}.bam
    samtools index ${prefix}_sorted.bam
    samtools flagstat ${prefix}_sorted.bam > ${prefix}_sorted.flagstat
    samtools idxstats ${prefix}_sorted.bam > ${prefix}_sorted.idxstats
    samtools stats ${prefix}_sorted.bam > ${prefix}_sorted.stats
    """
    }
} else {
  chFilteredBams
    .into{chFilteredBamsFinal;
          chGroupBamNameFeatCounts;
          chFilteredBamsSpikes}
  chFilteredFlagstat
    .into{chFilteredFlagstatFinal;
          chFilteredFlagstatSpikes}
  chFilteredStats
    .set{chFilteredStatsFinal}
}

// When using spikes, separating spike alignments from reference ones after filtering
if (params.spike){
  chFilteredFlagstatSpikes
    .filter{it[0][-6..-1] == '_spike'}
    .set{chFilteredFlagstatSpikes}

  chFilteredBamsSpikes
    .filter{it[0][-6..-1] == '_spike'}
    .set{chFilteredBamsSpikes}

  process genBinsScaleFactor{
    cache 'deep'
    publishDir "${params.outdir}/bigwig/10kbins"

    input:
    val filteredBamsAll from chFilteredBamsSpikes.toList()

    output:
    file '*.tab' into chDeseq2tab

    script:
    allBams=""
    for (List filteredBam : filteredBamsAll){
      allBams+=filteredBam[1][0]
      allBams+=" "
    }
    """
    multiBamSummary bins -b $allBams --binSize 10000 -o results.npz --outRawCounts readCounts_10kbins.tab
    """
  }

  process deseq2ScaleFactor{
    cache 'deep'
    publishDir "${params.outdir}/bigwig/scalefactors"

    input:
    file binsTab from chDeseq2tab

    output:
    file "*.sf" into chScaleFactor

    script:
    """
    getDESeqSF.R '${binsTab}'
    """
  }
}

// Preparing all filtered aligned reads for further analysis
chFilteredBamsFinal
  .filter{ it[0][-6..-1] != '_spike' }
  .into { chFilteredBamsMetrics;
      chFilteredBamsMacs1;
      chFilteredBamsMacs2;
      chFilteredBamsPhantompeakqualtools;
      chFilteredBamsDeeptoolsSingle;
      chFilteredBamsDeeptoolsCorrel;
      chFilteredBamsDeseq2;
      chFilteredBamsCounts }

chFilteredFlagstatFinal
  .into { chFilteredFlagstatMacs;
          chFilteredFlagstatMqc }

chFilteredStatsFinal
  .set { chFilteredStatsMqc }


/*
 * PhantomPeakQualTools QC
 */

if (!params.skipPpqt){
  process PPQT{
    tag "${prefix}"
    cache 'deep'
    publishDir "${params.outdir}/ppqt", mode: "copy"

    input:
    set val(prefix), file(filteredBams) from chFilteredBamsPhantompeakqualtools
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
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${filteredBams[0]}" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
    cp $sppCorrelationHeader ${prefix}_spp_correlation_mqc.tsv
    Rscript -e "load('${prefix}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${prefix}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
    awk -v OFS='\t' '{print "${prefix}", \$9}' ${prefix}.spp.out | cat $sppNSCHeader - > ${prefix}_spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${prefix}", \$10}' ${prefix}.spp.out | cat $sppRSCHeader - > ${prefix}_spp_rsc_mqc.tsv
    """
  }
}

/*
 * BigWig generation
 */
if (params.spike){
  chScaleFactor
    .splitCsv(header:false, sep:',')
    .map { row -> [row[0], row[1]]}
    .set{chScaleFactor}

  chFilteredBamsDeeptoolsSingle
    .filter{it[0][-6..-1] != '_spike'}
    .combine(chScaleFactor)
    .filter{it[0] == it[2]}
    .map { it -> it[0,1,3]}
    .set{chFilteredBamsDeeptoolsSingle}

  process bigWigGenerationScaled{
    tag "${prefix}"
    cache 'deep'
    publishDir "${params.outdir}/bigWig", mode: "copy"

    input:
    set val(prefix), file(filteredBams), val(normFactor) from chFilteredBamsDeeptoolsSingle

    output:
    set val(prefix), file('*_bw.bigwig') into chBigWig

    script:
    """
    bamCoverage -b ${filteredBams[0]} -o ${prefix}_bw.bigwig -p ${task.cpus} --scaleFactor ${normFactor}
    """
  }
} else {
  process bigWigGeneration{
    tag "${prefix}"
    cache 'deep'
    publishDir "${params.outdir}/bigWig", mode: "copy"

    input:
    set val(prefix), file(filteredBams) from chFilteredBamsDeeptoolsSingle

    output:
    set val(prefix), file('*_bw.bigwig') into chBigWig

    script:
    """
    bamCoverage -b ${filteredBams[0]} -o ${prefix}_bw.bigwig -p ${task.cpus}
    """
  }
}

/*
 * DeepTools QC
 */

if (!params.skipDeepTools){
  process deepToolsSingleQC{
    tag "${prefix}"
    cache 'deep'
    publishDir "${params.outdir}/deepTools/singleBam", mode: "copy"

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

  process deepToolsCorrelQC{
    cache 'deep'
    publishDir "${params.outdir}/deepTools/multipleBams"

    input:
    val filteredBamsAll from chFilteredBamsDeeptoolsCorrel.toList()

    output:
    file "bams_correlation.pdf" into chDeeptoolsCorrel
    file "bams_coverage.pdf" into chDeeptoolsCoverage
    file "bams_fingerprint.pdf" into chDeeptoolsFingerprint
    file "bams_correlation.tab" into chDeeptoolsCorrelMqc
    file "bams_coverage_raw.txt" into chDeeptoolsCovMqc
    file "bams_fingerprint*" into chDeeptoolsFingerprintMqc

    script:
    allBams=""
    allPrefixes=""
    for (List filteredBam : filteredBamsAll){
      allPrefixes+=filteredBam[0]
      allPrefixes+= " "
      allBams+=filteredBam[1][0]
      allBams+=" "
    }
    """
    multiBamSummary bins -b $allBams -o bams_summary.npz  -p ${task.cpus}
    plotCorrelation -in bams_summary.npz -o bams_correlation.pdf \\
            -c spearman -p heatmap -l $allPrefixes \\
            --outFileCorMatrix bams_correlation.tab

    plotCoverage -b $allBams -o bams_coverage.pdf -p ${task.cpus} -l $allPrefixes  --outRawCounts bams_coverage_raw.txt
    plotFingerprint -b $allBams -plot bams_fingerprint.pdf -p ${task.cpus} -l $allPrefixes --outRawCounts bams_fingerprint_raw.txt --outQualityMetrics bams_fingerprint_qmetrics.tab
    """
  }
}

/*
 * Peak calling index build
 */

if (!params.noInput && params.design){
  chFilteredBamsMacs1
    .join(chFilteredFlagstatMacs)
    .combine(chFilteredBamsMacs2)
    .set { chFilteredBamsMacs1 }

  chDesignControl
    .combine(chFilteredBamsMacs1)
    .filter { it[0] == it[5] && it[1] == it[8] }
    .map { it ->  it[2..-1] }
    .into { chGroupBamMacsSharp;
            chGroupBamMacsBroad;
            chGroupBamMacsVeryBroad;
            chGroupBamFeatCounts;
          }
} else if (params.noInput && params.design) {
  chDesignControl
    .combine(chFilteredBamsMacs1)
    .filter { it[0] == it[5]}
    .join(chFilteredFlagstatMacs)
    .map { it ->  it[2..-1] }
    .into { chGroupBamMacsSharp;
            chGroupBamMacsBroad;
            chGroupBamMacsVeryBroad;
            chGroupBamFeatCounts
          }
}

if (params.design){
  chGroupBamMacsSharp
    .filter { it[2] == 'sharp' }
    .set { chGroupBamMacsSharp }

  chGroupBamMacsBroad
    .filter { it[2] == 'broad' }
    .set { chGroupBamMacsBroad }

  chGroupBamMacsVeryBroad
    .filter { it[2] == 'very-broad' }
    .set { chGroupBamMacsVeryBroad }
}
/*
 * Peak calling WITH Input
 */

if (!params.skipPeakcalling && !params.noInput && params.design){
  // SHARP PEAKS
  process sharpMACS2{
    tag "${sampleID} - ${controlID}"
    publishDir path: "${params.outdir}/peakCalling/sharp", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsSharp
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputSharp
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.narrowPeak") into chMacsHomerSharp,
                                                                                                                 chMacsQcSharp,
                                                                                                                 chMacsIdrSharp
    //file "*igv.txt" into chMacs_igv_sharp
    file "*_mqc.tsv" into chMacsCountsSharp

    script:
    format = params.singleEnd ? "BAM" : "BAMPE"
    peaktypeMacs = "narrowPeak"
    """
    macs2 callpeak \\
      -t ${sampleBam[0]} \\
      -c ${controlBam[0]} \\
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

  // BROAD PEAKS
  process broadMACS2{
    tag "${sampleID} - ${controlID}"
    publishDir path: "${params.outdir}/peakCalling/broad", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    when:
    params.design

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsBroad
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputBroad
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.broadPeak") into chMacsHomerBroad,
                                                                                                                chMacsQcBroad,
                                                                                                                chMacsIdrBroad
    //file "*igv.txt" into chMacs_igv_broad
    file "*_mqc.tsv" into chMacsCountsBroad

    script:
    broad = "--broad --broad-cutoff ${params.broad_cutoff}"
    format = params.singleEnd ? "BAM" : "BAMPE"
    peaktypeMacs = "broadPeak"
    """
    macs2 callpeak \\
      -t ${sampleBam[0]} \\
      -c ${controlBam[0]} \\
      $broad \\
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

  // VERY BROAD PEAKS
  process veryBroadEpic2{
    tag "${sampleID} - ${controlID}"
    publishDir path: "${params.outdir}/peakCalling/very_broad", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    when:
    params.design

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat), val(controlID), file(controlBam) from chGroupBamMacsVeryBroad
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    // set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacs_output_vbroad ===>> Mandatory data for MQC in Sharp&Broad MACS2
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), val(controlID), file("*.broadPeak") into chMacsHomerVbroad,
                                                                                                                chMacsQcVbroad,
                                                                                                                chMacsIdrVbroad
    //file "*igv.txt" into chMacs_igv_vbroad
    file "*_mqc.tsv" into chMacsCountsVbroad

    script:
    peaktypeEpic = "broadPeak"
    """
    epic2 -t ${sampleBam[0]} \\
      -c ${controlBam[0]} \\
      -gn ${params.genome} \\
      -kd -a \\
      -o ${sampleID}_peaks.$peaktypeEpic
    cat ${sampleID}_peaks.$peaktypeEpic | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
    READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.$peaktypeEpic -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
    find * -type f -name "*.$peaktypeEpic" -exec echo -e "peakCalling/very_broad/"{}"\\t0,0,178" \\; > ${sampleID}_peaks.igv.txt
    """
  }
} else if (!params.skipPeakcalling && params.noInput && params.design){
  // SHARP PEAKS
  process sharpMACS2NoInput{
    tag "${sampleID}"
    publishDir path: "${params.outdir}/peakCalling/sharp", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    when:
    params.design

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat) from chGroupBamMacsSharp
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputSharp
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file("*.narrowPeak") into chMacsHomerSharp,
                                                                                                 chMacsQcSharp,
                                                                                                 chMacsIdrSharp
    //file "*igv.txt" into chMacs_igv_sharp
    file "*_mqc.tsv" into chMacsCountsSharp

    script:
    format = params.singleEnd ? "BAM" : "BAMPE"
    peaktypeMacs = "narrowPeak"
    """
    macs2 callpeak \\
      -t ${sampleBam[0]} \\
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

  // BROAD PEAKS
  process broadMACS2NoInput{
    tag "${sampleID}"
    publishDir path: "${params.outdir}/peakCalling/broad", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    when:
    params.design

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat) from chGroupBamMacsBroad
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacsOutputBroad
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chMacsHomerBroad,
                                                                                                chMacsQcBroad,
                                                                                                chMacsIdrBroad
    //file "*igv.txt" into chMacs_igv_broad
    file "*_mqc.tsv" into chMacsCountsBroad

    script:
    broad = "--broad --broad-cutoff ${params.broad_cutoff}"
    format = params.singleEnd ? "BAM" : "BAMPE"
    peaktypeMacs = "broadPeak"
    """
    macs2 callpeak \\
      -t ${sampleBam[0]} \\
      $broad \\
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

  // VERY BROAD PEAKS
  process veryBroadEpic2NoInput{
    tag "${sampleID}"
    publishDir path: "${params.outdir}/peakCalling/very_broad", mode: 'copy',
          saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else if (filename.endsWith(".igv.txt")) null
            else filename
          }
    cache 'deep'

    when:
    params.design

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file(sampleBam), file(sampleFlagstat) from chGroupBamMacsVeryBroad
    file peakCountHeader from chPeakCountHeader
    file fripScoreHeader from chFripScoreHeader

    output:
    // set val(sampleID), file("*.{bed,xls,gappedPeak,bdg}") into chMacs_output_vbroad ===>> Mandatory data for MQC in Sharp&Broad MACS2
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file("*.broadPeak") into chMacsHomerVbroad,
                                                                                                chMacsQcVbroad,
                                                                                                chMacsIdrVbroad
    //file "*igv.txt" into chMacs_igv_vbroad
    file "*_mqc.tsv" into chMacsCountsVbroad

    script:
    peaktypeEpic = "broadPeak"

    """
    epic2 -t ${sampleBam[0]} \\
      -gn ${params.genome} \\
      -kd -a \\
      -o ${sampleID}_peaks.$peaktypeEpic
    cat ${sampleID}_peaks.$peaktypeEpic | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
    READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.$peaktypeEpic -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
    find * -type f -name "*.$peaktypeEpic" -exec echo -e "peakCalling/very_broad/"{}"\\t0,0,178" \\; > ${sampleID}_peaks.igv.txt
    """
  }
}

/*
 * Irreproducible Discovery Rate
 */

if (!params.skipIdr && params.replicates  && params.design){
  chMacsIdrSharp
      .mix(chMacsIdrBroad, chMacsIdrVbroad)
      .into{chIdrSName;
        chIdrPeakfile}

  chIdrSName
      .map{ it -> it[0]}
      .unique()
      .set{chIdrSName}

  chIdrPeakfile
      .map{ it -> it[0, -1]}
      .set{chIdrPeakfile}

  process IDR{
    publishDir path: "${params.outdir}/IDR", mode: 'copy'

    input:
    val sampleNameAll from chIdrSName.toList()
    val peakFilesAll from chIdrPeakfile.toList()

    output:
    file "*.{narrowPeak,broadPeak}" into chIdr
    file "*_log.txt" into chMqcIdr

    script:
    """
    replicate_idr.py -sn ${sampleNameAll} -pf ${peakFilesAll}
    """

  }
}

/*
 * Peak annotation
 */

// HOMER


if (!params.skipPeakanno && !params.noInput && params.design){
  if (params.design){
  chMacsHomerSharp
    .mix(chMacsHomerBroad, chMacsHomerVbroad)
    .set{chMacs_homer}
  }

  process peakAnnoHomer{
    tag "${sampleID} - ${controlID}"
    publishDir path: "${params.outdir}/peak_annotation", mode: 'copy'

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), val(controlID), file (peakfile) from chMacs_homer
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
} else if (!params.skipPeakanno && params.noInput && params.design){
  process peakAnnoHomerNoInput{
    tag "${sampleID}"
    publishDir path: "${params.outdir}/peak_annotation", mode: 'copy'

    input:
    set val(sampleName), val(replicate), val(peaktype), val(sampleID), file (peakfile) from chMacs_homer
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
}

/*
 * Peak calling & annotation QC
 */

if (!params.skipPeakQC && !params.noInput && params.design){
  if (params.design){
    chMacsQcSharp
      .mix(chMacsQcBroad)
      .set{chMacsQc}
  }

  process peakQC{
    publishDir "${params.outdir}/peak_QC/", mode: 'copy'

    input:
    file peaks from chMacsQc.collect{ it[-1] }
    file annotations from chHomerPeakAnnotated.collect()
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
} else if (!params.skipPeakQC && params.noInput && params.design){
  process peakQCNoInput{
    publishDir "${params.outdir}/peak_QC/", mode: 'copy'

    input:
    file peaks from chMacsQc.collect{ it[-1] }
    file annotations from chHomerPeakAnnotated.collect()
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
}

/*
 * Feature counts
 */

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

/*
 * MultiQC
 */

// Retrieve software from environment
process getSoftwareVersions{
  publishDir path: "${params.outdir}/software_versions", mode: "copy"
  when:
    !params.skipMultiqc

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

// Workflow summary
process workflowSummaryMqc {
  when:
    !params.skipMultiqc

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
  !params.skipMultiqc

  input:
    file splan from chSplan.collect()
    file metadata from chMetadata.ifEmpty([])
    file multiqcConfig from chMultiqcConfig
    file ('software_versions/*') from softwareVersionsYaml.collect()
    file ('workflow_summary/*') from workflowSummaryYaml.collect()

    file ('fastQC/*') from chFastqcResults.collect().ifEmpty([])

    file ('filtering/markedBams/samtools_stats/*') from chMarkedBamSamstats.collect().ifEmpty([])
    file ('filtering/markedBams/samtools_stats/*') from chMarkedBamPicstats.collect().ifEmpty([])
    file ('filtering/filteredBams/samtools_stats/*') from chFilteredFlagstatMqc.collect().ifEmpty([])
    file ('filtering/filteredBams/samtools_stats/*') from chFilteredStatsMqc.collect().ifEmpty([])

    file ('preseq/*') from chPreseqStats.collect().ifEmpty([])

    file ('ppqt/*.spp.out') from chPpqtOutMqc.collect().ifEmpty([])
    file ('ppqt/*_mqc.tsv') from chPpqtCsvMqc.collect().ifEmpty([])

    file ('deepTools/singleBam/*') from chDeeptoolsSingle.collect().ifEmpty([])
    file ('deepTools/singleBam/*_corrected.tab') from chDeeptoolsSingleMqc.collect().ifEmpty([])
    file ("deepTools/multipleBams/bams_correlation.tab") from chDeeptoolsCorrelMqc.collect().ifEmpty([])
    file ("deepTools/multipleBams/bams_coverage_raw.txt") from chDeeptoolsCovMqc.collect().ifEmpty([])
    file ("deepTools/multipleBams/bams_fingerprint_*") from chDeeptoolsFingerprintMqc.collect().ifEmpty([])

    file ('peakCalling/sharp/*.xls') from chMacsOutputSharp.collect().ifEmpty([])
    file ('peakCalling/broad/*.xls') from chMacsOutputBroad.collect().ifEmpty([])
    file ('peakCalling/sharp/*') from chMacsCountsSharp.collect().ifEmpty([])
    file ('peakCalling/broad/*') from chMacsCountsBroad.collect().ifEmpty([])
    // file ('peakCalling/very_broad/*') from chMacsCountsVbroad.collect().ifEmpty([])

    file('peak_QC/*') from chPeakMqc.collect().ifEmpty([])

    file('featCounts/*') from chFeatCountsMqc.collect().ifEmpty([])


  output:
    file splan
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  metadata_opts = params.metadata ? "--metadata ${metadata}" : ""
  splan_opts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  modules_list = "-m custom_content -m fastqc -m samtools -m picard -m preseq -m phantompeakqualtools -m deeptools -m macs2 -m homer -m featureCounts"

  """
  mqc_header.py --name "Chip-seq" --version ${workflow.manifest.version} ${metadata_opts} > multiqc-config-header.yaml
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
