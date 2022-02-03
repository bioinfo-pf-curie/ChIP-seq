#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020-2021
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
========================================================================================
                                ChIP-seq DSL2
========================================================================================
 ChIP-seq Analysis Pipeline.
 https://gitlab.curie.fr/data-analysis/chip-seq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Initialize lintedParams and paramsWithUsage
NFTools.welcome(workflow, params)

// Use lintedParams as default params object
paramsWithUsage = NFTools.readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
params.putAll(NFTools.lint(params, paramsWithUsage))

// Run name
customRunName = NFTools.checkRunName(workflow.runName, params.name)

// Custom functions/variables
mqcReport = []

include {checkSpikeAlignmentPercent} from './lib/functions'
include {loadDesign} from './lib/functions'

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

params.fasta = NFTools.getGenomeAttribute(params, 'fasta')
params.chrsize = NFTools.getGenomeAttribute(params, 'chrsize')
params.spikeFasta = NFTools.getGenomeAttribute(params, 'fasta', genome=params.spike)
params.bwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex')
params.spikeBwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex', genome=params.spike)
params.bowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2Index')
params.spikeBowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2Index', genome=params.spike)
params.starIndex = NFTools.getGenomeAttribute(params, 'starIndex')
params.spikeStarIndex = NFTools.getGenomeAttribute(params, 'starIndex', genome=params.spike)
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')
params.geneBed = NFTools.getGenomeAttribute(params, 'geneBed')
params.blacklist = NFTools.getGenomeAttribute(params, 'blacklist')
params.effGenomeSize = NFTools.getGenomeAttribute(params, 'effGenomeSize')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

/*
==========================
 VALIDATE INPUTS
==========================
*/

if (params.genomes && params.spike && !params.genomes.containsKey(params.spike)) {
  exit 1, "The provided spike-in genome '${params.spike}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.aligner != 'bwa-mem' && params.aligner != 'star' && params.aligner != 'bowtie2' ) {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'bowtie2' or 'bwa-mem'"
}

if (!params.design) {
  log.info "=================================================================\n" +
            "INFO: No design file detected.\n" +
            "Peak calling and annotation will be skipped.\n" +
            "Please set up a design file '--design' to run these steps.\n" +
            "================================================================"
}

if (!params.effGenomeSize) {
  log.warn "=================================================================\n" +
            "WARNING! Effective Genome Size is not defined.\n" +
            "Peak calling and annotation will be skipped.\n" +
            "Please specify value for '--effGenomeSize' to run these steps.\n" +
            "================================================================"
}


/*
==========================
 BUILD CHANNELS
==========================
*/

// Genome Fasta file
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set{chFasta}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

if ( params.spikeFasta ){
  Channel
    .fromPath(params.spikeFasta, checkIfExists: true)
    .set{chFastaSpike}
}else{
  chFastaSpike = Channel.empty()
}

/*
 * Indexes
 */

if( params.starIndex && params.aligner == 'star' ){
  Channel
    .fromPath(params.starIndex)
    .ifEmpty { exit 1, "STAR index not found: ${params.starIndex}" }
    .set {chStarIndex}
}
else if ( params.bowtie2Index && params.aligner == 'bowtie2' ){
  Channel
    .fromPath("${params.bowtie2Index}")
    .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2Index}" }
    .set{chBowtie2Index}
}
else if ( params.bwaIndex && params.aligner == "bwa-mem" ){
  Channel
    .fromPath("${params.bwaIndex}")
    .ifEmpty { exit 1, "Bwa index not found: ${params.bwaIndex}" }
    .set{chBwaIndex}
}
else {
    exit 1, "No genome index specified!"
}

/*
 * Spike Indexes
 */

if( params.spikeStarIndex && params.aligner == 'star' ){
  Channel
    .fromPath(params.spikeStarIndex)
    .ifEmpty { exit 1, "STAR index not found: ${params.spikeStarIndex}" }
    .set {chSpikeStarIndex}
}
else if ( params.spikeBowtie2Index && params.aligner == 'bowtie2' ){
  Channel
    .fromPath("${params.spikeBowtie2Index}")
    .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bowtie2Index}" }
    .set{chSpikeBowtie2Index}
}
else if ( params.spikeBwaIndex && params.aligner == "bwa-mem" ){
  Channel
    .fromPath("${params.spikeBwaIndex}")
    .ifEmpty { exit 1, "Bwa index not found: ${params.spikeBwaIndex}" }
    .set{chSpikeBwaIndex}
}else{
  chSpikeStarIndex = Channel.empty()
  chSpikeBowtie2Index = Channel.empty()
  chSpikeBwaIndex = Channel.empty()
}

/*
 * Annotations
 */

if (params.geneBed) {
  Channel
    .fromPath(params.geneBed, checkIfExists: true)
    .ifEmpty {exit 1, "BED file ${params.geneBed} not found"}
    .set{chGeneBed}
}else{
  chGeneBed = Channel.empty()
}

if (params.blacklist) { 
  Channel
    .fromPath(params.blacklist, checkIfExists: true)
    .set {chBlacklist} 
}else{
  chBlacklist = Channel.empty()
}

if ( params.chrsize ){
  Channel
    .fromPath(params.chrsize, checkIfExists: true)
    .set{chChromSize}
}else{
  exit 1, "Chromosome size file not found: ${params.chrsize}"
}

if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .set{chGtf}
}else {
  exit 1, "GTF annotation file not specified!"
}

if (params.effGenomeSize){
  Channel
    .of(params.effGenomeSize)
    .set{ chEffGenomeSize }
}

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Run Name': customRunName,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Single-end' : params.singleEnd ?: null,
  'Fragment Size' : params.singleEnd ? params.fragmentSize : null,
  'Design' : params.design ?: null,
  'Genome' : params.genome,
  'Spike' : params.spike ?: null,
  'GTF Annotation' : params.gtf ?: null,
  'BED Annotation' : params.geneBed ?: null,
  'Blacklist' : params.blacklist ?: null,
  'Aligner' : params.aligner ?: null,
  'Keep Dups' : params.keepDups ?: null,
  'Min MapQ' : params.mapq ?: null,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)


/*
==============================
  LOAD INPUT DATA
==============================
*/

// TODO - start from BAM files

// Load raw reads
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)

// Make samplePlan if not available
chSplan = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

// Design
if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .set { chDesignFile }

  chDesignControl = loadDesign(params.design)

}else{
  chDesignControl = Channel.empty()
  chDesignFile = Channel.empty()
}

/*
==================================
           INCLUDE
==================================
*/ 

// Workflows
include { mappingBwaMemFlow as mappingBwaMemFlow } from './nf-modules/local/subworkflow/mappingBwaMem'
include { mappingStarFlow as mappingStarFlow } from './nf-modules/local/subworkflow/mappingStar'
include { mappingBowtie2Flow as mappingBowtie2Flow } from './nf-modules/local/subworkflow/mappingBowtie2'
include { bamFilteringFlow as bamFilteringFlowRef } from './nf-modules/local/subworkflow/bamFiltering'
include { bamFilteringFlow as bamFilteringFlowSpike } from './nf-modules/local/subworkflow/bamFiltering'

include { bamChipFlow }        from './nf-modules/local/subworkflow/bamChip' 
include { bamSpikesFlow }      from './nf-modules/local/subworkflow/bamSpikes' 

// Peak calling
include { peakCallingFlow }     from './nf-modules/local/subworkflow/peakcalling' 

// Processes
include { checkDesign } from './nf-modules/local/process/checkDesign'
include { fastqc } from './nf-modules/common/process/fastqc'
include { prepareAnnotation }   from './nf-modules/local/process/prepareAnnotation'
include { preseq } from './nf-modules/common/process/preseq'

include { featureCounts }       from './nf-modules/local/process/featureCounts'
include { getSoftwareVersions } from './nf-modules/common/process/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/outputDocumentation'
include { multiqc }             from './nf-modules/local/process/multiqc'


workflow {
  chVersions = Channel.empty()

  main:
  // Init MultiQC Channels
  chFastqcMqc = Channel.empty()
  chAlignedBamMqc = Channel.empty()
  chCompareBamsMqc = Channel.empty()
  chPreseqMqc = Channel.empty()
  chPeaksOutput=Channel.empty()
  chPeaksCountsMqc=Channel.empty()
  chFripResults=Channel.empty()
  chPeaksQCMqc=Channel.empty()

  // subroutines
  outputDocumentation(
    chOutputDocs,
    chOutputDocsImages
  )

  // Check design
  if (params.design){
    checkDesign(
      chDesignFile, 
      chSplan)
  }

  // PROCESS: fastqc
  if (!params.skipFastqc){
    fastqc(
      chRawReads
    )
    chFastqcMqc = fastqc.out.results.collect()
    chVersions = chVersions.mix(fastqc.out.versions)
  }


  //*******************************************
  // MAPPING

  // SUBWORKFLOW: STAR mapping
  if (params.aligner == "star"){
    mappingStarFlow(
      chRawReads,
      chStarIndex,
      chSpikeStarIndex,
    )
    chAlignedBam = mappingStarFlow.out.bam
    chVersions = chVersions.mix(mappingStarFlow.out.versions)
  }

  // SUBWORKFLOW: Bowtie2 mapping
  if (params.aligner == "bowtie2"){
    mappingBowtie2Flow(
      chRawReads,
      chBowtie2Index,
      chSpikeBowtie2Index,
    )
    chAlignedBam = mappingBowtie2Flow.out.bam
    chVersions = chVersions.mix(mappingBowtie2Flow.out.versions)
  }

  // SUBWORKFLOW: Bwa mapping
  if (params.aligner == "bwa-mem"){
    mappingBwaMemFlow(
      chRawReads,
      chBwaIndex.collect(),
      chSpikeBwaIndex,
    )
    chAlignedBam = mappingBwaMemFlow.out.bam
    chAlignedBamMqc = mappingBwaMemFlow.out.logs
    chAlignedFlagstat = mappingBwaMemFlow.out.flagstat
    chAlignedSpikeBam = mappingBwaMemFlow.out.spikeBam
    chCompareBamsMqc = mappingBwaMemFlow.out.compareBamsMqc
    chVersions = chVersions.mix(mappingBwaMemFlow.out.versions)
  }

  // Remove low mapping rate for spikes
  chCompareBamsMqc.join(chAlignedSpikeBam)
    .filter { prefix, logs, bam, bai -> checkSpikeAlignmentPercent(prefix, logs, params.spikePercentFilter) }
    .map { prefix, logs, bam, bai -> [ prefix, bam, bai ] }
    .set { chPassedSpikeBam }


  //*******************************************
  // CHIP WORFLOW

  //if (params.inputBam){
  //  qcFlow.out.chFastqcMqc = Channel.empty()
  //  mappingFlow.out.chMappingMqc = Channel.empty()
  //}

  if (!params.skipSaturation){ 
    preseq(
      chAlignedBam
    )
   chPreseqMqc = preseq.out.results.collect()
  }

  bamFilteringFlowRef(
    chAlignedBam
  )
  chVersions = chVersions.mix(bamFilteringFlowRef.out.versions)

  bamChipFlow(
    bamFilteringFlowRef.out.bam,
    chBlacklist,
    chGeneBed
  )
  chVersions = chVersions.mix(bamChipFlow.out.versions)


  //*******************************************
  // SPIKE WORKFLOW

  if (params.spike){
    bamFilteringFlowSpike(
      chPassedSpikeBam
    )
    chVersions = chVersions.mix(bamFilteringFlowSpike.out.versions)
   
    bamSpikesFlow(
      bamFilteringFlowRef.out.bam,
      bamFilteringFlowSpike.out.bam,
      chBlacklist
    )
    chVersions = chVersions.mix(bamSpikesFlow.out.versions)
  }

  //*********************************************
  // DOWNSTREAM ANALYSIS (DESIGN IS MANDATORY !)
 
  if (params.design){

    peakCallingFlow(
      bamFilteringFlowRef.out.bam,
      chDesignControl,
      chEffGenomeSize,
      chChromSize,
      chGtf,
      chFasta
    )
    chVersions = chVersions.mix(peakCallingFlow.out.versions)
    chPeaksOutput = peakCallingFlow.out.peaksOutput
    chPeaksCountsMqc = peakCallingFlow.out.peaksCountsMqc
    chFripResults = peakCallingFlow.out.fripResults
    chPeaksQCMqc = peakCallingFlow.out.peaksQCMqc
    chTSSFeatCounts = prepareAnnotation(chGeneBed.collect())
 
    featureCounts(
      bamFilteringFlowRef.out.bam.map{it -> it[1]}.collect(),
      chGeneBed.concat(chTSSFeatCounts)
    )
    chVersions = chVersions.mix(featureCounts.out.versions)
  }

  //*******************************************
  // MULTIQC

  if (!params.skipMultiQC){

    getSoftwareVersions(
      chVersions.unique().collectFile()
    )

    // Warnings
    chAlignedSpikeBam
      .join(chPassedSpikeBam, remainder: true)
      .filter{it -> it[3] == null}
      .flatMap{ it -> it[0] + ": Poor spike alignment rate. Sample ignored !"}
      .set{chWarnMapping}

    chWarnMapping
      .collectFile(name: 'warnings.txt', newLine: true)
      .set{chWarn}
  
    multiqc(
      customRunName, 
      chSplan.collect(),
      chMetadata.ifEmpty([]),
      chMultiqcConfig, 
      chDesignFile.collect().ifEmpty([]), 
      chFastqcMqc.collect().ifEmpty([]),
      chAlignedBamMqc.collect().ifEmpty([]),
      chCompareBamsMqc.collect().ifEmpty([]),
      chAlignedFlagstat.map{it->it[1]}.collect().ifEmpty([]),
      bamFilteringFlowRef.out.markdupMetrics.collect().ifEmpty([]),
      bamFilteringFlowRef.out.flagstat.map{it->it[1]}.collect().ifEmpty([]),
      chPreseqMqc.collect().ifEmpty([]),
       //  sortingFlow.out.chStatsMqc.collect().ifEmpty([]),
      bamChipFlow.out.fragmentsSize.collect().ifEmpty([]),
      bamChipFlow.out.ppqtOutMqc.collect().ifEmpty([]), 
      bamChipFlow.out.ppqtCsvMqc.collect().ifEmpty([]),
      bamChipFlow.out.deeptoolsProfileMqc.collect().ifEmpty([]),
      bamChipFlow.out.deeptoolsCorrelateMqc.collect().ifEmpty([]),
      bamChipFlow.out.deeptoolsFingerprintMqc.collect().ifEmpty([]),
      chPeaksOutput.collect().ifEmpty([]),
      chPeaksCountsMqc.collect().ifEmpty([]),
      chFripResults.collect().ifEmpty([]),
      chPeaksQCMqc.collect().ifEmpty([]),
      getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
      workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
      chWarn.collect().ifEmpty([])
    )
  }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
