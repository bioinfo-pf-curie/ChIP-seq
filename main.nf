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
params.spikeFasta = NFTools.getGenomeAttribute(params, 'fasta', params.spike)
params.bwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex')
params.spikeBwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex', params.spike)
params.bowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2Index')
params.spikeBowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2Index', params.spike)
params.starIndex = NFTools.getGenomeAttribute(params, 'starIndex')
params.spikeStarIndex = NFTools.getGenomeAttribute(params, 'starIndex', params.spike)
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')
params.geneBed = NFTools.getGenomeAttribute(params, 'geneBed')
params.blacklist = NFTools.getGenomeAttribute(params, 'blacklist')
params.effGenomeSize = NFTools.getGenomeAttribute(params, 'effGenomeSize')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$projectDir/docs/output.md")
chOutputDocsImages = file("$projectDir/docs/images/", checkIfExists: true)

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

chFasta              = params.fasta                 ? Channel.fromPath(params.fasta, checkIfExists: true).collect()                  : Channel.empty()
chFastaSpike         = params.spikeFasta            ? Channel.fromPath(params.spikeFasta, checkIfExists: true).collect()             : Channel.empty()
chChromSize          = params.chrsize               ? Channel.fromPath(params.chrsize, checkIfExists: true).collect()                : Channel.empty()
chEffGenomeSize      = params.effGenomeSize         ? Channel.of(params.effGenomeSize)                                               : Channel.empty()

chStarIndex          = params.starIndex             ? Channel.fromPath(params.starIndex, checkIfExists: true).collect()              : Channel.empty()
chBowtie2Index       = params.bowtie2Index          ? Channel.fromPath(params.bowtie2Index, checkIfExists: true).collect()           : Channel.empty()
chBwaIndex           = params.bwaIndex              ? Channel.fromPath(params.bwaIndex, checkIfExists: true).collect()               : Channel.empty()
chMappingIndex       = params.aligner == 'star'     ? chStarIndex : params.aligner == 'bowtie2' ? chBowtie2Index : chBwaIndex

chSpikeStarIndex     = params.spikeStarIndex        ? Channel.fromPath(params.spikeStarIndex, checkIfExists: true).collect()         : Channel.empty()
chSpikeBowtie2Index  = params.spikeBowtie2Index     ? Channel.fromPath(params.spikeBowtie2Index, checkIfExists: true).collect()      : Channel.empty()
chSpikeBwaIndex      = params.spikeBwaIndex         ? Channel.fromPath(params.spikeBwaIndex, checkIfExists: true).collect()          : Channel.empty()
chSpikeIndex         = params.aligner == 'star'     ? chSpikeStarIndex : params.aligner == 'bowtie2' ? chSpikeBowtie2Index : chSpikeBwaIndex

chGeneBed            = params.geneBed               ? Channel.fromPath(params.geneBed, checkIfExists: true).collect()                : Channel.empty()
chBlacklist          = params.blacklist             ? Channel.fromPath(params.blacklist, checkIfExists: true).collect()              : Channel.empty()
chGtf                = params.gtf                   ? Channel.fromPath(params.gtf, checkIfExists: true).collect()                    : Channel.empty()
chMetadata           = params.metadata              ? Channel.fromPath(params.metadata, checkIfExists: true).collect()               : Channel.empty()
chDesignFile         = params.design                ? Channel.fromPath(params.design, checkIfExists: true).collect()                 : Channel.empty()

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
  'EffGenomeSize' : params.effGenomeSize ?: null,
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

// Load raw reads
if (!params.bam){
  chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)
}else{
  chRawReads = Channel.empty()
  chInputBam = NFTools.getIntermediatesData(params.samplePlan, '.bam', params)
}

// Make samplePlan if not available
chSplan = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

// Design
chDesignControl = params.design ? loadDesign(params.design) : Channel.empty()

/*
==================================
           INCLUDE
==================================
*/ 

// Workflows
include { prepareAnnotationFlow } from './nf-modules/local/subworkflow/prepareAnnotation'
include { mappingFlow } from './nf-modules/local/subworkflow/mapping'
include { loadBamFlow } from './nf-modules/local/subworkflow/loadBam'
include { bamFilteringFlow as bamFilteringFlowRef } from './nf-modules/local/subworkflow/bamFiltering'
include { bamFilteringFlow as bamFilteringFlowSpike } from './nf-modules/local/subworkflow/bamFiltering'
include { bamChipFlow }      from './nf-modules/local/subworkflow/bamChip'
include { bamSpikesFlow }    from './nf-modules/local/subworkflow/bamSpikes'
include { peakCallingFlow }  from './nf-modules/local/subworkflow/peakcalling'

// Processes
include { checkDesign } from './nf-modules/local/process/checkDesign'
include { fastqc }      from './nf-modules/common/process/fastqc/fastqc'
include { trimGalore }  from './nf-modules/common/process/trimGalore/trimGalore'
include { preseq }      from './nf-modules/common/process/preseq/preseq'
include { featureCounts }       from './nf-modules/common/process/featureCounts/featureCounts'
include { samtoolsIndex }       from './nf-modules/common/process/samtools/samtoolsIndex'
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { multiqc }             from './nf-modules/local/process/multiqc'


workflow {
  chVersions = Channel.empty()

  main:
  // Init MultiQC Channels
  chFastqcMqc = Channel.empty()
  chPreseqMqc = Channel.empty()
  chPeaksOutput=Channel.empty()
  chPeaksCountsMqc=Channel.empty()
  chFripResults=Channel.empty()
  chPeaksQCMqc=Channel.empty()
  chTrimmingMqc = Channel.empty()

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

  // PROCESS: trimming
  if (params.trimming){
    trimGalore(
      chRawReads
    )
    chRawReads=trimGalore.out.fastq
    chTrimmingMqc=trimGalore.out.logs.map{it->it[1]}
    chVersions = chVersions.mix(trimGalore.out.versions)
  }

  // PROCESS: fastqc
  fastqc(
    chRawReads
  )
  chFastqcMqc = fastqc.out.results.collect()
  chVersions = chVersions.mix(fastqc.out.versions)


  //*******************************************
  // MAPPING

  if (!params.bam){
    mappingFlow(
      chRawReads,
      chMappingIndex.collect(),
      chSpikeIndex,
    )
    chAlignedBam = mappingFlow.out.bam
    chAlignedBamMqc = mappingFlow.out.logs
    chAlignedFlagstat = mappingFlow.out.flagstat
    chAlignedSpikeBam = mappingFlow.out.spikeBam
    chCompareBamsMqc = mappingFlow.out.compareBamsMqc
    chVersions = chVersions.mix(mappingFlow.out.versions)

    // Remove low mapping rate for spikes
    chCompareBamsMqc.join(chAlignedSpikeBam)
      .filter { prefix, logs, bam, bai -> checkSpikeAlignmentPercent(prefix, logs, params.spikePercentFilter) }
      .map { prefix, logs, bam, bai -> [ prefix, bam, bai ] }
      .set { chPassedSpikeBam }
  }else{
 
    loadBamFlow(
      chInputBam
    )

    chVersions = chVersions.mix(loadBamFlow.out.versions)
    chAlignedBam = loadBamFlow.out.bam
    chAlignedFlagstat = loadBamFlow.out.flagstat
    chAlignedBamMqc = Channel.empty()
    chAlignedSpikeBam = Channel.empty()
    chCompareBamsMqc = Channel.empty()
    chPassedSpikeBam = Channel.empty()     
  }

  //*******************************************
  // CHIP WORFLOW

  if (!params.skipSaturation){ 
    preseq(
      chAlignedBam
    )
    chPreseqMqc = preseq.out.results.collect()
    chVersions = chVersions.mix(preseq.out.versions)
  }

  bamFilteringFlowRef(
    chAlignedBam
  )
  chVersions = chVersions.mix(bamFilteringFlowRef.out.versions)

  bamChipFlow(
    bamFilteringFlowRef.out.bam,
    chBlacklist,
    chGeneBed,
    chEffGenomeSize
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
      chBlacklist,
      chEffGenomeSize
    )
    chVersions = chVersions.mix(bamSpikesFlow.out.versions)
  }

  //*********************************************
  // DOWNSTREAM ANALYSIS (DESIGN IS MANDATORY !)
 
  if (params.design && !params.skipPeakCalling){

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
  } 

  //*********************************************
  // COUNTS
  if (!params.skipFeatCounts){

    //SUB-WORKFLOW : prepare annotation files
    prepareAnnotationFlow(
      chGeneBed.collect()
    )

    featureCounts(
      bamFilteringFlowRef.out.bam.combine(prepareAnnotationFlow.out.saf)
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
      .flatMap{ it -> it[0].id + ": Poor spike alignment rate. Sample ignored !"}
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
      chTrimmingMqc.collect().ifEmpty([]),
      chFastqcMqc.collect().ifEmpty([]),
      chAlignedBamMqc.collect().ifEmpty([]),
      chCompareBamsMqc.map{it->it[1]}.collect().ifEmpty([]),
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
