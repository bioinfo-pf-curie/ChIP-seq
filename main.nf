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

nextflow.enable.dsl=2

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
  --singleEnd [bool]                 Specifies that the input is single end reads. Default: false
  --fragmentSize [int]               Estimated fragment length used to extend single-end reads. Default: 200
  --spike [str]                      Name of the genome used for spike-in analysis. Default: false
  --genomeAnnotationPath [file]      Path to genome annotation folder

  Annotation: If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --fasta [file]                     Path to Fasta reference
  --spikeFasta [file]                Path to Fasta reference for spike-in
  --geneBed [file]                   BED annotation file with gene coordinate.
  --gtf [file]                       GTF annotation file. Used in HOMER peak annotation
  --effGenomeSize [int]              Effective Genome size

  Alignment: If you want to modify default options or wish to overwrite any of the indexes given by the --genome field
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
  --blacklist [file]                 Path to black list regions (.bed). See the genome.config for details.
  --spikePercentFilter [float]       Minimum percent of reads aligned to spike-in genome. Default: 0.2

  Analysis:
  --noReadExtension [bool]           Do not extend reads to fragment length. Default: false
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
def genomeRef = params.genome

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

params.geneBed = genomeRef ? params.genomes[ genomeRef ].geneBed ?: false : false
if (params.geneBed) {
  Channel
    .fromPath(params.geneBed, checkIfExists: true)
    .ifEmpty {exit 1, "BED file ${geneBed} not found"}
    .set{chGeneBed}
}else{
  chGeneBed = Channel.empty()
}

params.blacklist = genomeRef ? params.genomes[ genomeRef ].blacklist ?: false : false
if (params.blacklist) { 
  Channel
    .fromPath(params.blacklist, checkIfExists: true)
    .set {chBlacklist} 
}else{
  chBlacklist = Channel.empty()
}


/***********************
 * Header and conf
 */

//Peak Calling
params.effGenomeSize = genomeRef ? params.genomes[ genomeRef ].effGenomeSize ?: false : false
if (!params.effGenomeSize) {
  log.warn "=================================================================\n" +
            "  WARNING! Effective Genome Size is not defined.\n" +
            "  Peak calling and annotation will be skipped.\n" +
            "  Please specify value for '--effGenomeSize' to run these steps.\n" +
            "================================================================"
}

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
      .set { rawReads }
  }else if (!params.singleEnd && !params.inputBam){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .set { rawReads }
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
      .set { rawReads }
   } else {
     Channel
       .from(params.readPaths)
       .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
       .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
       .set { rawReads }
   }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { rawReads }
}


/**************************
 * Make sample plan if not available
 */

if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .set {chSplan}
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .set{ chSplan}
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan }
  }
} else if(params.bamPaths){
  Channel
    .from(params.bamPaths)
    .collectFile() {
      item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
     }
    .set{ chSplan }
  params.aligner = false
} else {
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
      }     
      .set { chSplan}
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }     
      .set { chSplan }
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
    .set { chDesignCheck }

  chDesignControl = chDesignCheck

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
}

// Workflows
// QC : check design and factqc
include { qcFlow } from './nf-modules/subworkflow/qc'
// Alignment on reference genome
include { mappingFlow } from './nf-modules/subworkflow/mapping' 
// Spike-in and Sorting BAM files
include { sortingFlow } from './nf-modules/subworkflow/sorting' 
include { markdupFlow } from './nf-modules/subworkflow/markdup' 
include { bamsChipFlow } from './nf-modules/subworkflow/bamschip' 
include { bamsSpikesFlow } from './nf-modules/subworkflow/bamsspikes' 
// Peak calling
include { peakCallingFlow } from './nf-modules/subworkflow/peakcalling' 

// Processes
include { prepareAnnotation } from './nf-modules/processes/prepareAnnotation'
include { featureCounts } from './nf-modules/processes/featureCounts'
include { getSoftwareVersions } from './nf-modules/processes/getSoftwareVersions'
include { workflowSummaryMqc } from './nf-modules/processes/workflowSummaryMqc'
include { multiqc } from './nf-modules/processes/multiqc'
include { outputDocumentation } from './nf-modules/processes/outputDocumentation'

workflow {
    main:

      chTSSFeatCounts = prepareAnnotation(chGeneBed.collect())

      // subroutines
      outputDocumentation(
        chOutputDocs,
        chOutputDocsImages
      )

      // QC : check design and factqc
      qcFlow(
        chDesignCheck,
        chSplan,
        rawReads
      )
      chFastqcVersion = qcFlow.out.version
      chFastqcVersion.view()    
      chFastqcMqc = qcFlow.out.mqc

      // Alignment on reference genome
      mappingFlow(
	rawReads,
	chBwaIndex,
	chBt2Index,
	chStarIndex
      )
      chAlignReads = mappingFlow.out.bam
      chMappingMqc = mappingFlow.out.mqc
      chMappingMqc.view()

      if (params.inputBam){
        chFastqcMqc = Channel.empty()
        chMappingMqc = Channel.empty()
      }

      // Spike-in and Sorting BAM files
      sortingFlow(chAlignReads, useSpike)
      chStatsMqc = sortingFlow.out.statsMqc

      markdupFlow(
	sortingFlow.out.sortBams
      )
      chBams = markdupFlow.out.chFilteredBams
      chFlagstat = markdupFlow.out.chFilteredFlagstat
      //chBams.view()
      // Separate sample BAMs and spike BAMs
      chFlagstatMacs = Channel.empty()
      chFlagstatSpikes = Channel.empty()
      chFlagstatMacs = chFlagstat
      chFlagstat 
        .branch { prefix: it[0] =~ 'spike'}
        .set { chFlagstatSpikes }

      chBamsChip = Channel.empty()
      chBamsSpikes = Channel.empty() 
      chBamsChip = chBams
      chBams
       .branch { prefix: it[0] =~ 'spike'}
       .set { chBamsSpikes} 

      chBamsSpikes.view()
      chBamsChip.view()

      // Preparing all ChIP data for further analysis
      chBamsChip = chBamsChip.dump(tag:'cbams')

      // all ChIP analysis
      bamsChipFlow(
	chBamsChip,
	chBlacklist,
	chGeneBed
      )

      if (useSpike){
        // all Spikes analysis
        bamsSpikesFlow(
          chBamsSpikes
          chBamsChip,
          chBlacklist
        )
      }

      // /!\ From this point, 'design' is mandatory /!\

      // Peak calling
      peakCallingFlow(
        chBamsChip,
        chDesignControl,
        chNoInput,
        chFlagstatMacs
      )
 
      // Feature counts
      featureCounts(
        chBamsChip.map{items->items[1][0]}.collect(),
        chGeneBed.concat(chTSSFeatCounts)
      ) 
      
      // MultiQC
      getSoftwareVersions(
        qcFlow.out.version.first().ifEmpty([]),
        mappingFlow.out.chBwaVersion.first().ifEmpty([]),
        mappingFlow.out.chBowtie2Version.first().ifEmpty([]),
        mappingFlow.out.chStarVersion.first().ifEmpty([]),
        sortingFlow.out.chSamtoolsVersionBamSort.concat(markdupFlow.out.chSamtoolsVersionBamFiltering).first().ifEmpty([]),
        markdupFlow.out.chPicardVersion.first().ifEmpty([]),
        peakCallingFlow.out.chMacs2VersionMacs2Sharp.concat(peakCallingFlow.out.chMacs2VersionMacs2Broad).first().ifEmpty([]),
        peakCallingFlow.out.chEpic2Version.first().ifEmpty([]),
        markdupFlow.out.chPreseqVersion.first().ifEmpty([]),
        peakCallingFlow.out.chIdrVersion.first().ifEmpty([]),
        bamsChipFlow.out.chPPQTVersion.first().ifEmpty([]),
        bamsChipFlow.out.chDeeptoolsVersion.first().ifEmpty([]),
        featureCounts.out.version.first().ifEmpty([])
      ) 
      //workflowSummaryMqc( )
      //multiqc( )

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
    //println endWfSummary
    view(endWfSummary)
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

