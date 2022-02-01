/* 
 * all ChIP analysis 
 */

include { getFragmentSize } from '../process/getFragmentSize' 
include { PPQT } from '../process/PPQT' 
include { bigWig } from '../process/bigWig'
include { deepToolsComputeMatrix } from '../process/deepToolsComputeMatrix'
include { deepToolsCorrelationQC } from '../process/deepToolsCorrelationQC'
include { deepToolsFingerprint } from '../process/deepToolsFingerprint'

/***********************
 * Header and conf
 */

// PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

workflow bamChipFlow {

  take:
    bam // [prefix, bam, bai]
    blacklist
    geneBed

  main:
    chVersions = Channel.empty()

    // Fragment size
    if (!params.singleEnd){
      getFragmentSize(
        bam
      )
      chVersions = chVersions.mix(getFragmentSize.out.versions)
      chFragmentSize = getFragmentSize
    }else{
      chFragmentSize = Channel.empty()
    }

    // Phantom Peak
    if (!params.skipPPQT){
      PPQT(
        bam,
        chPpqtCorHeader,
        chPpqtNSCHeader,
        chPpqtRSCHeader
      )
      chVersions = chVersions.mix(PPQT.out.versions)
      chPpqtOut = PPQT.out.ppqtCsvMqc
      chPpqtCsv = PPQT.out.ppqtOutMqc
    }else{
      chPpqtOut = Channel.empty()
      chPpqtCsv = Channel.empty()
    }

    bigWig(
      bam,
      blacklist.collect().ifEmpty([])
    )

    if (!params.skipDeepTools){
      deepToolsComputeMatrix(
        bigWig.out.bigWig,
        geneBed.collect()
      ) 
      chVersions = chVersions.mix(deepToolsComputeMatrix.out.versions) 

      deepToolsCorrelationQC(
        bam.map{it[0]}.collect(),
        bam.map{it[1]}.collect(), 
        bam.map{it[2]}.collect(),
        blacklist.ifEmpty([])
      )
      chVersions = chVersions.mix(deepToolsCorrelationQC.out.versions) 

      deepToolsFingerprint(
        bam.map{it[0]}.collect(),
        bam.map{it[1]}.collect(),
        bam.map{it[2]}.collect()
      )
      chVersions = chVersions.mix(deepToolsFingerprint.out.versions) 
    }

  emit:
    fragmentsSize           = chFragmentSize
    ppqtOutMqc              = chPpqtCsv
    ppqtCsvMqc              = chPpqtOut
    deeptoolsProfileMqc     = deepToolsComputeMatrix.out.mqc
    deeptoolsCorrelateMqc   = deepToolsCorrelationQC.out.mqc
    deeptoolsFingerprintMqc = deepToolsFingerprint.out.mqc
    versions = chVersions
}

