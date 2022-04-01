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
chPpqtCorHeader = file("$projectDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$projectDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$projectDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

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
      chFragmentSize = getFragmentSize.out.results
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
    ppqtOutMqc              = !params.skipPPQT ? PPQT.out.ppqtOutMqc : Channel.empty()
    ppqtCsvMqc              = !params.skipPPQT ? PPQT.out.ppqtCsvMqc : Channel.empty()
    deeptoolsProfileMqc     = !params.skipDeepTools ? deepToolsComputeMatrix.out.mqc : Channel.empty()
    deeptoolsCorrelateMqc   = !params.skipDeepTools ? deepToolsCorrelationQC.out.mqc : Channel.empty()
    deeptoolsFingerprintMqc = !params.skipDeepTools ? deepToolsFingerprint.out.mqc : Channel.empty()
    versions                = chVersions
}

