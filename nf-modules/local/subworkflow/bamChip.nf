/* 
 * all ChIP analysis 
 */

include { collectInsertSizeMetrics } from '../../common/process/picard/collectInsertSizeMetrics' 

include { PPQT } from '../../local/process/PPQT'
include { getDepthNormFactor } from '../../local/process/getDepthNormFactor'

include { deeptoolsBamCoverage } from '../../common/process/deeptools/deeptoolsBamCoverage'
include { deeptoolsComputeMatrix } from '../../common/process/deeptools/deeptoolsComputeMatrix'
include { deeptoolsCorrelationQC } from '../../common/process/deeptools/deeptoolsCorrelationQC'
include { deeptoolsFingerprint } from '../../common/process/deeptools/deeptoolsFingerprint'

/***********************
 * Header and conf
 */

// PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

workflow bamChipFlow {

  take:
    bam // [meta, bam, bai]
    blacklist
    geneBed
    effGenomeSize

  main:
    chVersions = Channel.empty()

    // Fragment size
    if (!params.singleEnd){
      collectInsertSizeMetrics(
        bam
      )
      chVersions = chVersions.mix(collectInsertSizeMetrics.out.versions)
      chFragmentSize = collectInsertSizeMetrics.out.results
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

    getDepthNormFactor(
      bam
    )

    deeptoolsBamCoverage(
      bam.join(getDepthNormFactor.out.sf),
      blacklist.collect().ifEmpty([]),
      effGenomeSize.collect().ifEmpty([])
    )

    if (!params.skipDeepTools){
      deeptoolsComputeMatrix(
        deeptoolsBamCoverage.out.bigwig,
        geneBed.collect()
      ) 
      chVersions = chVersions.mix(deeptoolsComputeMatrix.out.versions) 

      deeptoolsCorrelationQC(
        bam.map{it[0].id}.collect(),
        bam.map{it[1]}.collect(), 
        bam.map{it[2]}.collect(),
        blacklist.ifEmpty([])
      )
      chVersions = chVersions.mix(deeptoolsCorrelationQC.out.versions) 

      deeptoolsFingerprint(
        bam.map{it[0].id}.collect(),
        bam.map{it[1]}.collect(),
        bam.map{it[2]}.collect()
      )
      chVersions = chVersions.mix(deeptoolsFingerprint.out.versions)
    }

  emit:
    fragmentsSize           = chFragmentSize
    ppqtOutMqc              = !params.skipPPQT ? PPQT.out.ppqtOutMqc : Channel.empty()
    ppqtCsvMqc              = !params.skipPPQT ? PPQT.out.ppqtCsvMqc : Channel.empty()
    deeptoolsProfileMqc     = !params.skipDeepTools ? deeptoolsComputeMatrix.out.mqc : Channel.empty()
    deeptoolsCorrelateMqc   = !params.skipDeepTools ? deeptoolsCorrelationQC.out.mqc : Channel.empty()
    deeptoolsFingerprintMqc = !params.skipDeepTools ? deeptoolsFingerprint.out.mqc : Channel.empty()
    versions                = chVersions
}

