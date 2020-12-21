/* 
 * all ChIP analysis 
 */

/* 
 * include requires tasks 
 */
include { PPQT } from '../processes/PPQT' 
include { bigWig } from '../processes/bigWig'
include { deepToolsComputeMatrix } from '../processes/deepToolsComputeMatrix'
// include { deepToolsCorrelationQC } from '../processes/deepToolsCorrelationQC'
// include { deepToolsFingerprint } from '../processes/deepToolsFingerprint'

/***********************
 * Header and conf
 */

//PPQT headers
chPpqtCorHeader = file("$baseDir/assets/ppqt_cor_header.txt", checkIfExists: true)
chPpqtNSCHeader = file("$baseDir/assets/ppqt_nsc_header.txt", checkIfExists: true)
chPpqtRSCHeader = file("$baseDir/assets/ppqt_rsc_header.txt", checkIfExists: true)

workflow bamsChipFlow {
    // required inputs
    take:
     chBamsChip 
     chBlacklist
     chGeneBed
    // workflow implementation
    main:

      PPQT(chBamsChip, chPpqtCorHeader, chPpqtNSCHeader, chPpqtRSCHeader)
      bigWig(chBamsChip, chBlacklist.collect())
      deepToolsComputeMatrix(bigWig.out.bigWig, chGeneBed.collect()) 

    // emit:
      chPpqtOutMqc = PPQT.out.ppqtOutMqc
      chPpqtCsvMqc = PPQT.out.ppqtCsvMqc
      chPPQTVersion = PPQT.out.version
}

