/***********************
 * Peak calling
 */

/* 
 * include requires tasks 
 */
include { sharpMACS2 } from '../processes/sharpMACS2' 
// include { broadMACS2 } from '../processes/broadMACS2'
// include { veryBroadEpic2 } from '../processes/veryBroadEpic2'
// include { peakAnnoHomer } from '../processes/peakAnnoHomer'
// include { peakQC } from '../processes/peakQC'
// include { IDR } from '../processes/IDR'

workflow peakCallingFlow {
    // required inputs
    take:
     chBamsChip 
     chDesignControl
     chNoInput
     chFlagstatMacs 
     chPeakCountHeader
     chFripScoreHeader
     chPeakAnnotationHeader
     
    // workflow implementation
    main:
      /*
      * Prepare channels
      */
      chBamsChip 
       .join(chFlagstatMacs)
       .combine(chNoInput.concat(chBamsChip))
       .set { chBamsChip }  

       chDesignControl
       .combine(chBamsChip)
       .filter { it[0] == it[5] && it[1] == it[8] }
       .map { it ->  it[2..-1] }
       .dump(tag:'peakCall')
       .set { chGroupBamMacs }  
 

      sharpMACS2(chGroupBamMacs.filter { it[2] == 'sharp' }, chPeakCountHeader.collect(), chFripScoreHeader.collect())

    // emit:
}

