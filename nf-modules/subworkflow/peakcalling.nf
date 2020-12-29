/***********************
 * Peak calling
 */

/* 
 * include requires tasks 
 */
include { sharpMACS2 } from '../processes/sharpMACS2' 
include { broadMACS2 } from '../processes/broadMACS2'
// include { veryBroadEpic2 } from '../processes/veryBroadEpic2'
// include { peakAnnoHomer } from '../processes/peakAnnoHomer'
// include { peakQC } from '../processes/peakQC'
// include { IDR } from '../processes/IDR'

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .set { chPeakCountHeader }

Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .set { chFripScoreHeader }

Channel
  .fromPath("$baseDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }

workflow peakCallingFlow {
    // required inputs
    take:
     chBamsChip 
     chDesignControl
     chNoInput
     chFlagstatMacs 
     
    // workflow implementation
    main:

      /*
       * Prepare channels
       */
      if (params.design){
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
      
     //  chGroupBamMacs
      //  .filter { it[2] == 'very-broad' }
      //  .dump(tag:'peakCall')
      //  .set { chGroupBamMacsVeryBroad }
    }
   // else{
    //  chGroupBamMacsSharp=Channel.empty()
    //  chGroupBamMacsBroad=Channel.empty()
    //  chGroupBamMacsVeryBroad=Channel.empty()
    //}
 
      /*
       * MACS2 - sharp mode
       */
      
      // chGroupBamMacs.filter { it[2] == 'sharp' }.view()
      sharpMACS2(
	chGroupBamMacs.filter { it[2] == 'sharp' },
	chPeakCountHeader.collect(),
	chFripScoreHeader.collect()
      )
      
      /*
       * MACS2  - Broad
       */
      broadMACS2(
        chGroupBamMacs.filter { it[2] == 'broad' },
        chPeakCountHeader.collect(),
        chFripScoreHeader.collect()
      )


    emit:

      version = broadMACS2.out.version // channel: [ path("v_macs2.txt") ]
}

