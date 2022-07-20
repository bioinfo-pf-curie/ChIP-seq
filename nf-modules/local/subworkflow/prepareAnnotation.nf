/* 
 * Prepare annotation files
 */

include { extractTSS } from '../../local/process/extractTSS' 
include { bed2saf } from '../../local/process/bed2saf'

workflow prepareAnnotationFlow {

  take:
  bed

  main:
  chVersions = Channel.empty()
 
  extractTSS(
    bed
  )
  //chVersions = chVersions.mix(extractTSS.out.versions)
    
  bed2saf(
    bed.concat(extractTSS.out.tss)
  )
  //chVersions = chVersions.mix(bed2saf.out.versions)

  emit:
  tss = extractTSS.out.tss
  saf = bed2saf.out.saf
  versions = chVersions
}

