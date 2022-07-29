/* 
 * Spikes-in analysis 
 */

include { deeptoolsMultiBamSummary } from '../../common/process/deeptools/deeptoolsMultiBamSummary'
include { getSpikeScalingFactor } from '../process/getSpikeScalingFactor'
include { deeptoolsBamCoverage } from '../../common/process/deeptools/deeptoolsBamCoverage'

workflow bamSpikesFlow {

  // required inputs
  take:
   bamsChip // [meta,bam,bai]
   bamsSpikes 
   blacklist
   effGenomeSize

  // workflow implementation
  main:
    chVersions = Channel.empty()

    deeptoolsMultiBamSummary(
       bamsSpikes.map{it[1]}.collect(),
       bamsSpikes.map{it[2]}.collect()
    )
    chVersions = chVersions.mix(deeptoolsMultiBamSummary.out.versions)

    getSpikeScalingFactor(
      deeptoolsMultiBamSummary.out.output
    )
    chVersions = chVersions.mix(getSpikeScalingFactor.out.versions)

    getSpikeScalingFactor.out.tabSF 
      .splitCsv(header:false, sep:',')
      .map { row -> [row[0], row[1]]}
      .set{chScaleFactor}
    
    bamsChip
      .combine(chScaleFactor)
      .filter{it[0].id == it[3]}
      .map { it -> it[0,1,2,4]}
      .set{chBigWigScaleFactor}

    deeptoolsBamCoverage(
       chBigWigScaleFactor,
       blacklist,
       effGenomeSize
    ) 
    chVersions = chVersions.mix(deeptoolsBamCoverage.out.versions)

   emit:
     bigwig = deeptoolsBamCoverage.out.bigwig
     versions = chVersions
}

