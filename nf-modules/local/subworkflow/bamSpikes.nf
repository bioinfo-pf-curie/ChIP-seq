/* 
 * Spikes-in analysis 
 */

/* 
 * include requires tasks 
 */
include { deepToolsMultiBamSummary } from '../process/deepToolsMultiBamSummary'
include { getSpikeScalingFactor } from '../process/getSpikeScalingFactor'
include { bigWigSpikeNorm } from '../process/bigWigSpikeNorm'

workflow bamSpikesFlow {

  // required inputs
  take:
   bamsChip // [prefix,bam,bai]
   bamsSpikes 
   blacklist

  // workflow implementation
  main:
    chVersions = Channel.empty()

    deepToolsMultiBamSummary(
       bamsSpikes.map{it[1]}.collect(),
       bamsSpikes.map{it[2]}.collect()
    )
    chVersions = chVersions.mix(deepToolsMultiBamSummary.out.versions)

    getSpikeScalingFactor(
      deepToolsMultiBamSummary.out.output
    )
    chVersions = chVersions.mix(getSpikeScalingFactor.out.versions)

    getSpikeScalingFactor.out.tabSF 
      .splitCsv(header:false, sep:',')
      .map { row -> [row[0], row[1]]}
      .set{chScaleFactor}
    
    bamsChip
      .combine(chScaleFactor)
      .filter{it[0] == it[3]}
      .map { it -> it[0,1,2,4]}
      .set{chBigWigScaleFactor}

    bigWigSpikeNorm(
       chBigWigScaleFactor,
       blacklist.collect().ifEmpty([])
    ) 
    chVersions = chVersions.mix(bigWigSpikeNorm.out.versions)

   emit:
     bigwig = bigWigSpikeNorm.out.bigwig
     versions = chVersions
}

