/* 
 * all Spikes analysis 
 */

/* 
 * include requires tasks 
 */
include { getSpikeCountPerBin } from '../process/getSpikeCountPerBin' 
include { getSpikeScalingFactor } from '../process/getSpikeScalingFactor'
include { bigWigSpikeNorm } from '../process/bigWigSpikeNorm'

workflow bamsSpikesFlow {
    // required inputs
    take:
     chBamsSpikes 
     chBamsChip
     chBlacklist
    // workflow implementation
    main:

      getSpikeCountPerBin(
         chBamsSpikes.map{it[1][0]}.collect(),
         chBamsSpikes.map{it[1][1]}.collect()
      )

      getSpikeScalingFactor(getSpikeCountPerBin.out)

      getSpikeScalingFactor.out.tabSF 
        .splitCsv(header:false, sep:',')
        .map { row -> [row[0], row[1]]}
        .set{chScaleFactor}

      
      chBamsChip
        .combine(chScaleFactor)
        .filter{it[0] == it[2]}
        .map { it -> it[0,1,3]}
        .set{chBigWigScaleFactor}

      bigWigSpikeNorm(
         chBigWigScaleFactor,
         chBlacklist.collect().ifEmpty([])
      ) 

     emit:
     chBigWigSF = bigWigSpikeNorm.out.bigWigSF // channel: [ tuple val(prefix), path('*.bigwig') ]
}

