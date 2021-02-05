/****************
 * Spike-in
 */

// Merging, if necessary reference aligned reads and spike aligned reads
process compareRefSpike{
     tag "${sample}"
     label 'compbam'
     label 'minCpu'
     label 'medMem'
     publishDir "${params.outDir}/spike", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename }

     input:
     tuple val(sample), path('unsortedBamSpike'), path('unsortedBamRef') 

     output:
     tuple val(sample), path("*_clean.bam")                          , emit: refBams 
     tuple val(sampleSpike), path("*.mqc"), path("*_clean_spike.bam"), emit: spikeBams
     path('*.mqc')                                                   , emit: mappingSpikeMqc
     path('*.log')

     script:
     sampleSpike = sample + '_spike'
     """
     samtools sort $unsortedBamRef -n -@ ${task.cpus} -T ${sample}_ref -o ${sample}_ref_sorted.bam
     samtools sort $unsortedBamSpike -n -@ ${task.cpus} -T ${sample}_spike -o ${sample}_spike_sorted.bam
     compareAlignments.py -a ${sample}_ref_sorted.bam -b ${sample}_spike_sorted.bam -oa ${sample}_clean.bam -ob ${sample}_clean_spike.bam 2> ${sample}_compareAlignments.log
     """
}
