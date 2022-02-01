/*
 * Compare reference and spike-in alignment
 */

// Merging, if necessary reference aligned reads and spike aligned reads
process compareBams{
     tag "${prefix}"
     label 'compbam'
     label 'minCpu'
     label 'medMem'
     publishDir "${params.outDir}/spike", mode: 'copy',
              saveAs: {filename ->
              if (filename.indexOf(".log") > 0) "logs/$filename"
              else filename }

     input:
     tuple val(prefix), path('unsortedBamRef'), path('unsortedBamSpike')
     val(genome)
     val(spike)

     output:
     tuple val(prefix), path("*${refSuffix}"), emit: refBam
     tuple val(prefix), path("*${spikeSuffix}"), emit: spikeBam
     path("*mqc"), emit: mqc

     script:
     refSuffix = "clean_${genome}.bam"
     spikeSuffix = "clean_${spike}.bam"
     """
     samtools sort $unsortedBamRef -n -@ ${task.cpus} -T ${prefix}_ref -o ${prefix}_ref_sorted.bam
     samtools sort $unsortedBamSpike -n -@ ${task.cpus} -T ${prefix}_spike -o ${prefix}_spike_sorted.bam
     compareAlignments.py -a ${prefix}_ref_sorted.bam -b ${prefix}_spike_sorted.bam -oa ${prefix}_${refSuffix} -ob ${prefix}_${spikeSuffix} 2> ${prefix}_compareAlignments.log
     """
}
