/*
 * MACS2  - Broad
 */

process broadMACS2{
  tag "${sampleID} - ${controlID}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/broad", mode: 'copy',
    saveAs: { filename ->
            if (filename.endsWith(".tsv")) "stats/$filename"
            else filename
            }

  when:
  !params.skipPeakCalling && params.effGenomeSize

  input:
  tuple val(group), val(replicate), val(peaktype), val(sampleID), path(sampleBam), path(sampleFlagstat), val(controlID), path(controlBam)
  path peakCountHeader
  path fripScoreHeader

  output:
  path("*.xls"), emit:macsOutputBroad
  tuple val(group), val(replicate), val(peaktype), val(sampleID), path("*.broadPeak"), emit: peaksMacsBroad
  path "*_mqc.tsv", emit: macsCountsBroad
  path("v_macs2.txt"), emit: version 

  script:
  broad = "--broad --broad-cutoff ${params.broadCutoff}"
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  """
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    ${broad} \\
    -f $format \\
    -g $params.effGenomeSize \\
    -n $sampleID \\
    --SPMR --trackline --bdg \\
    --keep-dup all
  cat ${sampleID}_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.broadPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
}

