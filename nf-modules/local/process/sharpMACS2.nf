/*
 * MACS2 - sharp mode
 */

process sharpMACS2{
  tag "${sampleID} - ${controlID}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/sharp", mode: 'copy',
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
  path("*.xls"), emit: macsOutputSharp
  tuple val(group), val(replicate), val(peaktype), val(sampleID), path("*.narrowPeak"), emit:  peaksMacsSharp
  path "*_mqc.tsv", emit: macsCountsSharp
  path("v_macs2.txt"), emit: version

  script:
  format = params.singleEnd ? "BAM" : "BAMPE"
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  """
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  macs2 callpeak \\
    -t ${sampleBam[0]} \\
    ${ctrl} \\
    -f $format \\
    -g $params.effGenomeSize \\
    -n $sampleID \\
    --SPMR --trackline --bdg \\
    --keep-dup all
  cat ${sampleID}_peaks.narrowPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.narrowPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
}

