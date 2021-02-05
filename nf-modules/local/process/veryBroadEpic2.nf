/*
 * EPIC2 - very broad
 */

process veryBroadEpic2{
  tag "${sampleID} - ${controlID}"
  label 'epic2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/very-broad", mode: 'copy',
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
  path chromsize

  output:
  tuple val(group), val(replicate), val(peaktype), val(sampleID), path("*.broadPeak"), emit: peaksEpic
  path ("*.out"), emit: epicOutput
  path "*_mqc.tsv", emit: macsCountsVbroad
  path("v_epic2.txt"), emit: version

  script:
  ctrl = controlID != 'NO_INPUT' ? "-c ${controlBam[0]}" : ''
  opts = (params.singleEnd && params.fragmentSize > 0) ? "--fragment-size ${params.fragmentSize}" : ""
  """
  epic2 --version &> v_epic2.txt
  epic2 -t ${sampleBam[0]} \\
    ${ctrl} \\
    --chromsizes ${chromsize} \\
    --effective-genome-fraction ${params.effGenomeSize} \\
    -a \\
    --bin-size 200 \\
    --gaps-allowed 3 \\
    --false-discovery-rate-cutoff 0.05 \\
    -o ${sampleID}_epic.out

  echo "track type=broadPeak name=\"${sampleID}\" description=\"${sampleID}\" nextItemButton=on" > ${sampleID}_peaks.broadPeak
  awk -v id="${sampleID}" 'NF<10{fc=0;pval=-1;qval=-1}NF==10{fc=\$10;pval=\$4;qval=\$9}NR>1{OFS="\t"; print \$1,\$2,\$3,id"_peak_"NR,\$5,".",fc,pval,qval}' ${sampleID}_epic.out >> ${sampleID}_peaks.broadPeak
  cat ${sampleID}_peaks.broadPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${sampleID}", \$1 }' | cat $peakCountHeader - > ${sampleID}_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${sampleBam[0]} -b ${sampleID}_peaks.broadPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sampleID}", a/\$1}' | cat $fripScoreHeader - > ${sampleID}_peaks.FRiP_mqc.tsv
  """
}

