/*
 * FRiP - Fraction of reads in Peaks
 */

process frip{
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  tag("${meta.id}")

  input:
  tuple val(meta), path(bam), path(stats), path(peaks)
  path(fripScoreHeader)

  output:
  path("*tsv"), emit: output
  path("versions.txt"), emit: versions

  script:
  """
  echo "BEDtools"\$(intersectBed 2>&1 | grep "Version" | cut -f2 -d:) > versions.txt
  READS_IN_PEAKS=\$(intersectBed -a ${bam} -b ${peaks} -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $stats | awk -v a="\$READS_IN_PEAKS" '{printf "${meta.id}\\t%.2f\\n", a/\$1}' | cat $fripScoreHeader - > ${peaks.baseName}_FRiP.tsv
  """
}

