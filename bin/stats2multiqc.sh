#!/bin/bash

splan=$1
design=$2
aligner=$3
is_pe=$4

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_uniquely_mapped_reads,Percent_of_uniquely_mapped_reads,Number_of_multiple_mapped_reads,Percent_of_multiple_mapped_reads,Number_of_duplicates,Percent_of_duplicates,Normalized_strand_correlation,Relative_strand_correlation,Fraction_of_reads_in_peaks" > mqc.stats

for sample in $all_samples
do

#SAMPLE NAME
  sname=$(grep "$sample" $splan | awk -F, '{print $2}')

#ALIGNMENT
  if [ $aligner == "bowtie2" ]; then
    nb_reads=$(grep "reads;" mapping/${sample}_mapping.log | sed 's/ .*//')
    nb_uniq_reads=$(grep "exactly" mapping/${sample}_mapping.log | awk '{print $1}')
    perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mult_reads=$(grep ">1" mapping/${sample}_mammping.log | awk '{print $1}')
    perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  elif [ $aligner == "bwa-mem" ]; then
    nb_reads=$(grep 'Total' mapping/${sample}_mapping.log | awk -F "=" '{print $2}' | sed 's/ //g')
    nb_uniq_reads=$(grep 'Uniquely' mapping/${sample}_mapping.log | awk -F "=" '{print $2}'| sed 's/ //g')
    perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mult_reads=$(grep 'Multi' mapping/${sample}_mapping.log | awk -F "=" '{print $2}'| sed 's/ //g')
    perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  elif [ $aligner == "star" ]; then
    nb_reads=$(grep "Number of input reads" mapping/${sample}_mapping.log | cut -d"|" -f 2 | sed -e 's/\t//g')
    nb_uniq_reads=$(grep "Uniquely.*number" mapping/${sample}_mapping.log | awk '{print $NF}')
    perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mult_reads=$(grep "Number.*multiple" mapping/${sample}_mapping.log | awk '{print $NF}')
    perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  fi
  
#PICARD
  nb_dups=$(tail -3 picard/${sample}.MarkDuplicates.metrics.txt | sed '/^$/d' | awk '{print $6 + $7}')
  perc_dups=$(tail -3 picard/${sample}.MarkDuplicates.metrics.txt | sed '/^$/d' | awk '{print ($6 + $7)*100 / ($3 + $4)}')

#PPQT
  if [[ -e ppqt/${sample}.spp.out && -e ppqt/${sample}_spp_nsc_mqc.tsv ]]; then
      frag_length=$(awk '{print $3}' ppqt/${sample}.spp.out | sed 's/,.*//')
      nsc=$(grep "$sample" ppqt/${sample}_spp_nsc_mqc.tsv | awk '{print $2}')  
      rsc=$(grep "$sample" ppqt/${sample}_spp_rsc_mqc.tsv | awk '{print $2}')
  else
      frag_length='NA'
      nsc='NA'
      rsc='NA'
  fi

#PeakCalling 
  peaktype=$(grep "^$sample" ${design} | awk -F, '{print $5}') 
  if [ -e peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv ]; then
      frip=$(grep "$sample" peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv | awk '{print $2}')
  else
      frip='NA'
  fi

#To file
  echo -e ${sample},${sname},${nb_reads},${frag_length},${nb_mapped},${perc_mapped},${nb_uniq_reads},${perc_uniq_reads},${nb_mult_reads},${perc_mult_reads},${nb_dups},${perc_dups},${nsc},${rsc},${frip} >> mqc.stats
done

