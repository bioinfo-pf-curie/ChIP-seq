#!/bin/bash

splan=$1
design=$2
aligner=$3
is_pe=$4

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_uniquely_mapped_reads,Percentage_of_uniquely_mapped_reads,Number_of_multiple_mapped_reads,Percentage_of_multiple_mapped_reads,Number_of_duplicates,Percentage_of_duplicates,Normalized_strand_correlation,Relative_strand_correlation,Fraction_of_reads_in_peaks" > mqc.stats

for sample in $all_samples
do

#SAMPLE NAME
  sname=$(grep "$sample" $splan | awk -F, '{print $2}')

#ALIGNMENT
  if [ $aligner == "bowtie2" ]; then
    nb_reads=$(grep "reads;" alignment/reference/$sample.log | sed 's/ .*//')
    nb_uniq_reads=$(grep "exactly" alignment/reference/$sample.log | awk '{print $1}')
    perc_uniq_reads=$(grep "exactly" alignment/reference/$sample.log | awk '{print substr($2, 2, length($2) - 3)}')
    nb_mult_reads=$(grep ">1" alignment/reference/$sample.log | awk '{print $1}')
    perc_mult_reads=$(grep ">1" alignment/reference/$sample.log | awk '{print substr($2, 2, length($2) - 3)}')
  elif [ $aligner == "bwa-mem" ]; then
    nb_reads=$(grep 'Total' alignment/reference/$sample.log | awk '{print $4}')
    nb_uniq_reads=$(grep 'Uniquely' alignment/reference/$sample.log | awk '{print $5}')
    perc_uniq_reads=$(echo "scale=2 ; 100 * $nb_uniq_reads / $nb_reads" | bc)
    nb_mult_reads=$(grep 'Multiply' alignment/reference/$sample.log | awk '{print $5}')
    perc_mult_reads=$(echo "scale=2 ; 100 * $nb_mult_reads / $nb_reads" | bc)
  elif [ $aligner == "star" ]; then
    nb_reads=$(grep "Uniquely.*%" alignment/reference/$sample.log | awk '{print substr($NF, 1, length($NF)-1)}')  
    nb_uniq_reads=$(grep "Uniquely.*number" alignment/reference/$sample.log | awk '{print $NF}')
    perc_uniq_reads=$(grep "Uniquely.*number" alignment/reference/$sample.log | awk '{print $NF}')
    nb_mult_reads=$(grep "Number.*multiple" alignment/reference/$sample.log | awk '{print $NF}')
    perc_mult_reads=$(grep "%.*multiple" alignment/reference/$sample.log | awk '{print substr($NF, 1, length($NF)-1)}')
  fi
  
#PICARD
  nb_dups=$(tail -3 picard/${sample}.MarkDuplicates.metrics.txt | sed '/^$/d' | awk '{print $6 + $7}')
  perc_dups=$(tail -3 picard/${sample}.MarkDuplicates.metrics.txt | sed '/^$/d' | awk '{print ($6 + $7)*100 / ($3 + $4)}')

#PPQT
  frag_length=$(awk '{print $3}' ppqt/${sample}.spp.out | sed 's/,.*//')
  nsc=$(grep "$sample" ppqt/${sample}_spp_nsc_mqc.tsv | awk '{print $2}')  
  rsc=$(grep "$sample" ppqt/${sample}_spp_rsc_mqc.tsv | awk '{print $2}')

#PeakCalling 
  peaktype=$(grep "^$sample" ${design} | awk -F, '{print $5}') 
  if [ -z "$peaktype" ]; then
    frip='NA'
  else
    frip=$(grep "$sample" peakCalling/${peaktype}/stats/${sample}_peaks.FRiP_mqc.tsv | awk '{print $2}')
  fi

#To file
  echo -e ${sample},${sname},${nb_reads},${frag_length},${nb_uniq_reads},${perc_uniq_reads},${nb_mult_reads},${perc_mult_reads},${nb_dups},${perc_dups},${nsc},${rsc},${frip} >> mqc.stats
done

