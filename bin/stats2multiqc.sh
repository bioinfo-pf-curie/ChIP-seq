#!/bin/bash

function usage {
    echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -d DESIGN -a ALIGNER [-p][-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "stat2multiqc.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -s SAMPLE_PLAN"
    echo "   -d DESIGN"
    echo "   -a ALIGNER"
    echo "   [-p]: paired-end mode"
    echo "   [-h]: help"
    exit;
}

while getopts "s:d:a:ph" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
	d) design=$OPTARG;;
	a) aligner=$OPTARG;;
	p) is_pe=1;;
	h) help ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    usage
	    exit 1
	    ;;
    esac
done

if  [[ -z $splan ]]; then
    usage
    exit
fi

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_uniquely_mapped_reads,Percent_of_uniquely_mapped_reads,Number_of_multiple_mapped_reads,Percent_of_multiple_mapped_reads,Number_of_duplicates,Percent_of_duplicates,Normalized_strand_correlation,Relative_strand_correlation,Fraction_of_reads_in_peaks" > mqc.stats

for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample" $splan | awk -F, '{print $2}')

    #ALIGNMENT
    if [ $aligner == "bowtie2" ]; then
	nb_reads=$(grep "reads;" mapping/${sample}_bowtie2.log | sed 's/ .*//')
	if [[ $is_pe == "1" ]]; then
           nb_uniq_reads=$(grep "concordantly exactly" mapping/${sample}_bowtie2.log | awk '{print $1}')
           perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
           nb_mult_reads=$(grep "concordantly >1" mapping/${sample}_bowtie2.log | awk '{print $1}')
           perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
           nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
           perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	else 
	    nb_uniq_reads=$(grep "exactly" mapping/${sample}_bowtie2.log | awk '{print $1}')
	    perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	    nb_mult_reads=$(grep ">1" mapping/${sample}_bowtie2.log | awk '{print $1}')
	    perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	    nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
	    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	fi
    elif [ $aligner == "bwa-mem" ]; then
	nb_unmapped=$(grep 'Unmapped' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	nb_uniq_reads=$(grep 'Uniquely' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	nb_mult_reads=$(grep 'Multi' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	nb_reads=$(( $nb_unmapped + $nb_uniq_reads + $nb_mult_reads ))
	perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    elif [ $aligner == "star" ]; then
	nb_reads=$(grep "Number of input reads" mapping/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	nb_uniq_reads=$(grep "Uniquely.*number" mapping/${sample}Log.final.out | awk '{print $NF}')
	perc_uniq_reads=$(echo "${nb_uniq_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	nb_mult_reads=$(grep "Number.*multiple" mapping/${sample}Log.final.out | awk '{print $NF}')
	perc_mult_reads=$(echo "${nb_mult_reads} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	nb_mapped=$(($nb_uniq_reads + $nb_mult_reads))
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_reads='NA'
	nb_uniq_reads='NA'
	perc_uniq_reads='NA'
	nb_mult_reads='NA'
	perc_mult_reads='NA'
	nb_mapped='NA'
	perc_mapped='NA'
    fi
  
    #PICARD
    if [[ -e mapping/${sample}.MarkDuplicates.metrics.txt ]]; then
	if [[ $is_pe == "1" ]]; then
	    nb_dups=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	else
	    nb_dups=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	fi
	perc_dups=$(echo "${nb_dups} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_dups='NA'
	perc_dups='NA'
    fi

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
    if [[ ! -z $design ]];
    then
	peaktype=$(grep "^$sample" ${design} | awk -F, '{print $5}') 
	if [ -e peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv ]; then
	    frip=$(grep "$sample" peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv | awk '{print $2}')
	else
	    frip='NA'
	fi
    else
	frip='NA'
    fi

    #To file
    echo -e ${sample},${sname},${nb_reads},${frag_length},${nb_mapped},${perc_mapped},${nb_uniq_reads},${perc_uniq_reads},${nb_mult_reads},${perc_mult_reads},${nb_dups},${perc_dups},${nsc},${rsc},${frip} >> mqc.stats
done

