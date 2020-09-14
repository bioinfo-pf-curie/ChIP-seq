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

is_pe=0
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

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_hq_mapped_reads,Percent_of_hq_mapped_reads,Number_of_lq_mapped_reads,Percent_of_lq_mapped_reads,Number_of_duplicates,Percent_of_duplicates,Normalized_strand_correlation,Relative_strand_correlation,Fraction_of_reads_in_peaks" > mqc.stats

for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample," $splan | awk -F, '{print $2}')

    #ALIGNMENT
    if [ $aligner == "bowtie2" ]; then
	nb_reads=$(grep "reads;" mapping/${sample}_bowtie2.log | sed 's/ .*//')
    elif [ $aligner == "bwa-mem" ]; then
	nb_reads=$(grep 'Total' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	tail -n +3 mapping/${sample}_bwa.log > mapping/${sample}_bwa.mqc
    elif [ $aligner == "star" ]; then
	nb_reads=$(grep "Number of input reads" mapping/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
    fi

    #Mapping stats (always in reads - so must be converted for PE)
    #These statistics are calculated after spike cleaning but before filtering
    nb_mapped=$(awk -F, '$1=="Mapped"{print $2}' mapping/${sample}_mappingstats.mqc)
    nb_mapped_hq=$(awk -F, '$1=="HighQual"{print $2}' mapping/${sample}_mappingstats.mqc)
    nb_mapped_lq=$(awk -F, '$1=="LowQual"{print $2}' mapping/${sample}_mappingstats.mqc)

    if [[ $is_pe == 1 ]]; then
	nb_mapped=$(( $nb_mapped / 2))
	nb_mapped_hq=$(( $nb_mapped_hq / 2))
	nb_mapped_lq=$(( $nb_mapped_lq / 2))
    fi
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_hq=$(echo "${nb_mapped_hq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_lq=$(echo "${nb_mapped_lq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    #PICARD
    if [[ -e mapping/${sample}.MarkDuplicates.metrics.txt ]]; then
	nb_dups_pair=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	nb_dups_single=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	nb_dups_optical=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	nb_dups=$(($nb_dups_pair + $nb_dups_single + $nb_dups_optical))
	perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
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
    echo -e ${sample},${sname},${nb_reads},${frag_length},${nb_mapped},${perc_mapped},${nb_mapped_hq},${perc_mapped_hq},${nb_mapped_lq},${perc_mapped_lq},${nb_dups},${perc_dups},${nsc},${rsc},${frip} >> mqc.stats
done

