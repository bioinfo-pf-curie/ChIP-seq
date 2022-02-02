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

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Number_reads_after_filt,Percent_reads_after_filt,Number_of_duplicates,Percent_of_duplicates,Normalized_strand_correlation,Relative_strand_correlation,Fraction_of_reads_in_peaks"

for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample," $splan | awk -F, '{print $2}')

    #ALIGNMENT
    if [ $aligner == "bowtie2" ]; then
	nb_frag=$(grep "reads;" mapping/${sample}_bowtie2.log | sed 's/ .*//')
	if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	else
            nb_reads=$nb_frag
	fi
    elif [ $aligner == "bwa-mem" ]; then
	# bwa.log file is in reads number (not pairs)
	nb_reads=$(grep 'Total' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	if [[ $is_pe == 1 ]]; then
	    nb_frag=$(( $nb_reads / 2 ))
	else
	    nb_frag=$nb_reads
	fi
	tail -n +3 mapping/${sample}_bwa.log > mapping/${sample}_bwa.mqc
    elif [ $aligner == "star" ]; then
	nb_frag=$(grep "Number of input reads" mapping/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	else
            nb_reads=$nb_frag
	fi
    fi

    #Mapping stats (always in reads - so must be converted for PE)
    #These statistics are calculated after spike cleaning but before filtering
    nb_mapped=$(grep "mapped (" mapping/${sample}*.flagstats | awk '{print $1}')
    nb_filter=$(grep "mapped (" filtering/${sample}*.flagstats | awk '{print $1}')
    #nb_mapped=$(awk -F, '$1=="Mapped"{print $2}' mapping/${sample}_mappingstats.mqc)
    #nb_mapped_hq=$(awk -F, '$1=="HighQual"{print $2}' mapping/${sample}_mappingstats.mqc)
    #nb_mapped_lq=$(awk -F, '$1=="LowQual"{print $2}' mapping/${sample}_mappingstats.mqc)
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_filter=$(echo "${nb_filter} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    #perc_mapped_hq=$(echo "${nb_mapped_hq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    #perc_mapped_lq=$(echo "${nb_mapped_lq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    #PICARD
    if [[ -e filtering/${sample}.markDups_metrics.txt ]]; then
	nb_dups_pair=$(grep -a2 "## METRICS" filtering/${sample}.markDups_metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	nb_dups_single=$(grep -a2 "## METRICS" filtering/${sample}.markDups_metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	nb_dups_optical=$(grep -a2 "## METRICS" filtering/${sample}.markDups_metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	nb_dups=$(( $nb_dups_pair * 2 + $nb_dups_single + $nb_dups_optical ))
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
	if [ $(ls peakCalling/${sample}*FRiP.tsv 2>/dev/null | wc -l) -gt 0 ]; then
	    frip=''
	    for i in $(ls peakCalling/${sample}*FRiP.tsv); do
		frip=$(grep "$sample" $i | awk '{print $2}')"|${frip}"
	    done
	    frip=$(echo $frip | sed -e 's/|$//')
	else
	    frip='NA'
	fi
    else
	frip='NA'
    fi

    #To file
    echo -e ${sample},${sname},${nb_frag},${frag_length},${nb_mapped},${perc_mapped},${nb_filter},${perc_filter},${nb_dups},${perc_dups},${nsc},${rsc},${frip}
done

