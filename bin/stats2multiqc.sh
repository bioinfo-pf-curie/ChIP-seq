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
    echo "   -g GENOME"
    echo "   [-p]: paired-end mode"
    echo "   [-h]: help"
    exit;
}

is_pe=0
while getopts "s:d:a:g:ph" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
	d) design=$OPTARG;;
	a) aligner=$OPTARG;;
        g) genome=$OPTARG;;
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
n_header=0
for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample," $splan | awk -F, '{print $2}')
    header="Sample_id,Sample_name"
    output="${sample},${sname}"

    ##trimming
    if [[ -e trimming/${sample}.fastq.gz_trimming_report.txt ]]; then
        nb_frag=$(grep "Total reads processed" trimming/${sample}.fastq.gz_trimming_report.txt | awk '{print $NF}' | sed -e 's/,//g')
        nb_reads=$nb_frag
        nb_trimmed=$(grep "Reads with adapters" trimming/${sample}.fastq.gz_trimming_report.txt | awk '{print $4}' | sed -e 's/,//g')
        perc_trimmed=$(echo "${nb_trimmed} ${nb_frag}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_frag,Percent_trimmed"
        output+=",${nb_frag},${perc_trimmed}"
    elif [[ -e trimming/${sample}_1.fastq.gz_trimming_report.txt && -e trimming/${sample}_2.fastq.gz_trimming_report.txt ]]; then
        nb_frag=$(grep "Total reads processed" trimming/${sample}_1.fastq.gz_trimming_report.txt | awk '{print $NF}' | sed -e 's/,//g')
        nb_reads=$(( $nb_frag * 2 ))
        nb_trimmed_r1=$(grep "Reads with adapters" trimming/${sample}_1.fastq.gz_trimming_report.txt | awk '{print $4}' | sed -e 's/,//g')
        nb_trimmed_r2=$(grep "Reads with adapters" trimming/${sample}_2.fastq.gz_trimming_report.txt | awk '{print $4}' | sed -e 's/,//g')
        perc_trimmed=$(echo "${nb_trimmed_r1} ${nb_trimmed_r2} ${nb_frag}" | awk ' { printf "%.*f",2,($1+$2)*100/($3*2) } ')
        header+=",Number_of_frag,Percent_trimmed"
        output+=",${nb_frag},${perc_trimmed}"
    else
      #ALIGNMENT
      if [[ $aligner == "bowtie2" && -e mapping/${sample}_${genome}_bowtie2.log ]]; then
        nb_frag=$(grep "reads;" mapping/${sample}_${genome}_bowtie2.log | sed 's/ .*//')
	  if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	  else
            nb_reads=$nb_frag
	  fi
      elif [[ $aligner == "bwa-mem" && -e mapping/${sample}_${genome}_bwa.log ]]; then
	# bwa.log file is in reads number (not pairs)
        nb_reads=$(grep 'Total' mapping/${sample}_${genome}_bwa.log | awk -F "\t" '{print $2}')
	if [[ $is_pe == 1 ]]; then
	    nb_frag=$(( $nb_reads / 2 ))
	else
	    nb_frag=$nb_reads
	fi
	tail -n +3 mapping/${sample}_${genome}_bwa.log > mapping/${sample}_bwa.mqc
      elif [[ $aligner == "star" && -e mapping/${sample}_${genome}Log.final.out ]]; then
	nb_frag=$(grep "Number of input reads" mapping/${sample}_${genome}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	else
            nb_reads=$nb_frag
	fi
      else
        nb_reads=$(grep "total" mapping/${sample}_${genome}.flagstats | awk '{print $1}')
	if [[ $is_pe == 1 ]]; then
	    nb_frag=$(( $nb_reads / 2 ))
	else
	    nb_frag=$nb_reads
	fi
      fi
      header+=",Number_of_frag"
      output+=",${nb_frag}"
    fi

    #Mapping stats (always in reads - so must be converted for PE)
    #These statistics are calculated after spike cleaning but before filtering
    #Note that the mapped line (second) includes both primary+secondary alignment
    if [[ $is_pe == 1 ]]; then
	nb_paired_mapped=$(grep "with itself and mate mapped" mapping/${sample}_${genome}.flagstats | awk '{print $1}')
	nb_single_mapped=$(grep "singletons" mapping/${sample}_${genome}.flagstats | awk '{print $1}')
	nb_mapped=$(( $nb_paired_mapped + $nb_single_mapped ))
	nb_paired_filter=$(grep "with itself and mate mapped" filtering/${sample}_${genome}_filtered.flagstats | awk '{print $1}')
	nb_single_filter=$(grep "singletons" filtering/${sample}_${genome}_filtered.flagstats | awk '{print $1}')
	nb_filter=$(( $nb_paired_filter + $nb_single_filter ))
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	perc_filter=$(echo "${nb_filter} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_aligned_reads,Percent_of_aligned_reads,Number_reads_after_filt,Percent_reads_after_filt"
	output+=",${nb_mapped},${perc_mapped},${nb_filter},${perc_filter}"
    else
	nb_mapped=$(grep "primary mapped (" mapping/${sample}_${genome}.flagstats | awk '{print $1}')
	nb_filter=$(grep "primary mapped (" filtering/${sample}_${genome}.flagstats | awk '{print $1}')
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
	perc_filter=$(echo "${nb_filter} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_aligned_reads,Percent_of_aligned_reads,Number_reads_after_filt,Percent_reads_after_filt"
        output+=",${nb_mapped},${perc_mapped},${nb_filter},${perc_filter}"
    fi

    #SPIKE
    perc_spike='NA'
    if [[ -e mapping/${sample}_ref_bamcomp.mqc ]]; then
	nb_spike=$(grep spike mapping/${sample}_ref_bamcomp.mqc | awk -F"\t" '{print $2}')
	perc_spike=$(echo "${nb_spike} ${nb_frag}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_spike,Percent_of_spike"
        output+=",${nb_spike},${perc_spike}"
    fi

    #PICARD
    if [[ -e filtering/${sample}_${genome}_markDups_metrics.txt ]]; then
	nb_dups_pair=$(grep -a2 "## METRICS" filtering/${sample}_${genome}_markDups_metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	nb_dups_single=$(grep -a2 "## METRICS" filtering/${sample}_${genome}_markDups_metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	nb_dups_optical=$(grep -a2 "## METRICS" filtering/${sample}_${genome}_markDups_metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	nb_dups=$(( $nb_dups_pair * 2 + $nb_dups_single + $nb_dups_optical * 2 ))
	perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
        header+=",Number_of_duplicates,Percent_of_duplicates"
        output+=",${nb_dups},${perc_dups}"
    fi

    #Fragment size
    if [[ -e fragSize/${sample}_insert_size_metrics.txt ]]; then
      frag_length=$(grep -A2 "## METRIC" fragSize/${sample}_insert_size_metrics.txt | tail -n 1 | awk '{print $1}')
      header+=",Fragment_length"
      output+=",${frag_length}"
    fi

    #PPQT
    if [[ -e ppqt/${sample}.spp.out && -e ppqt/${sample}_spp_nsc_mqc.tsv ]]; then
	if [[ ${frag_length} == 'NA' ]]; then
	  frag_length=$(awk '{print $3}' ppqt/${sample}.spp.out | sed 's/,.*//')
          header+=",Fragment_length"
          output+=",${frag_length}"
	fi
	nsc=$(grep "$sample" ppqt/${sample}_spp_nsc_mqc.tsv | awk '{print $2}')  
	rsc=$(grep "$sample" ppqt/${sample}_spp_rsc_mqc.tsv | awk '{print $2}')
        header+=",Normalized_strand_correlation,Relative_strand_correlation"
        output+=",${nsc},${rsc}"
    fi

    #PeakCalling 
    if [[ ! -z $design ]];
    then
      frip='NA'
      if [ $(ls peakCalling/${sample}*FRiP.tsv 2>/dev/null | wc -l) -gt 0 ]; then
        frip=''
        for i in $(ls peakCalling/${sample}*FRiP.tsv); do
          frip=$(grep "$sample" $i | awk '{print $2}')"|${frip}"
        done
        frip=$(echo $frip | sed -e 's/|$//')
      fi
      header+=",Fraction_of_reads_in_peaks"
      output+=",${frip}"
    fi

    if [ $n_header == 0 ]; then
        echo -e $header
        n_header=1
    fi
    echo -e $output
done

