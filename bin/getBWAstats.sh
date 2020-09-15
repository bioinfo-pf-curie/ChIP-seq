#!/bin/bash 

## BWA mapping statistics
## Note that the statistics are in number of reads (not pairs)

function usage {
    echo -e "usage : getBWAstats.sh -i INPUT [-p PROC] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "getBWAstats.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i INPUT"
    echo "   [-p] PROC"
    echo "   [-h]: help"
    exit;
}

proc=1
while getopts "i:p:h" OPT
do
    case $OPT in
        i) input=$OPTARG;;
	p) proc=$OPTARG;;
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


##########################################################
nb_pairs=$(samtools view -h $input | head -n1000 | samtools view -@ $proc -f 0x1 -c -)
tot=$(samtools view -@ $proc -F 0x100 -F 0x800 -c $input)
mapped=$(samtools view -@ $proc -F 0x4 -F 0x100 -F 0x800 -c $input)
unmapped=$(( $tot - $mapped))
uniq=$(samtools view -@ $proc -q 1 -F 0x4 -F 0x100 -F 0x800 $input | grep -v XA:Z | grep -v SA:Z | wc -l )
multi=$(( $mapped - $uniq ))

if [[ ${nb_pairs} -gt 0 ]]; then
    paired=$(samtools view -@ $proc -F 0x4 -F 0x100 -F 0x800 -f 0x1 -c $input)
    single=$(( $mapped - $paired ))
    uniq_paired=$(samtools view -@ $proc -q 1 -F 0x4 -F 0x100 -F 0x800 -f 0x1 $input | grep -v XA:X | grep -v SA:Z | wc -l)
    uniq_single=$(( $uniq - $uniq_paired ))
    multi_paired=$(( $paired - $uniq_paired ))
    multi_single=$(( $single - $uniq_single ))
    unmapped=$(samtools view -@ $proc -f 12 -c $input)

    tot=$(( $tot / 2 ))
    mapped=$(( $mapped / 2))
    paired=$(( $paired / 2 ))
    uniq_paired=$(( $uniq_paired / 2 ))
    multi_paired=$(( $multi_paired / 2 ))
    unmapped=$(( $unmapped / 2 ))

    echo -e "Total\t${tot}" > $2
    echo -e "Mapped\t${mapped}" >> $2
    echo -e "PE neither mate aligned\t${unmapped}" >> $2
    echo -e "PE mapped uniquely\t${uniq_paired}" >> $2
    echo -e "PE one mate mapped uniquely\t${uniq_single}" >> $2
    echo -e "PE multi mapped\t${multi_paired}" >> $2
    echo -e "PE one mate multi\t${multi_single}" >> $2
else
    echo -e "Total\t${tot}" > $2
    echo -e "Mapped\t${mapped}" >> $2
    echo -e "Unmapped\t${unmapped}" >> $2
    echo -e "Uniquely mapped reads\t${uniq}" >> $2
    echo -e "Multi mapped reads\t${multi}" >> $2
fi
