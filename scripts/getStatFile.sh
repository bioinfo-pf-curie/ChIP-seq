## getStatFile.sh
##
## Copyright (c) 2017 Institut Curie                               
## Author(s):  Aurélie Teissandier
## Contact: aurelie.teissandier@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details


VERSION=0.0.1

function usage {
    echo -e "usage : ./getStatFile.sh -f FORWARD -b BAM -c CONFIG [-r REVERSE] [-s SAMPLE_ID]"
    echo -e "Use option -h|--help for more information"
}

function version {
    echo -e "version $VERSION"
    exit
}


function help {
    usage;
    echo
    echo "getStatFile.sh $VERSION"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -f FORWARD: input forward fastq file"
    echo "   -b BAM: bam file of aligned reads"
    echo "   -c CONFIG: configuration file for RNA processing"
    echo "   [-r REVERSE]: input reverse fastq file"
    echo "   [-s SAMPLE_ID]: biosample ID"
    echo "   [-h]: help"
    echo "   [-v]: version"
    exit;
}

while getopts "f:b:c:r:x:g:s:hv" OPT
do
    case $OPT in
        f) FORWARD=$OPTARG;;
	b) BAM=$OPTARG;;
	c) CONF=$OPTARG;;
	r) REVERSE=$OPTARG;;
	s) SAMPLE_ID=$OPTARG;;
        v) version ;;
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


if [[ -z $CONF  || -z $BAM ]]; then
    usage
    exit 1
fi

## Load config file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/utils.inc.sh
read_config $CONF

## sampleID
if [ ! -z ${FORWARD} ]; then
    if [ ! -z ${REVERSE} ]; then
	fastqID=$(basename $(printf "%s\n%s\n" "${FORWARD}}" "${REVERSE}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' ))
    else
	fastqID=$(basename ${FORWARD} | sed -e 's/.fastq\(.gz\)//')
    fi

    #cluster
    if [[ $FORWARD =~ "gz" ]];then
	nb_cluster=$(($(zcat ${FORWARD} | wc -l)/4))
    elif [[ $FORWARD =~ "fastq" ]];then
	nb_cluster=$(($(wc -l < ${FORWARD})/4))
    else
	die "ERROR : Wrong file type in input for fastq file; file: ${FORWARD}"
    fi


    if [ ! -z ${REVERSE} ]; then
	if [[ $REVERSE =~ "gz" ]];then
            nb_reverse=$(($(zcat ${REVERSE} | wc -l)/4))
	elif [[ $REVERSE =~ "fastq" ]];then
	    nb_reverse=$(($(wc -l < ${REVERSE})/4))
	else
	    die "ERROR : Wrong file type in input for fastq file; file: ${REVERSE}"  
	fi 
    else
	nb_reverse=0
    fi
    nb_reads=$((${nb_cluster} + ${nb_reverse}))
else
    nb_reads=NA
    nb_reads=NA
fi

#################
##
## Mapping Stats
##
#################

    
#https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
aligned=$(${SAMTOOLS_PATH}/samtools view -F 4 -c ${BAM})
ubam=$(${SAMTOOLS_PATH}/samtools view ${BAM} | awk '{if($0 ~! /^@/){if($0 ~ /XS:i:/){split($12,a,":");split($13,b,":");if(a[3] > b[3]){print}}else{print} }}' | wc -l )
mbam=$(($aligned - $ubam))

## verbose
echo -e "Sample_identifier,Biological_identifier,Number_of_cluster,Number_of_reads,Number_of_aligned_reads,Number_of_uniquely_aligned_reads,Number_of_multiple_aligned_reads"
echo -e ${fastqID},${SAMPLE_ID},${nb_cluster},${nb_reads},${aligned},${ubam},${mbam}



