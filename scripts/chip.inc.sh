## chip.inc.sh
##
## Copyright (c) 2017 Institut Curie
## Author(s): Aurélie Teissandier
## Contact: aurelie.teissandier@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details


## $1 = input files
## $2 = output dir
bowtie2_func()
{

    if [[ -z ${BOWTIE2_IDX_PATH} ]]; then
		die "indexes file not set. Exit"
    fi

    local out=$2/mapping
    mkdir -p ${out}

    #need to figure out --sam output ${bowtie_sam}
    prefix=$(basename $1 | sed -e 's/.fastq\(.gz\)*//')
    bowtie2_sam=${out}/$(basename $1 | sed -e 's/.fastq\(.gz\)*/.sam/')
    bowtie2_bam=${out}/$(basename $1 | sed -e 's/.fastq\(.gz\)*/.bam/')

    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
	cmd_in="<(gzip -cd ${inputs[0]})"
    elif [[ ${#inputs[@]} -eq 2 ]]; then
	cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 <(gzip -cd ${inputs[1]})"
    else
		die "Bowtie2 -  found more than two input files !"
    fi

    local cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_OPTS} -p 4 -x ${BOWTIE2_IDX_PATH} ${cmd_in} -S ${bowtie2_sam}"
    exec_cmd ${cmd}

    local cmd="${SAMTOOLS_PATH}/samtools view -bS ${bowtie2_sam} | ${SAMTOOLS_PATH}/samtools sort -@4 -T ${out}/${prefix} -o ${bowtie2_bam} -  "
    exec_cmd ${cmd}

    local cmd="rm $bowtie2_sam"
    exec_cmd ${cmd}
}


## $1 = input files
## $2 = output dir
rmDup_func()
{

	local out=$2/mapping 
	rmdup_bam=$(basename $1 ".bam")
	
	local cmd="${JAVA_PATH}/java -jar $PICARD_PATH/MarkDuplicates.jar I=$1 O=$out/${rmdup_bam}_noDup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=$out/${rmdup_bam}_metric"
	exec_cmd ${cmd}
	
	local cmd="${SAMTOOLS_PATH}/samtools index $out/${rmdup_bam}_noDup.bam"
    exec_cmd ${cmd}
   
}

## $1 = input files
## $2 = output dir
## $3 = label
fragSize_func()
{

	local out=$2/export 
	mkdir -p ${out}
	local cmd="${DEEPTOOLS_PATH}/bamPEFragmentSize --bamfiles $1 --histogram $out/bamPEFragmentSize.png --numberOfProcessors 4 --samplesLabel $3 > $2/logs/reportPEsize.log"
	exec_cmd ${cmd}
}

## $1 = input files
## $2 = output dir
## $3 = path to bedGraphToBigWig
bw_func()
{

	local out=$2/export 
	name=$(basename $1 "_noDup.bam")
	mkdir -p ${out}
	PATH=$3:$PATH
	local cmd="${DEEPTOOLS_PATH}/bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors 4 --normalizeUsingRPKM "
	exec_cmd ${cmd}
}

mapping_stat(){

    local output=$4/export
    mkdir -p $output
    inputs=($1)

    if [[ ${#inputs[@]} -eq 1 ]]; then
	cmd_input="-f ${inputs[0]}"

    elif [[ ${#inputs[@]} -eq 2 ]]; then
	cmd_input="-f ${inputs[0]} -r ${inputs[1]}"

    fi

    outfile=$(basename ${inputs[0]} | sed -e 's/.fastq\(.gz\)*/.stats/')

    if [ ! -z ${SAMPLE_ID} ]; then
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 -s ${SAMPLE_ID} > $output/$outfile"
    else
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 > $output/$outfile"
    fi
    exec_cmd ${cmd}
}



