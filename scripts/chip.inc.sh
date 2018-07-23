## chip.inc.sh
##
## Copyright (c) 2017 Institut Curie
## Author(s): Aurélie Teissandier
## Contact: aurelie.teissandier@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

check_env()
{
    if [ -z ${NB_PROC} ]; then
	export NB_PROC=2
    fi
    
    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/activate /bioinfo/local/build/Centos/envs_conda/deeptools_3.1.0
    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/activate /bioinfo/local/build/Centos/envs_conda/epic_0.2.9
}

isPEexperiment()
{
    nbpairs=$(samtools view -c -f 1 $1)
    if [[ $nbpairs == 0 ]]; then
	echo "0"
    else
	echo "1"
    fi
}


## $1 = input files
## $2 = output dir
bowtie2_func()
{
    check_env
    
    if [[ -z ${BOWTIE2_IDX_PATH} ]]; then
		die "indexes file not set. Exit"
    fi

    local out=$2/mapping
    mkdir -p ${out}
    local log=$2/logs
    mkdir -p ${log}

    #need to figure out --sam output ${bowtie_sam}
    prefix=$(basename $1 | sed -e 's/.fastq\(.gz\)*//')
    build=$(basename ${BOWTIE2_IDX_PATH})
    bowtie2_sam=${out}/$(basename $1 | sed -e 's/.fastq\(.gz\)*/_'${build}'.sam/')
    bowtie2_bam=${out}/$(basename $1 | sed -e 's/.fastq\(.gz\)*/_'${build}'.bam/')

    log=$log/mapping_${build}.log

    echo -e "Running Bowtie2 mapping on $build ..."
    echo -e "Logs: $log"
    echo


    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
	cmd_in="<(gzip -cd ${inputs[0]})"
    elif [[ ${#inputs[@]} -eq 2 ]]; then
	cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 <(gzip -cd ${inputs[1]})"
    else
	die "Bowtie2 - found more than two input files !"
    fi

    local cmd="bowtie2 ${BOWTIE2_OPTS} -p ${NB_PROC} -x ${BOWTIE2_IDX_PATH} ${cmd_in} -S ${bowtie2_sam}"
    exec_cmd ${cmd} > ${log} 2>&1

    local cmd="samtools view -bS ${bowtie2_sam} | samtools sort -@${NB_PROC} -T ${out}/${prefix} -o ${bowtie2_bam} -  "
    exec_cmd ${cmd} >> ${log} 2>&1

    local cmd="rm $bowtie2_sam"
    exec_cmd ${cmd} >> ${log} 2>&1
}

isPairedBam()
{
    nb_paired=$(samtools view -c -f 1 $1)
    if [[ $nb_paired -gt 0 ]]; then
	return 1
    else
	return 0
    fi
}

bam2fastq_func()
{ 
    check_env
    local bam=$1
    local filters=$3

    local out=$2/mapping
    mkdir -p ${out}
    local log=$2/logs
    mkdir -p ${log}
    log=$log/bam2fastq.log

    echo -e "Transform Bam to Fastq ..."
    echo -e "Logs: $log"
    echo

    prefix=$(basename $bam | sed -e 's/.bam//')
    if [[ $(isPairedBam $bam) == 1 ]]; then
	##Extract reads following filtering options
	local cmd="samtools view -b $filters $bam | bamToFastq -i stdin -fq ${out}/${prefix}_R1.fastq -fq2 ${out}/${prefix}_R2.fastq"
    else
	##Extract reads following filtering options
	local cmd="samtools view -b $filters $bam | bamToFastq -i stdin -fq ${out}/${prefix}_R1.fastq"
    fi
    exec_cmd ${cmd} >> ${log} 2>&1
}

## $1 = input files
## $2 = output dir
rmDup_func()
{

    check_env
    local out=$2/mapping 
    mkdir -p ${out}
    local log=$2/logs
    mkdir -p ${log}
    log=$log/rmdup.log

    echo -e "Remove duplicates ..."
    echo -e "Logs: $log"
    echo

    rmdup_bam=$(basename $1 ".bam")
    local cmd="java -jar ${PICARD_PATH}/MarkDuplicates.jar I=$1 O=$out/${rmdup_bam}_noDup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=$out/${rmdup_bam}_metric"
    exec_cmd ${cmd} > ${log} 2>&1
    
    local cmd="samtools index $out/${rmdup_bam}_noDup.bam"
    exec_cmd ${cmd} >> ${log} 2>&1
}

## $1 = input files
## $2 = output dir
## $3 = label
fragSize_func()
{
    check_env
    local out=$2/mapping 
    mkdir -p ${out}
    local log=$2/logs
    mkdir -p ${log}
    log=$log/fragsize.log
    
    echo -e "Infer fragment size ..."
    echo -e "Logs: $log"
    echo

    local cmd="bamPEFragmentSize --bamfiles $1 --histogram $out/bamPEFragmentSize.png --numberOfProcessors ${NB_PROC} --samplesLabel $3 > $2/logs/reportPEsize.log"
    exec_cmd ${cmd} > ${log} 2>&1
}

## $1 = input files
## $2 = output dir
## $3 = path to bedGraphToBigWig
bw_func()
{
    check_env
    local out=$2/tracks
    mkdir -p ${out}
    local log=$2/logs
    mkdir -p ${log}
    log=$log/bamCoverage.log

    echo -e "Generate bigwig file(s) ..."
    echo -e "Logs: $log"
    echo

    if [[ -z ${EFFECTIVE_GENOME_SIZE} ]]; then die "BigWig- EFFECTIVE_GENOME_SIZE is not defined. Stop."; fi

    name=$(basename $1 "_noDup.bam")
 
    if [[ ! -z ${ENCODE_BLACKLIST} && -e ${ENCODE_BLACKLIST} ]]; then
	local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM --blackListFileName ${ENCODE_BLACKLIST}"
    else
	local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM "
    fi
    exec_cmd ${cmd} > ${log} 2>&1
}

mapping_stat()
{
    check_env
    local output=$4/stats
    mkdir -p $output
    local log=$4/logs
    mkdir -p ${log}
    log=$log/mapping_stat.log
 
    echo -e "Run mapping stats ..."
    echo -e "Logs: $log"
    echo

    inputs=($1)

    if [[ $inputs != "" ]]; then
	if [[ ${#inputs[@]} -eq 1 ]]; then
	    cmd_input="-f ${inputs[0]}"

	elif [[ ${#inputs[@]} -eq 2 ]]; then
	    cmd_input="-f ${inputs[0]} -r ${inputs[1]}"
	    
	fi
	outfile=$(basename ${inputs[0]} | sed -e 's/.fastq\(.gz\)*/.stats/')
    else
	cmd_input=""
	if [[ -e $3 ]]; then
	    outfile=$(basename ${3} | sed -e 's/bam/stats/')
	else
	    die "No fastq file(s) nor BAM file specified in mapping stats. Stop"
	fi
    fi

    if [ ! -z ${SAMPLE_ID} ]; then
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 -s ${SAMPLE_ID} > $output/$outfile"
    else
	cmd="bash ${SCRIPTS_PATH}/getStatFile.sh $cmd_input -c $2 -b $3 > $output/$outfile"
    fi
    exec_cmd ${cmd} > ${log} 2>&1
}


frip()
{ 
    check_env
    local output=$3/qc
    mkdir -p $output
    local log=$3/logs
    mkdir -p ${log}
    log=$log/frip.log
 
    echo -e "calculatin FRIP ..."
    echo -e "Logs: $log"
    echo

    local bam=$1
    local peaks=$2

    if [[ ! -e $bam || ! -e $peaks ]]; then
	die "FRIP - input files not found !"
    fi

    prefix=$(basename $bam | sed -e 's/.bam//')
    if [ -z "${DRY_RUN+x}" ]; then
	tot=$(samtools view -F4 -c $bam)
	inpeaks=$(intersectBed -a $bam -b $peaks | samtools view -F4 -c -)
	echo -e $(basename $bam)"\t$tot\t$inpeaks" > $output/${prefix}_frip.stats
    fi
}


chipenrich_func()
{
    check_env
    local output=$2/qc
    mkdir -p $output
    local log=$2/logs
    mkdir -p ${log}
    log=$log/chipenrich.log

    echo -e "Calculate ChIP enrichment ..."
    echo -e "Logs: $log"
    echo

    local inputs_bam=($1)
    prefix=$(basename ${inputs_bam[0]} | sed -e 's/.bam//')

    cmd="plotFingerprint -b $1 \
--smartLabels --minMappingQuality 0 --skipZeros --binSize 500 --numberOfSamples 50000 \
-T \"ChIP enrichment against Input\" --plotFile ${output}/${prefix}_fingerprints.png \
--outRawCounts ${output}/${prefix}_fingerprints.tab -p ${NB_PROC}"
    exec_cmd $cmd > $log 2>&1  
}




## TODO - extract the fragment size
## TODO - and use no model option
macs_func()
{
    check_env
    local bam=$1
    local control=$2
    local odir=$3
    local mode=$4

    local output=$3/peaks
    mkdir -p $output
    local log=$3/logs
    mkdir -p ${log}
    log=$log/peak_calling.log
 
    echo -e "Run MACS2 peak calling ..."
    echo -e "Logs: $log"
    echo

    local name=$(basename $1 | sed -e 's/.bam//')
    if [[ -z ${EFFECTIVE_GENOME_SIZE} ]]; then die "Peak calling - EFFECTIVE_GENOME_SIZE is not defined. Stop."; fi
    if [[ ! -z ${SAMPLE_ID} ]]; then MACS_OPTS="$MACS_OPTS -n ${SAMPLE_ID}_macs2";fi


    if [[ ${mode} == "TF" ]]; then
	cmd="${MACS_PATH}/macs2 callpeak \
    -g ${EFFECTIVE_GENOME_SIZE} \
    -t ${bam} \
    -c ${control} \
    -n ${name} \
    ${MACS_OPTS} \
    --outdir ${odir}/peaks/ --tempdir ${odir}/peaks/"

    elif [[ ${mode} == "BROAD" ]]; then
	cmd="${MACS_PATH}/macs2 callpeak \
    -g ${EFFECTIVE_GENOME_SIZE} \
    -t ${bam} \
    -c ${control} \
    -n ${name} \
    --broad \
    ${MACS_OPTS} \
    --outdir ${odir}/peaks/ --tempdir ${odir}/peaks/"
    fi
    exec_cmd ${cmd} > ${log} 2>&1
}


epic_func()
{
    check_env
    local bam=$1
    local control=$2
    local odir=$3

    local output=$3/peaks
    mkdir -p $output
    local log=$3/logs
    mkdir -p ${log}
    log=$log/peak_calling_epic.log

    echo -e "Run EPIC peak calling ..."
    echo -e "Logs: $log"
    echo

    local prefix=$(basename $bam | sed -e 's/.bam//')
    local chipbed=$(echo $bam | sed -e 's/bam/bed/')
    cmd="bamToBed -i $bam > $chipbed"
    exec_cmd ${cmd} > ${log} 2>&1

    local controlbed=$(echo $control | sed -e 's/bam/bed/')
    cmd="bamToBed -i $control > $controlbed"
    exec_cmd ${cmd} >> ${log} 2>&1

    cmd="epic --treatment ${chipbed} --control ${controlbed} --number-cores ${NB_PROC} ${EPIC_OPTS} \
     --genome ${GENOME} --chromsizes ${CHROMSIZES} \
     --bed ${output}/${prefix}_results.bed --outfile ${output}/${prefix}_results.out "
    exec_cmd ${cmd} >> ${log} 2>&1

    if [ ! -z ${SAMPLE_ID} ]; then
	local trackline="track name='${SAMPLE_ID}_EPIC' description='${SAMPLE_ID}_EPIC' useScore=1"
	cmd="sed -i.bak 1i\"$trackline\" ${output}/${prefix}_results.bed"
	exec_cmd ${cmd} >> ${log} 2>&1
	rm ${output}/${prefix}_results.bed.bak
    fi

}


heatmap_func()
{
check_env

local bigwig=$1
local bedfiles=$2
local odir=$3

local output=$3/heatmaps
mkdir -p $output
local log=$3/logs
mkdir -p ${log}
log=$log/heatmaps.log

echo -e "Run Heatmaps on functional annotations ..."
echo -e "Logs: $log"

#for bed in ${bedfiles}
#do
#echo "- ${bed}"
#prefix=$(basename ${bigwig} | sed -e 's/.bw//')_$(basename ${bed} | sed -e 's/.bed//')
prefix=$(basename ${bigwig} | sed -e 's/.bw//')
 
cmd="computeMatrix reference-point -S ${bigwig} -R ${bedfiles} \
--referencePoint TSS --beforeRegionStartLength 2500 --afterRegionStartLength 2500 -skipZeros \
--outFileName ${odir}/heatmaps/${prefix}.mat.gz -p ${NB_PROC}"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="plotHeatmap --matrixFile ${odir}/heatmaps/${prefix}.mat.gz --outFileName ${odir}/heatmaps/${prefix}_heatmap.png \
 --colorMap Blues --alpha 0.8 \
--legendLocation upper-right --heatmapWidth 5 --dpi 300"
exec_cmd ${cmd} >> $log 2>&1
#done
echo

}

