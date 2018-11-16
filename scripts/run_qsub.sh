#!/bin/bash
## Nicolas Servant
## Run the pipeline for a set of sample using PBS -t option

cd $PBS_O_WORKDIR

## run:
## qsub -v SAMPLE_PLAN='./SAMPLE_PLAN',ODIR='/data/tmp/test_multi',CONFIG='CONFIG_star' run_qsub.sh

## Set PATHS
SCRIPTS_PATH=`dirname $0`
ABS_SCRIPTS_PATH=`cd "$SCRIPTS_PATH"; pwd`
BIN_PATH="$ABS_SCRIPTS_PATH/../bin/"

if [[ -z ${SAMPLE_PLAN} || -z ${ODIR} || -z ${CONFIG} ]];then
    echo -e "Error - Missing parameter(s). Stop."
    exit 1
fi

function getControlSample 
{

    local sample_id=$2
    forward=$(awk -F"," -v i=${sample_id} '$1==i{print $3}' ${SAMPLE_PLAN})
    reverse=$(awk -F"," -v i=${sample_id} '$1==i{print $4}' ${SAMPLE_PLAN})

    echo "F=$forward"
    echo "R=$reverse"

    if [[ ! -z ${reverse} && ! -z ${forward} ]]; then
	echo -e "-F $forward -R $reverse"
    elif [[ ! -z ${forward} ]]; then
	echo -e "-F $forward"
    else
	echo -e  "Control not found for $2 sample. Exit"
	exit 1
    fi
}


id=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $1}' ${SAMPLE_PLAN})
sampleid=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $2}' ${SAMPLE_PLAN})
forward=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $3}' ${SAMPLE_PLAN})
reverse=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $4}' ${SAMPLE_PLAN})
control=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $5}' ${SAMPLE_PLAN})
peakmode=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $6}' ${SAMPLE_PLAN})

mkdir -p ${ODIR}/${id}


echo -e "--------------------"
echo -e "Running the ChIP-seq pipeline on cluster mode (${VERSION})"
when=$(date +%Y-%m-%d)
echo -e "Date: ${when}"
where=$(hostname)
echo -e "Host: ${where}"
echo
echo -e "Id: ${id}"
echo -e "Sample_id: ${bioid}"
echo -e "Forward: ${forward}"
echo -e "Reverse:  ${reverse}"
echo -e "Control=${control}"
echo -e "Peakcalling=${peakmode}"
echo -e "Ouput: ${ODIR}"
echo -e "Config: ${CONFIG}"
echo -e "Log: ${ODIR}/${id}/chippip.log"

if [[ ! -z ${control} && ! -z ${peakmode} ]]; 
then
    opts=$(getControlSample ${control})
    opts="$opts -p ${peakmode}"
fi

## Run RNA pipeline per sample
if [ ! -z "$reverse" ]; then
    cmd="${BIN_PATH}/ChIPpip -f ${forward} -r ${reverse} -o ${ODIR}/${id} -c ${CONFIG} -s ${sampleid} ${opts} > ${ODIR}/${id}/chippip.log "
else
    cmd="${BIN_PATH}/ChIPpip -f ${forward} -o ${ODIR}/${id} -c ${CONFIG} -s ${sampleid} ${opts} > ${ODIR}/${id}/chippip.log"
fi
echo -e "Cmd: ${cmd}"
echo -e "--------------------"
echo
eval $cmd


