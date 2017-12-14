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
    exit -1
fi

id=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $1}' ${SAMPLE_PLAN})
forward=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $3}' ${SAMPLE_PLAN})
reverse=$(awk -F"," -v i=${PBS_ARRAYID} 'NR==i{print $4}' ${SAMPLE_PLAN})

mkdir -p ${ODIR}/${id}

echo "SAMPLE_PLAN=${SAMPLE_PLAN}"
echo "ODIR=${ODIR}"
echo "CONFIG=${CONFIG}"
echo "id=${id}"
echo "reverse=${reverse}"
echo "forward=${forward}"

## Run RNA pipeline per sample
if [ ! -z "$reverse" ]; then
    echo "Run ${id} in PE mode ..."
    echo "--$reverse--"
    ${BIN_PATH}/ChIPpip -f ${forward} -r ${reverse} -o ${ODIR}/${id} -c ${CONFIG} -s ${id} > ${ODIR}/${id}/rnapip.log 
else
    echo "Run ${id} in SE mode ..."
    ${BIN_PATH}/ChIPpip -f ${forward} -o ${ODIR}/${id} -c ${CONFIG} -s ${id} > ${ODIR}/${id}/rnapip.log
fi
