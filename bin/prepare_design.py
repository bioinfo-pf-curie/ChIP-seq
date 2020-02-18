#!/usr/bin/env python
import argparse
import csv

def argsParse():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("i", metavar="INPUT_DESIGN", help="Enter a valid "
                                                          "design csv file")
    parser.add_argument("m", metavar="OUTPUT_TO_MAP", help="Enter a valid "
                                                           "csv file name")
    parser.add_argument("o", metavar="OUTPUT_DESIGN", help="Enter a valid "
                                                           "csv file name")
    parser.add_argument("se", metavar="SINGLE_END", help="Is data SE or PE ?")

    args = parser.parse_args()
    inputDesign = args.i
    outputToMap = args.m
    outputDesign = args.o
    singleEnd = args.se
    return inputDesign, outputToMap, outputDesign, singleEnd


def prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd):
    """
    Generates two distinct files. A file destined to the alignment part, and
    a design file that will be used later in the pipeline to associate sample
    and control inputs for peak calling.
    """
    dictDesign = {
        'sampleID': [],
        'sampleName': [],
        'read1': [],
        'read2': [],
        'isInput': [],
		'replicate': [],
		'peaktype': []
    }
    with open(inputDesign, 'r') as inpFile:
        lines = csv.reader(inpFile)
        for line in lines:
            dictDesign['sampleID'].append(line[0])
            dictDesign['sampleName'].append(line[1][:-2])
            dictDesign['replicate'].append(line[1][-1])
            dictDesign['read1'].append(line[2])
            if singleEnd == False:
                dictDesign['read2'].append(line[3])
            else:
                dictDesign['read2'].append('')
            if 'Input' in line[1]:
                dictDesign['isInput'].append('INPUT')
                dictDesign['peaktype'].append('')
            else:
                dictDesign['isInput'].append('')
                mark = line[1].rsplit('_', 1)[0]
                if (mark == "H3K4me3"):
                    dictDesign['peaktype'].append('sharp')
                elif ((mark == "H3K27me3") or (mark == "H3K9me3")):
                    dictDesign['peaktype'].append('broad')
                elif (mark == "H3K9me2"):
                    dictDesign['peaktype'].append('very-broad')
    with open(outputToMap, 'w') as toMapFile:
        header = 'SAMPLE_ID,FASTQ_R1,FASTQ_R2\n'
        toMapFile.write(header)
        for sNumber in range(len(dictDesign['sampleID'])):
            sampleID = dictDesign['sampleID'][sNumber]
            fastqR1 = dictDesign['read1'][sNumber]
            fastqR2 = dictDesign['read2'][sNumber]
            toMapFile.write(sampleID + ',' + fastqR1 + ',' + fastqR2 + '\n')

    list_control = []
    list_name_control = []
    list_replicate = []
    sNumber = 0
    max_sample = len(dictDesign['sampleID'])-1
    while sNumber <= max_sample:
        if dictDesign['isInput'][sNumber] == 'INPUT':
            list_control.append(dictDesign['sampleID'][sNumber])
            list_name_control.append(dictDesign['sampleName'][sNumber])
            list_replicate.append(dictDesign['replicate'][sNumber])
            dictDesign['sampleID'].pop(sNumber)
            dictDesign['sampleName'].pop(sNumber)
            dictDesign['read1'].pop(sNumber)
            dictDesign['read2'].pop(sNumber)
            dictDesign['isInput'].pop(sNumber)
            dictDesign['replicate'].pop(sNumber)
            dictDesign['peaktype'].pop(sNumber)
            max_sample -= 1
        else:
            sNumber += 1
    with open(outputDesign, 'w') as designFile:
        header = 'SAMPLE_ID,CONTROL_ID,SAMPLENAME,REPLICATE,PEAKTYPE\n'
        designFile.write(header)
        for replicateNumber in range(len(list_replicate)):
            for sNumber in range(len(dictDesign['sampleID'])):
                if (list_replicate[replicateNumber] \
                    == dictDesign['replicate'][sNumber] \
                    and dictDesign['sampleName'][sNumber].rsplit('_', 1)[1] \
                    == list_name_control[replicateNumber].rsplit('_', 1)[1]):
                    sampleID = dictDesign['sampleID'][sNumber]
                    ctrl_sample = list_control[replicateNumber]
                    samplename = dictDesign['sampleName'][sNumber]
                    replicate = list_replicate[replicateNumber]
                    peaktype = dictDesign['peaktype'][sNumber]
                    designFile.write(sampleID + ',' \
                                     + ctrl_sample + ',' \
                                     + samplename + ',' \
                                     + replicate + ',' \
                                     + peaktype + '\n')

if __name__ == "__main__":
    inputDesign, outputToMap, outputDesign, singleEnd = argsParse()
    prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd)