#!/usr/bin/env python
import argparse
import csv

def argsParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("i", metavar="INPUT_DESIGN", help="Enter a valid "
                                                          "design csv file")
    parser.add_argument("m", metavar="OUTPUT_TO_MAP", help="Enter a valid "
                                                           "csv file name")
    parser.add_argument("o", metavar="OUTPUT_DESIGN", help="Enter a valid "
                                                           "csv file name")
    parser.add_argument("se", metavar="SINGLE_END", help="Enter a boolean")                                                       

    args = parser.parse_args()
    inputDesign = args.i
    outputToMap = args.m
    outputDesign = args.o
    singleEnd = args.se
    return inputDesign, outputToMap, outputDesign, singleEnd


def prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd):
    dictDesign = {
        'sampleID': [],
        'method': [],
        'read1': [],
        'read2': [],
        'isInput': []
    }
    with open(inputDesign, 'r') as inpFile:
        lines = csv.reader(inpFile)
        for line in lines:
            dictDesign['sampleID'].append(line[0])
            dictDesign['method'].append(line[1].rsplit('_', 1)[1])
            dictDesign['read1'].append(line[2])
            if singleEnd == False:
                dictDesign['read2'].append(line[3])
            else:
                dictDesign['read2'].append('')
            if 'Input' in line[1]:
                dictDesign['isInput'].append('INPUT')
            else:
                dictDesign['isInput'].append('')
    with open(outputToMap, 'w') as toMapFile:
        header = 'SAMPLE_ID,FASTQ_R1,FASTQ_R2\n'
        toMapFile.write(header)
        for sampleNumber in range(len(dictDesign['sampleID'])):
            sampleID = dictDesign['sampleID'][sampleNumber]
            fastqR1 = dictDesign['read1'][sampleNumber]
            fastqR2 = dictDesign['read2'][sampleNumber]
            toMapFile.write(sampleID + ',' + fastqR1 + ',' + fastqR2 + '\n')
    with open(outputDesign, 'w') as outFile:
        header = 'SAMPLE_ID,CONTROL_ID,METHOD\n'
        outFile.write(header)
        for method in set(dictDesign['method']):
            for sampleNumber in range(len(dictDesign['sampleID'])):
                if dictDesign['method'][sampleNumber] == method:
                    if dictDesign['isInput'][sampleNumber] == 'INPUT':
                        ctrlSample = dictDesign['sampleID'][sampleNumber]
                        continue
                    else:
                        inputSample = dictDesign['sampleID'][sampleNumber]
                        outFile.write(inputSample + ',' + ctrlSample + \
                        ',' + method + '\n')


if __name__ == "__main__":
    inputDesign, outputToMap, outputDesign, singleEnd = argsParse()
    prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd)