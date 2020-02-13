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
    parser.add_argument("se", metavar="SINGLE_END", help="Is data SE or PE ?")

    args = parser.parse_args()
    inputDesign = args.i
    outputToMap = args.m
    outputDesign = args.o
    singleEnd = args.se
    return inputDesign, outputToMap, outputDesign, singleEnd


def prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd):
    dictDesign = {
        'sampleID': [],
        'sampleName': [],
        'read1': [],
        'read2': [],
        'isInput': [],
        'mark': [],
		'peaktype': []
    }
    with open(inputDesign, 'r') as inpFile:
        lines = csv.reader(inpFile)
        for line in lines:
            dictDesign['sampleID'].append(line[0])
            dictDesign['sampleName'].append(line[1].rsplit('_', 1)[1])
            dictDesign['read1'].append(line[2])
            if singleEnd == False:
                dictDesign['read2'].append(line[3])
            else:
                dictDesign['read2'].append('')
            if 'Input' in line[1]:
                dictDesign['isInput'].append('INPUT')
                dictDesign['mark'].append('')
                dictDesign['peaktype'].append('')
            else:
                dictDesign['isInput'].append('')
                mark = line[1].rsplit('_', 1)[0]
                dictDesign['mark'].append(mark)
                if (mark == "H3K4me3"):
                    dictDesign['peaktype'].append('sharp')
                elif ((mark == "H3K27me3") or (mark == "H3K9me3")):
                    dictDesign['peaktype'].append('broad')
                elif (mark == "H3K9me2"):
                    dictDesign['peaktype'].append('very-broad')
    with open(outputToMap, 'w') as toMapFile:
        header = 'SAMPLE_ID,FASTQ_R1,FASTQ_R2\n'
        toMapFile.write(header)
        for sampleNumber in range(len(dictDesign['sampleID'])):
            sampleID = dictDesign['sampleID'][sampleNumber]
            fastqR1 = dictDesign['read1'][sampleNumber]
            fastqR2 = dictDesign['read2'][sampleNumber]
            toMapFile.write(sampleID + ',' + fastqR1 + ',' + fastqR2 + '\n')
    with open(outputDesign, 'w') as outFile:
        header = 'SAMPLE_ID,CONTROL_ID,SAMPLENAME,MARK,PEAK_TYPE\n'
        outFile.write(header)
        for sampleName in set(dictDesign['sampleName']):
            for sampleNumber in range(len(dictDesign['sampleID'])):
                if dictDesign['sampleName'][sampleNumber] == sampleName:
                    print(dictDesign['sampleName'][sampleNumber])
                    if dictDesign['isInput'][sampleNumber] == 'INPUT':
                        ctrlSample = dictDesign['sampleID'][sampleNumber]
                        continue
                    else:
                        inputSample = dictDesign['sampleID'][sampleNumber]
                        mark = dictDesign['mark'][sampleNumber]
                        peaktype = dictDesign['peaktype'][sampleNumber]
                        outFile.write(inputSample + ',' + ctrlSample + \
                        ',' + sampleName + ',' + mark + ',' + peaktype + '\n')


if __name__ == "__main__":
    inputDesign, outputToMap, outputDesign, singleEnd = argsParse()
    prepareDesign(inputDesign, outputToMap, outputDesign, singleEnd)