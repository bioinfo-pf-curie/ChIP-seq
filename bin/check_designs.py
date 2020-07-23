#!/usr/bin/env python

import argparse
import csv
import sys
import re
import os

def argsParse():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--design", dest="design", help="Design file (csv)", default=None)
    parser.add_argument("-s", "--sampleplan", dest="sampleplan", help="SamplePlan file (csv)")
    parser.add_argument("--singleEnd", help="Specify that input reads are single-end", action="store_true")
    parser.add_argument("--baseDir", help="Base dir if needed", default=".")
    parser.add_argument("--bam", help="Specify that input files are BAM files", action="store_true")
    args = parser.parse_args()
    inputDesign = args.design
    inputData = args.sampleplan
    singleEnd = args.singleEnd
    baseDir = args.baseDir
    inputBam = args.bam

    print(inputBam)

    return inputDesign, inputData, singleEnd, baseDir, inputBam


def check_designs(inputDesign, inputData, isSingleEnd, baseDir, isInputBam):
    dict_design = {
        'SAMPLEID': [],
        'CONTROLID': [],
        'SAMPLENAME': [],
        'GROUP': [],
        'PEAKTYPE': []
    }
    if (isSingleEnd):
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
        }
    else:
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
            'FASTQR2': []
        }

    ### Checks for design file
    if inputDesign is not None:
        with open(inputDesign, 'r') as designFile:
            lines = csv.reader(designFile)
            header = next(lines)
            for i in range(0, len(header)):
                try:
                    header[i] == [*dict_design][i]
                except:
                    print('Design file columns are not valid, should be : {}'
                          .format([*dict_design]))
                    sys.exit(1)
            # Fill dict to check all input design data
            for sample in lines:
                print(sample)
                dict_design['SAMPLEID'].append(sample[0])
                dict_design['CONTROLID'].append(sample[1])
                dict_design['SAMPLENAME'].append(sample[2])
                dict_design['GROUP'].append(sample[3])
                dict_design['PEAKTYPE'].append(sample[4])
            # Check if samples and controls are correctly separated
            for ID in dict_design['CONTROLID']:
                if ID in dict_design['SAMPLEID']:
                    print('The sample {} is both qualified as a control and a sample'
                        .format(ID))
                    sys.exit(1)
            # Check if peaktypes are correct for every sample
            peaktype_list = ['sharp', 'broad', 'very-broad']
            index = 0
            for peaktype in dict_design['PEAKTYPE']:
                if not peaktype in peaktype_list:
                    print('Peaktype for {} is invalid, can be : {}'
                        .format(dict_design['SAMPLEID'][index], 
                        ', '.join(peaktype_list)))
                index += 1
    ### Checks for sampleplan file
    with open(inputData, 'r') as dataFile:
        lines = csv.reader(dataFile)
        # Fill dict to check all input sample data
        for sample in lines:
            dict_reads['SAMPLEID'].append(sample[0])
            dict_reads['SAMPLENAME'].append(sample[1])
            if sample[2][0] != '/':
                readfile = baseDir + '/' + sample[2]
                dict_reads['FASTQR1'].append(readfile)
            else:
                dict_reads['FASTQR1'].append(sample[2])
            if not isSingleEnd:
                if sample[3][0] != '/':
                    readfile = baseDir + '/' + sample[3]
                    dict_reads['FASTQR2'].append(readfile)
                else:
                    dict_reads['FASTQR2'].append(sample[3])
        if inputDesign is not None:
            # Check if there is a missing ID in the design or sample file
            for ID in dict_reads['SAMPLEID']:
                if not (ID in dict_design['SAMPLEID'] or 
                        ID in dict_design['CONTROLID']):
                    print('Sample plan contains an ID that is not in the '
                        'design file ({})'.format(ID))
                    sys.exit(1)
            for ID in dict_design['SAMPLEID']:
                if not ID in dict_reads['SAMPLEID']:
                    print('Missing sample in design file ({})'.format(ID))
                    sys.exit(1)
            # Check if there is an input control given in the design
            if not '' in dict_design['CONTROLID']:
                for ID in dict_design['CONTROLID']:
                    if not ID in dict_reads['SAMPLEID']:
                        print('Missing input control in design file ({})'.format(ID))
                        sys.exit(1)
            # Check paths to files
        for samplePath in dict_reads['FASTQR1']:
            if not os.path.exists(samplePath):
                print('The path to {} does not lead to a file'
                      .format(os.path.basename(samplePath)))
                sys.exit(1)
        # Check file extensions to match fastq or sam/bam ones
            if(isInputBam == False):
                if not ((samplePath.endswith('fq.gz'))
                    or (samplePath.endswith('fastq.gz'))
                    or (samplePath.endswith('fastq'))
                    or (samplePath.endswith('fq'))):
                    print('The file {} is not a fastq file'
                          .format(os.path.basename(samplePath)))
                    sys.exit(1)
            else:
                if not  ((samplePath.endswith('sam'))
                    or  (samplePath.endswith('bam'))):
                    print('The file {} is not a sam/bam file'
                          .format(os.path.basename(samplePath)))
                    sys.exit(1)


if __name__ == '__main__':
    inputDesign, inputData, singleEnd, baseDir, inputBam = argsParse()
    check_designs(inputDesign, inputData, singleEnd, baseDir, inputBam)
