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
    parser.add_argument("D", metavar="INPUT_DESIGN", help="Enter a valid "
                                                        "design csv file")
    parser.add_argument("R", metavar="INPUT_READS_FILE", help="Enter a valid "
                                                        "reads csv file name")
    parser.add_argument("SE", metavar="SINGLE_END", help="Is data SE (True) or"
                                                        " PE (False) ?")
    parser.add_argument("BD", metavar="BASE DIR", help="Base dir if needed")

    args = parser.parse_args()
    inputDesign = args.D
    inputReads = args.R
    singleEnd = args.SE
    baseDir = args.BD
    return inputDesign, inputReads, singleEnd, baseDir


def check_designs(inputDesign, inputReads, singleEnd, baseDir):
    dict_design_keys = ['SAMPLEID', 'CONTROLID', 'SAMPLENAME', 'REPLICATE', 'PEAKTYPE']

    dict_design = {
        'SAMPLEID': [],
        'CONTROLID': [],
        'SAMPLENAME': [],
        'REPLICATE': [],
        'PEAKTYPE': []
    }
    if (singleEnd):
        dict_reads_keys = ['SAMPLEID', 'SAMPLENAME', 'FASTQR1']
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
        }
    else:
        dict_reads_keys = ['SAMPLEID', 'SAMPLENAME', 'FASTQR1', 'FASTQR2']
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
            'FASTQR2': []
        }
    ### Checks for design file
    with open(inputDesign, 'r') as designFile:
        lines = csv.reader(designFile)
        header = next(lines, None)
        # Check if design header is correct
        if header != dict_design_keys:
            print('Design file header is wrong, should be : {}'
                  .format(','.join(dict_design_keys)))
            sys.exit(1)
        # Fill dict to check all input design data
        for sample in lines:
            dict_design['SAMPLEID'].append(sample[0])
            dict_design['CONTROLID'].append(sample[1])
            dict_design['SAMPLENAME'].append(sample[2])
            dict_design['REPLICATE'].append(sample[3])
            dict_design['PEAKTYPE'].append(sample[4])
        # Check if samples and controls are correctly separated
        for ID in dict_design['CONTROLID']:
            if ID in dict_design['SAMPLEID']:
                print('The sample {} is both qualified as a control and a sample'
                      .format(ID))
                sys.exit(1)
        # Check if sample name does not contain replicate number
        for sampleName in dict_design['SAMPLENAME']:
            if re.search('-[0-9]$', sampleName):
                print('Sample name should not contain the replicate number ({}).'
                      ' Remove \'{}\''.format(sampleName, sampleName[-2:]))
        index = 0
        # Check if replicate number is an integer
        for replicateNumber in dict_design['REPLICATE']:
            try:
                int(replicateNumber)
            except:
                print('Replicate number of {} is not a valid number'
                      .format(dict_design['SAMPLEID'][index]))
                sys.exit(1)
            index += 1
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
    with open(inputReads, 'r') as readsFile:
        lines = csv.reader(readsFile)
        header = next(lines, None)
        # Check if sampleplan header is correct
        if header != dict_reads_keys:
            print('Sampleplan file header is wrong, should be : {}'
                  .format(','.join(dict_reads_keys)))
            sys.exit(1)
        # Fill dict to check all input sample data
        for sample in lines:
            dict_reads['SAMPLEID'].append(sample[0])
            dict_reads['SAMPLENAME'].append(sample[1])
            if sample[2][0] != '/':
                readfile = baseDir + '/' + sample[2]
                dict_reads['FASTQR1'].append(readfile)
            else:
                dict_reads['FASTQR1'].append(sample[2])
            if not singleEnd:
                if sample[3][0] != '/':
                    readfile = baseDir + '/' + sample[3]
                    dict_reads['FASTQR2'].append(readfile)
                else:
                    dict_reads['FASTQR2'].append(sample[3])
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
        # Check file extensions to match fastq ones
            if not ((samplePath.endswith('fq.gz'))
                or (samplePath.endswith('fastq.gz'))
                or (samplePath.endswith('fastq'))
                or (samplePath.endswith('fq'))):
                print('The file {} is not a fastq file'
                      .format(os.path.basename(samplePath)))
                sys.exit(1)


if __name__ == '__main__':
    inputDesign, inputReads, singleEnd, baseDir = argsParse()
    check_designs(inputDesign, inputReads, singleEnd, baseDir)