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
    args = parser.parse_args()
    inputDesign = args.design
    inputData = args.sampleplan
    singleEnd = args.singleEnd
    return inputDesign, inputData, singleEnd

def loadSamplePlan(inputFile, isSingleEnd=False):
    """
    Load SamplePlan file with sampleId,sampelName,fastqR1,[fastqr2]
    """
    dictSamplePlan={'SAMPLEID':[], 'SAMPLENAME':[], 'FASTQR1':[]}
    if not isSingleEnd:
        dictSamplePlan['FASTQR2']=[]

    with open(inputFile, 'r') as dataFile:
        lines = csv.reader(dataFile)
        for sample in lines:
            dictSamplePlan['SAMPLEID'].append(sample[0])
            dictSamplePlan['SAMPLENAME'].append(sample[1])
            dictSamplePlan['FASTQR1'].append(sample[2])
            if not isSingleEnd:
                dictSamplePlan['FASTQR2'].append(sample[3])
    return(dictSamplePlan)


def loadDesign(inputFile, headers):
    """
    Load Design file using the defined headers
    """
    dictDesign = dict.fromkeys(headers, '')
    with open(inputFile, 'r') as designFile:
        lines = csv.reader(designFile)
        for sample in lines:
            for i in range(len(headers)):
                if dictDesign[headers[i]]=='':
                    dictDesign[headers[i]]=[]
                else:
                    dictDesign[headers[i]].append(sample[i])
    return(dictDesign)


def checkHeaders(inputDesign, headerDict):
    """
    Check headers on the design file
    """
    ### Checks for design file
    with open(inputDesign, 'r') as designFile:
        lines = csv.reader(designFile)
        header = next(lines)
        for i in range(0, len(header)):
            try:
                if not header[i] == [*headerDict][i]:
                    raise()
            except:
                print('\nError: Headers are not valid, should be : {}'
                      .format([*headerDict]))
                sys.exit(1)


def checkColumnContent(column, values):
    """
    Check the content of a column
    """
    for val in column:
        if not val in values:
            print('\nError: The value \'{}\' is invalid, should be : {}'
                  .format(val, [*values]))
            sys.exit(1)

def checkColumnsMatch(col1, col2, exclusive=False):
    """
    Check that values in col1 are (not) in col2
    """
    ## Remove empty values from col1/col2
    col1 = [i for i in col1 if i]
    col2 = [i for i in col2 if i]

    match=[]
    for ID in col1:
        if ID in col2:
            if exclusive:
                print('\nError: The value {} cannot be set in two columns'
                      .format(ID))
                sys.exit(1)
            else:
                if not ID in match:
                    match.append(ID)

    #print("\n")
    #print(set(col1))
    #print(set(col2))

    if not exclusive and len(set(col1).difference(match)) != 0:
        print('\nError: Values {} are not found in {}'
              .format(set(col1).difference(match), set(col2)))
        sys.exit(1)


if __name__ == '__main__':

    designHeader=['SAMPLEID', 'CONTROLID', 'SAMPLENAME', 'GROUP', 'PEAKTYPE']

    ## Get args
    inputDesign, inputSamplePlan, isSingleEnd = argsParse()
    
    ## Load SamplePlan
    print("[SAMPLEPLAN] Load data ", end='...')
    dictSamplePlan=loadSamplePlan(inputSamplePlan, isSingleEnd)
    print("ok") 

    ## Check Design headers
    print("[DESIGN] Check headers ", end='...')
    checkHeaders(inputDesign, designHeader)
    print("ok")

    ## Load Design
    print("[DESIGN] Load data ", end="...")
    dictDesign=loadDesign(inputDesign, designHeader)
    print("ok")

    ## Checks for design file
    print("[DESIGN] Check peak type content ", end='...')
    checkColumnContent(dictDesign['PEAKTYPE'], ['sharp', 'broad', 'very-broad'])
    print("ok")

    ## Check that a sample is not a control
    print("[DESIGN] Check samples/controls IDs ", end='...') 
    checkColumnsMatch(dictDesign['SAMPLEID'], dictDesign['CONTROLID'], exclusive=True)
    print("ok")

    ## Check that all samples from samplePlan are also in the design file (and the reverse)
    print("[DESIGN] Check samples matches between samplePlan and design ", end='...')
    checkColumnsMatch(dictSamplePlan['SAMPLEID'], dictDesign['SAMPLEID'] + dictDesign['CONTROLID'])
    checkColumnsMatch(dictDesign['SAMPLEID'] + dictDesign['CONTROLID'], dictSamplePlan['SAMPLEID'])
    print("ok")


