#!/usr/bin/env python
import os
import argparse
import re

def argsParse():
    """
    Parsing all sample names and peakfiles to perform IDR on multiple peak 
    files for a given sample name
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-sn", nargs='+', metavar="SAMPLENAMES", help="Enter a"
                                                " valid list of sample files")
    parser.add_argument("-pf", nargs='+', metavar="PEAKFILES", help="Enter a"
                                                " valid list of peak files")

    args = parser.parse_args()
    sampleNames = args.sn
    str_peakFiles = args.pf
    return sampleNames, str_peakFiles


def IDR(sampleNames, str_peakFiles):
    """
    Computes IDR on each sample name that has replicates in the experiment
    """
    sampleNames = ''.join(sampleNames).strip('[]').split(',')

    str_peakFiles = ''.join(str_peakFiles)
    peakFiles = [[ID for ID in sample.strip(" []").split(",")] for sample 
                  in str_peakFiles.strip('[]').split("],")]

    sampleNames_list = []
    peakFiles_list = []
    peakType_list = []

    for sampleName in sampleNames:
        sampleNames_list.append(sampleName)
        peakfile_str = ""
        for peakFile in peakFiles:
            if peakFile[0] == sampleName:
                peakfile_str += peakFile[1] + ' '
            if peakFile[1].rsplit('.',1)[1] not in peakType_list:
                peakType_list.append(peakFile[1].rsplit('.',1)[1])
        peakFiles_list.append(peakfile_str)
    index = 0
    for sampleName in sampleNames:
        peaktype = re.search("\.[a-zA-Z]*Peak", peakFiles_list[index])
        peaktype = peaktype.group()[1:]
        if peakFiles_list[index].count(peaktype) > 1:
            os.system("idr --samples " + peakFiles_list[index]
                      + '--input-file-type ' + peaktype
                      + ' -o ' + sampleName + '_IDR.' + peaktype
                      + ' -l ' + sampleName + '_log.txt')
        index += 1

if __name__ == '__main__':
    sampleNames, str_peakFiles = argsParse()
    IDR(sampleNames, str_peakFiles)

 