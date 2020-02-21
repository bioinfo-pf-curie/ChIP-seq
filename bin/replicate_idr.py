#!/usr/bin/env python
import os
import argparse

def argsParse():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-sn", nargs='+', metavar="SAMPLE_NAMES", help="Enter a"
                                                " valid list of sample files")
    parser.add_argument("-pf", nargs='+', metavar="PEAK_FILES", help="Enter a"
                                                " valid list of peak files")

    args = parser.parse_args()
    sampleNames = args.sn
    str_peakFiles = args.pf
    return sampleNames, str_peakFiles


def IDR(sampleNames, str_peakFiles):
    sampleNames = ''.join(sampleNames).strip('[]').split(',')

    str_peakFiles = ''.join(str_peakFiles)
    peakFiles = [[ID for ID in sample.strip(" []").split(",")] for sample in str_peakFiles.strip('[]').split("],")]
    peaktype = peakFiles[0][1].rsplit('.', 1)[1]
    sampleNames_list = []
    peakFiles_list = []
    for sampleName in sampleNames:
        sampleNames_list.append(sampleName)
        peakfile_str = ""
        for peakFile in peakFiles:
            if peakFile[0] == sampleName:
                peakfile_str += peakFile[1] + ' '
        peakFiles_list.append(peakfile_str)
    index = 0
    for sampleName in sampleNames:
        os.system("idr --samples " + peakFiles_list[index]
                    + '--input-file-type ' + peaktype
                    + ' -o ' + sampleName + '_IDR.' + peaktype
                    + ' -l ' + sampleName + '_log.txt')
        index += 1

if __name__ == '__main__':
    sampleNames, str_peakFiles = argsParse()
    IDR(sampleNames, str_peakFiles)