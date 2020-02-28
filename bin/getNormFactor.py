#!/usr/bin/env python
import argparse
import os
import csv

def argsParse():
    """
    Parsing input flagstats and output csv
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", nargs='+', metavar="INPUT_FLAGSTATS", 
                                help="Enter a valid list of flagstat files")
    parser.add_argument("-O", metavar="OUTPUT_CSV", help="Enter a valid "
                                                        "csv file name")

    args = parser.parse_args()
    str_inputFlagstats = args.F
    outputCsv = args.O
    return str_inputFlagstats, outputCsv


def normFactors(str_inputFlagstats, outputCsv):
    """
    This script aims at calculating normalization factors for every sample
    of a ChIP-Seq experiment using spike-in analysis.
    """
    str_inputFlagstats = ''.join(str_inputFlagstats)
    inputFlagstats = [[flag for flag in flagfile.strip(" []").split(",")] for flagfile 
                  in str_inputFlagstats.strip('[]').split("],")]
    dictNormFactors = {
        'sampleID': [],
        'nbReads': [],
        'normFactor': []
    }
    
    for index in range(len(inputFlagstats)):
        dictNormFactors['sampleID'].append(os.path.basename(
                                                inputFlagstats[index][1])
                                                .rsplit('_',1)[0])
        with open(inputFlagstats[index][1], 'r') as flagfile:
            nbReads = flagfile.readline().rsplit(' + ',2)[0]
            dictNormFactors['nbReads'].append(int(nbReads))

    minReads = min(dictNormFactors['nbReads'])
    for nbReads in dictNormFactors['nbReads']:
        normFactor = str(round((minReads / nbReads), 2))[2:]
        dictNormFactors['normFactor'].append(normFactor)

    colnames = ['sampleID', 'nbReads', 'normFactor']
    with open(outputCsv, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(colnames)
        writer.writerows(zip(*[dictNormFactors[key] for key in colnames]))


if __name__ == '__main__':
    str_inputFlagstats,outputCsv = argsParse()
    normFactors(str_inputFlagstats, outputCsv)




