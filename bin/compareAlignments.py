#!/usr/bin/env python
import argparse
import os
import pysam
import re
import sys

def argsParse():
    """
    Parsing input & outputs BAM files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", metavar="INPUT_BAM_REF", help="Enter a valid "
                                                        "BAM file")
    parser.add_argument("-s", metavar="INPUT_BAM_SPIKE", help="Enter a valid "
                                                        "BAM file")
    parser.add_argument("-o", metavar="OUTPUT_BAM_REF", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-os", metavar="OUTPUT_BAM_SPIKE", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-se", metavar="SINGLE_END", help="Is data SE (True) or"
                                                        " PE (False) ?")
    args = parser.parse_args()
    inputBamRef = args.i
    inputBamSpike = args.s
    outputBamRef = args.o
    outputBamSpike = args.os
    singleEnd = args.se
    return inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd

def compareTwoAlignments (read1, read2):
    if read1.mapping_quality > read2.mapping_quality:
        return(read1)
    elif read1.mapping_quality < read2.mapping_quality:
        return(read2)
    else:
        return(None)


## Remove everything after "/" or " " in read's name
def get_read_name(read):
    name = read.query_name
    return re.split('/| ', name)[0]


## Get max Mapping Quality
def getMaxMapQ(reads):
    maxMapQ=-1
    for r in reads:
        if r.mapping_quality > maxMapQ:
            maxMapQ=r.mapping_quality
    return maxMapQ


## Compare sets of alignments
def compareAlignments(set1, set2, filterNonPrimary=True):

    if filterNonPrimary:
        for r in set1:
            if r.is_secondary or r.is_supplementary:
                set1.remove(r)
        for r in set2:
            if r.is_secondary or r.is_supplementary:
                set2.remove(r)

    if getMaxMapQ(set1) > getMaxMapQ(set2):
        #print (get_read_name(set1[0]),"----BAM1") 
        return(set1)
    elif getMaxMapQ(set1) < getMaxMapQ(set2):
        #print (get_read_name(set1[0]),"----BAM2") 
        return(set2)
    else:
        #print (get_read_name(set1[0]),"----UNKNOWN") 
        return None


def compareBams(bam1, bam2, obam1, obam2, singleEnd):
    
    reads_counter=0
    reads_bam1only=0
    reads_bam2only=0
    reads_bam1rescue=0
    reads_bam2rescue=0
    reads_unmapped=0
    reads_ambiguous=0
    startInit = None
    r1Aln=[]
    r2Aln=[]

    oambi = "ambiguousReads.txt"
        
    with pysam.Samfile(bam1, "rb") as hr1,  pysam.Samfile(bam1, "rb") as hr1Mm, pysam.Samfile(bam2, "rb") as hr2, pysam.Samfile(bam2, "rb") as hr2Mm:

        out1 = pysam.AlignmentFile(obam1, "wb", template=hr1)
        out2 = pysam.AlignmentFile(obam2, "wb", template=hr2)
        out3 = open(oambi, 'w')

        for r1, r1Mm, r2, r2Mm in zip(hr1.fetch(until_eof=True), hr1Mm.fetch(until_eof=True), hr2.fetch(until_eof=True), hr2Mm.fetch(until_eof=True)):

            ## First initialization - iterate Mm + 1
            if startInit is None:
                r1Mm = next(hr1Mm)
                r2Mm = next(hr2Mm)
                startInit=1

            if get_read_name(r1) == get_read_name(r2):
                if r1.is_unmapped == False:
                    r1Aln.append(r1)
                if r2.is_unmapped == False:
                    r2Aln.append(r2)

                ## Deal with multi-Hits
                while True:
                    ## Same reads found again - Can also is.secondary or is.supplementary ?
                    if get_read_name(r1) == get_read_name(r1Mm):
                        r1Aln.append(r1Mm)
                        try:
                            r1 = next(hr1)
                            r1Mm = next(hr1Mm)
                        except StopIteration:
                            break
                    else:
                        break

                while True:
                    if get_read_name(r2) == get_read_name(r2Mm):
                        r2Aln.append(r2Mm)
                        try:
                            r2 = next(hr2)
                            r2Mm = next(hr2Mm)
                        except StopIteration:
                            break
                    else:
                        break

                ## r1/r2 are now on the last Multi reads
                ## r1Mm/r2Mm are on the next read

                ## Compare two sets of alignments
                if len(r1Aln) == 0 and len(r2Aln) == 0:
                    reads_unmapped += 1

                elif len(r1Aln) > 0 and len(r2Aln) == 0:
                    reads_bam1only +=1
                    for r in r1Aln:
                        out1.write(r)

                elif len(r1Aln) == 0 and len(r2Aln) > 0:
                    reads_bam2only +=1
                    for r in r2Aln:
                        out2.write(r)

                elif len(r1Aln) > 0 and len(r2Aln) > 0:
                    bestAlign = compareAlignments(r1Aln, r2Aln)
                    if bestAlign == r1Aln:
                        reads_bam1rescue +=1
                        for r in r1Aln:
                            out1.write(r)
                    elif bestAlign == r2Aln:
                        reads_bam2rescue +=1
                        for r in r2Aln:
                            out2.write(r)
                    elif bestAlign is None:
                        reads_ambiguous += 1
                        out3.write(get_read_name(r1) + "\n")
                
                ## Drop list elements
                r1Aln.clear()
                r2Aln.clear()
            
            else:
                print("Check that BAM files have the same read names and are sorted in the same way [", get_read_name(r1), "!=", get_read_name(r2),"]")
                sys.exit(1)       

    ## Last occurence of r1/r2
    if get_read_name(r1) != get_read_name(r1Mm) and get_read_name(r2) != get_read_name(r2Mm):
        if r1Mm.is_unmapped == False:
            r1Aln.append(r1Mm)
        if r2.is_unmapped == False:                                                        r2Aln.append(r2Mm)
        if len(r1Aln) == 0 and len(r2Aln) == 0:
            reads_unmapped += 1

        elif len(r1Aln) > 0 and len(r2Aln) == 0:
            reads_bam1only +=1
            for r in r1Aln:
                out1.write(r)
            
        elif len(r1Aln) == 0 and len(r2Aln) > 0:
            reads_bam2only +=1
            for r in r2Aln:
                out2.write(r)
            
        elif len(r1Aln) > 0 and len(r2Aln) > 0:
            bestAlign = compareAlignments(r1Aln, r2Aln)
            if bestAlign == r1Aln:
                reads_bam1rescue +=1
                for r in r1Aln:
                    out1.write(r)
            elif bestAlign == r2Aln:
                reads_bam2rescue +=1
                for r in r2Aln:
                    out2.write(r)
            elif bestAlign is None:
                reads_ambiguous += 1
                out3.write(get_read_name(r1) + "\n")
  
    
    hr1.close()
    hr1Mm.close()
    hr2.close()
    hr2Mm.close()

    logName = os.path.basename(inputBamRef).rsplit('_',1)[0] + '_log.txt'
    with open(logName, 'w') as logFile:
        logFile.write('Reads on BAM1 only : ' + str(reads_bam1only) + '\n')
        logFile.write('Reads rescue on BAM1 : ' + str(reads_bam1rescue) + '\n')
        logFile.write('Reads on BAM2 only : ' + str(reads_bam2only) + '\n')
        logFile.write('Reads rescue on BAM2 : ' + str(reads_bam2rescue) + '\n')
        logFile.write('Reads unmapped : ' + str(reads_unmapped) + '\n')
        logFile.write('Ambiguous reads : ' + str(reads_ambiguous) + '\n')

     
if __name__ == '__main__':
    inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd = argsParse()
    compareBams(inputBamRef,inputBamSpike,outputBamRef,outputBamSpike,singleEnd)
