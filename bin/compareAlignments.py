#!/usr/bin/env python
import argparse
import os
import pysam
import re
import sys
import time

def argsParse():
    """
    Parsing input & outputs BAM files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--bam1", help="Alignment file 1 (.bam)")
    parser.add_argument("-b", "--bam2", help="Alignment file 2 (.bam)")
    parser.add_argument("-oa", "--out1", help="Output file 1")
    parser.add_argument("-ob", "--out2", help="Output file 2")
    parser.add_argument("-c", "--comp", help="How to compare alignments (mapq, AS, NM)", default="AS")
    parser.add_argument("--debug",  help="Activate debug mode", action="store_true")

    args = parser.parse_args()
    return (args)


"""
Remove everything after "/" or " " in read's name
"""
def get_read_name(read):
    name = read.query_name
    return re.split('/| ', name)[0]


"""
Get min score on a list of alignment
"""
def getMinScore(reads, score="AS"):
    minScore=100000

    for r in reads:
        if score == "mapq":
            s = r.mapping_quality
        else:
            s = r.get_tag(score)
        if s < minScore:
            minScore = s
        return minScore

"""
Split a list of aligment into R1/R2 reads
"""
def splitPE(x):
    R1=list()
    R2=list()
    for r in x:
        if r.is_read1:
            R1.append(r)
        elif r.is_read2:
            R2.append(r)
    return R1,R2

"""
Compare sets of alignments
"""
def compareAlignments(set1, set2, filterNonPrimary=True, score="AS", debug=False):

    ## Remove unmapped reads
    g1 = [x for x in set1 if not x.is_unmapped]
    g2 = [x for x in set2 if not x.is_unmapped]

    if len(g1) == 0 and len(g2) == 0:
        return(-2)
    if len(g1) > 0 and len(g2) == 0:
        return(1)
    elif len(g1) == 0 and len(g2) > 0:
        return(2)

    ## Remove secondary alignments
    if filterNonPrimary:
        g1 = [x for x in g1 if not x.is_secondary or not x.is_supplementary]
        g2 = [x for x in g2 if not x.is_secondary or not x.is_supplementary]

    if getMinScore(set1, score=score) > getMinScore(set2, score=score):
        return(1)
    elif getMinScore(set1, score=score) < getMinScore(set2, score=score):
        return(2)
    else:
        return(0)

"""
Compare two BAM files
"""
def compareBams(bam1, bam2, out1, out2, score="mapq", debug=False):
    
    # Set up counters for report                                                                                                                                                                            
    counter_alignments = {}
    counter_alignments["bam1"] = 0
    counter_alignments["bam2"] = 0
    counter_alignments["unmapped"] = 0
    counter_alignments["ambiguous"] = 0
    counter_alignments["total"] = 0

    startInit = None
    isPE = False
    g1Aln=[]
    g2Aln=[]

    with pysam.Samfile(bam1, "rb") as hg1,  pysam.Samfile(bam1, "rb") as hg1Mm, pysam.Samfile(bam2, "rb") as hg2, pysam.Samfile(bam2, "rb") as hg2Mm:
        out1 = pysam.AlignmentFile(out1, "wb", template=hg1)
        out2 = pysam.AlignmentFile(out2, "wb", template=hg2)
        out3 = open("ambiguousReads.txt", 'w')

        for g1, g1Mm, g2, g2Mm in zip(hg1.fetch(until_eof=True), hg1Mm.fetch(until_eof=True), hg2.fetch(until_eof=True), hg2Mm.fetch(until_eof=True)):

            ## First initialization - iterate Mm + 1
            if startInit is None:
                g1Mm = next(hg1Mm)
                g2Mm = next(hg2Mm)
                startInit=1
                ## Check if data are paired
                if g1.is_paired and g2.is_paired:
                    isPE=True
                    if debug:
                        print ("[debug] Data are paired")

            if get_read_name(g1) == get_read_name(g2):
                counter_alignments['total'] += 1
                g1Aln.append(g1)
                g2Aln.append(g2)

                ## Deal with multi-Hits and/or PE data
                ## push in g1Aln/g2Aln until the read name changes
                while True:
                    if get_read_name(g1) == get_read_name(g1Mm):
                        g1Aln.append(g1Mm)
                        try:
                            g1 = next(hg1)
                            g1Mm = next(hg1Mm)
                        except StopIteration:
                            break
                    else:
                        break

                while True:
                    if get_read_name(g2) == get_read_name(g2Mm):
                        g2Aln.append(g2Mm)
                        try:
                            g2 = next(hg2)
                            g2Mm = next(hg2Mm)
                        except StopIteration:
                            break
                    else:
                        break

                ## g1/g2 are now on the last Multi reads
                ## g1Mm/g2Mm are on the next read name

                ## Compare two sets of alignments
                if debug:
                    print("[debug] ------------------------")
                    for r in g1Aln: 
                        print("[debug] BAM1:", r.query_name, r.reference_name, r.reference_start)
                    for r in g2Aln:
                        print("[debug] BAM2:", r.query_name, r.reference_name, r.reference_start)

                if isPE:
                    g1PE = splitPE(g1Aln)
                    g2PE = splitPE(g2Aln)

                    flagR1 = compareAlignments(g1PE[0], g2PE[0], score=score, debug=debug)
                    flagR2 = compareAlignments(g1PE[1], g2PE[1], score=score, debug=debug) 
                    flag = flagR1 + flagR2

                    if flag == -4:
                        counter_alignments["unmapped"] += 1
                    elif flag == 2 or flag == -1:
                        counter_alignments["bam1"] += 1
                        for r in g1Aln:
                            out1.write(r)
                    elif flag == 4 and flag == 0:
                        counter_alignments["bam2"] += 1
                        for r in g2Aln:
                            out2.write(r)
                    else:
                        counter_alignments["ambiguous"] += 1 
                        out3.write(get_read_name(g1) + "\n")
                    
                else:
                    flag = compareAlignments(g1Aln, g2Aln, score=score, debug=debug)
                    if flag == -2:
                        counter_alignments["unmapped"] += 1
                    elif flag == 1:
                        counter_alignments["bam1"] += 1
                        for r in g1Aln:
                            out1.write(r)  
                    elif flag == 2:
                        counter_alignments["bam2"] += 1
                        for r in g2Aln:
                            out2.write(r)
                    else:
                        counter_alignments["ambiguous"] += 1
                        out3.write(get_read_name(g1) + "\n")

                if debug:
                    print("[debug] FLAG:", flag)
                    print("\n")

                ## Drop list elements
                g1Aln.clear()
                g2Aln.clear()
            
            else:
                print("Check that BAM files have the same read names and are sorted in the same way [", get_read_name(g1), "!=", get_read_name(g2),"]")
                sys.exit(1)       

    ## Last occurence of g1/g2
    if get_read_name(g1) != get_read_name(g1Mm) and get_read_name(g2) != get_read_name(g2Mm):
        g1Aln.append(g1Mm)
        g2Aln.append(g2Mm)

        if isPE:
            g1PE = splitPE(g1Aln)
            g2PE = splitPE(g2Aln)

            flagR1 = compareAlignments(g1PE[0], g2PE[0], score=score, debug=debug)
            flagR2 = compareAlignments(g1PE[1], g2PE[1], score=score, debug=debug)
            flag = flagR1 + flagR2
            
            if flag == -4:
                counter_alignments["unmapped"] += 1
            elif flag == 2 or flag == -1:
                counter_alignments["bam1"] += 1
                for r in g1Aln:
                    out1.write(r)
            elif flag == 4 and flag == 0:
                counter_alignments["bam2"] += 1
                for r in g2Aln:
                    out2.write(r)
            else:
                counter_alignments["ambiguous"] += 1
                out3.write(get_read_name(r1) + "\n")   

        else:
            flag = compareAlignments(g1Aln, g2Aln, score=score, debug=debug)
            if flag == -2:
                counter_alignments["unmapped"] += 1
            elif flag == 1:                                                                                                                                                                            
                counter_alignments["bam1"] += 1
                for r in g1Aln:
                    out1.write(r)
            elif flag == 2:
                counter_alignments["bam2"] += 1
                for r in g2Aln:
                    out2.write(r)
            else:
                counter_alignments["ambiguous"] += 1
                out3.write(get_read_name(r1) + "\n")
 
            if debug:
                print("[debug] FLAG:", flag)
                print("\n")
  
    hg1.close()
    hg1Mm.close()
    hg2.close()
    hg2Mm.close()

    logName = os.path.basename(bam1).rsplit('_',1)[0] + '_bamcomp.log'
    with open(logName, 'w') as logFile:
        logFile.write("Ref BAM\t" + bam1 + "\n")
        logFile.write("Spike BAM\t" + bam2 + "\n")
        if isPE:
            logFile.write('Data type\tPaired-end\n')
        else:
            logFile.write('Data type\tSingle-end\n')
        logFile.write("Score\t" + score + "\n")
        logFile.write('Reads on ref\t' + str(counter_alignments["bam1"]) + '\n')
        logFile.write('Reads on spike\t' + str(counter_alignments["bam2"]) + '\n')
        logFile.write('Reads unmapped\t' + str(counter_alignments["unmapped"]) + '\n')
        logFile.write('Ambiguous reads\t' + str(counter_alignments["ambiguous"]) + '\n')

if __name__ == '__main__':
    args = argsParse()
    compareBams(bam1=args.bam1, bam2=args.bam2, out1=args.out1, out2=args.out2, 
                score=args.comp, debug=args.debug)
