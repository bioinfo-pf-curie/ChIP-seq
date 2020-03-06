#!/usr/bin/env python
import argparse
import csv
import sys
import re
import os
import pysam

AS_TAG = "XX"

def argsParse():
    """
    Parsing input & outputs BAM files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-I", metavar="INPUT_BAM", help="Enter a valid "
                                                        "BAM file")
    parser.add_argument("-OR", metavar="OUTPUT_BAM_REF", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-OS", metavar="OUTPUT_BAM_SPIKE", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-SE", metavar="SINGLE_END", help="Is data SE (True) or"
                                                        " PE (False) ?")
    # parser.add_argument("-M", metavar="METHOD", help="Comparison method")
    args = parser.parse_args()
    inputBam = args.I
    outputBamRef = args.OR
    outputBamSpike = args.OS
    singleEnd = args.SE
    # method = args.M
    return inputBam, outputBamRef, outputBamSpike, singleEnd


def write_single_alignment(alignment,outbamFile):
    '''
    Write the alignment in the correct outbam file

    alignment: alignment [PYSAM object]
    all_outbams: all the different output BAM files [dict]
    counter_alignments: counters for the report [dict]
    '''

    if (alignment.get_tag(AS_TAG) != "CF"): # mapped alignment
        outbamFile.write(alignment)
        return
    else:
        return
    


def calc_scores(list_alignments,comparison):
    '''
    This function calculates scores of different alignments from a list
    The smaller the better:
      AS : Alignment Score (-c/--comparison 1)
       -> take absolute value
      NM : Number of mismatch (-c/--comparison 2)

    list_alignments: alignments to calculate score from [list]
    comparison: chosen parameter to compare (AS or NM) [int]

    RETURN scores associated to alignments (same index) [list]
    '''
    scores = []
    for alignment in list_alignments:
            if (alignment.flag >= 512):
                scores.append(1000)
            elif ((alignment.flag >= 4) and (str(bin(alignment.flag))[-3] == '1')):
                scores.append(1000)
            elif (comparison == 1):
                scores.append(abs(alignment.get_tag("AS")))
            elif (comparison == 2):
                scores.append(alignment.get_tag("NM"))
    return scores


def comp_two_single_alignments(list_alignments,scores):
    '''
    This function compares two alignments on scores where smaller is better
    and returns the best alignment with a TAG (AS_TAG):
        0: Equally mapped on both parents
        1: Preferentially mapped on one parent
        2: Ambiguous

    list_alignments: 2 alignments to compare [list]
    scores: 2 scores associated to alignments (same index) [list]
    '''
    alignment1 = list_alignments[0]
    alignment2 = list_alignments[1]
 
    # If both alignments unmapped
    if ((scores[0] == 1000) & (scores[1] == 1000)):
        alignment1.set_tag(AS_TAG,"UA")
        return alignment1

    # Reads mapped at the same position on the genome
    if (scores[0] < scores[1]): # Best score from parent 1
        alignment1.set_tag(AS_TAG,1)
        return alignment1
    elif (scores[0] > scores[1]): # Best score from parent 2
        alignment2.set_tag(AS_TAG,1)
        return alignment2
    else: # Same score
        if ((alignment1.pos == alignment2.pos) and (alignment1.rname == alignment2.rname)):
            # No allelic information on this read
            alignment1.set_tag(AS_TAG,"UA")
            return alignment1
        else:
            # Ambiguous/Conflictual case
            alignment1.set_tag(AS_TAG,"CF")
            return alignment1 


def comp_three_single_alignments(list_alignments,scores):
    '''
    This function compares three alignments on scores where smaller is better
    and returns the best alignment with a TAG (AS_TAG):
        0: Equally mapped on both parents
        1: Preferentially mapped on one parent
        2: Ambiguous

    list_alignments: 3 alignments to compare [list]
    scores: 3 scores associated to alignments (same index) [list]
    '''
    
    if (max(scores) == min(scores)): # Ambiguous case : at least 2 positions with same score
        alignment = list_alignments[0]
        alignment.set_tag(AS_TAG,"CF")
        return alignment
    else:
        # Get first minimum score
        alignment1 = list_alignments[scores.index(min(scores))]
        score1 = scores[scores.index(min(scores))]
        # Remove alignment and score from lists
        del list_alignments[scores.index(min(scores))]
        del scores[scores.index(min(scores))]
        alignment2 = list_alignments[scores.index(min(scores))]
        score2 = scores[scores.index(min(scores))]
        return comp_two_single_alignments([alignment1,alignment2],[score1,score2])


def comp_single_alignments(list_alignments,comparison):
    '''
    Function to select and return the best position for a alignment either:

    list_alignments: different lines for the alignment [list]
    comparison: parameter chosen for the comparison [int]
    '''
    scores = calc_scores(list_alignments,comparison)
    if (len(list_alignments) > 2):
    # Three alignments have been reported
        alignment = comp_three_single_alignments(list_alignments,scores)
    elif (len(list_alignments) > 1):
    # If only two alignments were reported
        alignment = comp_two_single_alignments(list_alignments,scores)
    else:
        alignment = list_alignments[0]
        if ((alignment.flag >= 4) and (str(bin(alignment.flag))[-3] == '1')):
            alignment.set_tag(AS_TAG,"UA")
        else:
            alignment.set_tag(AS_TAG,1)
    return alignment


def sepMetagenome(inputBam, outputBamRef, outputBamSpike, singleEnd):
    if singleEnd:
        seq_type = 1
    else:
        seq_type = 2

    samInput = pysam.AlignmentFile(inputBam, 'rb')
    samOutputRef = pysam.AlignmentFile(outputBamRef, 'wb', template=samInput)
    samOutputSpike = pysam.AlignmentFile(outputBamSpike, 'wb', template=samInput)

    counter_alignments = {
        'all': 0,
        'reference': 0,
        'spike': 0,
        'unmapped': 0,
        'ambiguous': 0,
    }

    comparison = 2
    prev_alignments = []

    for alignment in samInput.fetch(until_eof = True):
        counter_alignments['all'] += 1            
        if (len(prev_alignments) == 0):
            prev_alignments.append(alignment)
        else:
            if (alignment.query_name == prev_alignments[0].query_name):
                prev_alignments.append(alignment)
            else:
                # First we need to process the previous group of alignments
                if (seq_type == 1): # Single-end
                    # Select best position
                    selected_alignment = comp_single_alignments(prev_alignments,comparison)
                    refName = str(samInput.get_reference_name(selected_alignment.reference_id))
                    # Write alignment in output file
                    if selected_alignment.get_tag(AS_TAG) == 'UA':
                        counter_alignments['unmapped'] += 1
                    elif selected_alignment.get_tag(AS_TAG) == 'CF':
                        counter_alignments['ambiguous'] += 1
                    elif refName[0:3] == 'chr':
                        write_single_alignment(selected_alignment,samOutputRef)
                        counter_alignments['reference'] += 1
                    else:
                        write_single_alignment(selected_alignment,samOutputSpike)
                        counter_alignments['spike'] += 1
                prev_alignments = []
                # Then we process the alignment as it is a first alignment
                prev_alignments.append(alignment)
    samInput.close()
    samOutputRef.close()
    samOutputSpike.close()

    logName = inputBam.rsplit('.',1)[0] + '_log.txt'
    with open(logName, 'w') as logFile:
        logFile.write('Total number of reads processed : ' + str(counter_alignments['all']) + '\n')
        logFile.write('Reads mapped on reference genome : ' + str(counter_alignments['reference']) + '\n')
        logFile.write('Reads mapped on spike genome : ' + str(counter_alignments['spike']) + '\n')
        logFile.write('Reads unmapped : ' + str(counter_alignments['unmapped']) + '\n')
        logFile.write('Ambiguous reads : ' + str(counter_alignments['ambiguous']) + '\n')

if __name__ == '__main__':
    inputBam, outputBamRef, outputBamSpike, singleEnd = argsParse()
    sepMetagenome(inputBam, outputBamRef, outputBamSpike, singleEnd)