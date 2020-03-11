#!/usr/bin/env python
import argparse
import os
import pysam

AS_TAG = "XX"

def argsParse():
    """
    Parsing input & outputs BAM files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-IR", metavar="INPUT_BAM_REF", help="Enter a valid "
                                                        "BAM file")
    parser.add_argument("-IS", metavar="INPUT_BAM_SPIKE", help="Enter a valid "
                                                        "BAM file")
    parser.add_argument("-OR", metavar="OUTPUT_BAM_REF", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-OS", metavar="OUTPUT_BAM_SPIKE", help="Enter a valid "
                                                        "BAM file name")
    parser.add_argument("-SE", metavar="SINGLE_END", help="Is data SE (True) or"
                                                        " PE (False) ?")
    args = parser.parse_args()
    inputBamRef = args.IR
    inputBamSpike = args.IS
    outputBamRef = args.OR
    outputBamSpike = args.OS
    singleEnd = args.SE
    return inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd


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
    selected = ''
    # If both alignments unmapped
    if ((scores[0] == 1000) & (scores[1] == 1000)):
        alignment1.set_tag(AS_TAG,"UA")
        return alignment1

    # Reads mapped at the same position on the genome
    if (scores[0] < scores[1]): # Best score from parent 1
        alignment1.set_tag(AS_TAG,1)
        selected = 'reference'
        return alignment1, selected
    elif (scores[0] > scores[1]): # Best score from parent 2
        alignment2.set_tag(AS_TAG,1)
        selected = 'spike'
        return alignment2, selected
    else: # Same score
        if ((alignment1.pos == alignment2.pos) and (alignment1.rname == alignment2.rname)):
            # No allelic information on this read
            alignment1.set_tag(AS_TAG,"UA")
            selected = 'unmapped'
            return alignment1, selected
        else:
            # Ambiguous/Conflictual case
            alignment1.set_tag(AS_TAG,"CF")
            selected = 'ambiguous'
            return alignment1, selected

def compareRefSpike(inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd):
    if singleEnd:
        seq_type = 1
    else:
        seq_type = 2

    samInputRef = pysam.AlignmentFile(inputBamRef, 'rb')
    samInputSpike = pysam.AlignmentFile(inputBamSpike, 'rb')
    samOutputRef = pysam.AlignmentFile(outputBamRef, 'wb', template=samInputRef)
    samOutputSpike = pysam.AlignmentFile(outputBamSpike, 'wb', template=samInputSpike)

    counter_alignments = {
        'all': 0,
        'reference': 0,
        'spike': 0,
        'unmapped': 0,
        'ambiguous': 0,
    }

    comparison = 2
    prev_alignments = []

    for alignRef, alignSpike in zip(samInputRef.fetch(until_eof = True), samInputSpike.fetch(until_eof = True)):
        counter_alignments['all'] += 1
        pair_reads = [alignRef, alignSpike]
        scores = calc_scores(pair_reads, comparison)
        # First we need to process the previous group of alignments
        if (seq_type == 1): # Single-end
            # Select best position
            selected_alignment, selected = comp_two_single_alignments(pair_reads,scores)
            # Write alignment in output file
            if selected == 'unmapped':
                counter_alignments['unmapped'] += 1
            elif selected == 'ambiguous':
                counter_alignments['ambiguous'] += 1
            elif selected == 'reference':
                write_single_alignment(selected_alignment,samOutputRef)
                counter_alignments['reference'] += 1
            elif selected == 'spike':
                write_single_alignment(selected_alignment,samOutputSpike)
                counter_alignments['spike'] += 1
        # Then we process the alignment as it is a first alignment
    samInputRef.close()
    samInputSpike.close()
    samOutputRef.close()
    samOutputSpike.close()

    logName = os.path.basename(inputBamRef).rsplit('_',1)[0] + '_log.txt'
    with open(logName, 'w') as logFile:
        logFile.write('Total number of reads processed : ' + str(counter_alignments['all']) + '\n')
        logFile.write('Reads mapped on reference genome : ' + str(counter_alignments['reference']) + '\n')
        logFile.write('Reads mapped on spike genome : ' + str(counter_alignments['spike']) + '\n')
        logFile.write('Reads unmapped : ' + str(counter_alignments['unmapped']) + '\n')
        logFile.write('Ambiguous reads : ' + str(counter_alignments['ambiguous']) + '\n')

if __name__ == '__main__':
    inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd = argsParse()
    compareRefSpike(inputBamRef, inputBamSpike, outputBamRef, outputBamSpike, singleEnd)