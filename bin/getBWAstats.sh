#!/bin/bash 

## BWA mapping statistics

tot=$(samtools view -c -F 0x100 $1)
mapped=$(samtools view -F 0x4 -F 0x100 -c $1)
unmapped=$(( $tot - $mapped))
uniq=$(samtools view -q 1 -F 0x4 -F 0x100 $1 | grep -v XA:Z | grep -v SA:Z | wc -l )
multi=$(( $mapped - $uniq ))

echo -e "Unmapped\t${unmapped}" >> $2
echo -e "Uniquely mapped reads\t${uniq}" >> $2
echo -e "Multi mapped reads\t${multi}" >> $2
