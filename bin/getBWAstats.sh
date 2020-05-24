#!/bin/bash 

## BWA mapping statistics

tot=$(samtools view -c -F 0x100 $1)
mapped=$(samtools view -F 0x4 -F 0x100 -c $1)
uniq=$(samtools view -q 1 -F 0x4 -F 0x100 $1 | grep -v XA:Z | grep -v SA:Z | wc -l )
multi=$(( $tot - $uniq ))

echo -e "Total reads= ${tot}" > $2
echo -e "Mapped reads= ${mapped}" >> $2
echo -e "Uniquely mapped reads= ${uniq}" >> $2
echo -e "Multi mapped reads= ${multi}" >> $2
