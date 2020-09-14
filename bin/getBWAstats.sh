#!/bin/bash 

## BWA mapping statistics
## Note that the statistics are in number of reads (not pairs)

nb_pairs=$(samtools view -h $1 | head -n1000 | samtools view -f 0x1 -c -)

tot=$(samtools view -c -F 0x100 -F 0x800 $1)
mapped=$(samtools view -F 0x4 -F 0x100 -F 0x800 -c $1)
unmapped=$(( $tot - $mapped))
uniq=$(samtools view -q 1 -F 0x4 -F 0x100 -F 0x800 $1 | grep -v XA:Z | grep -v SA:Z | wc -l )
multi=$(( $mapped - $uniq ))

if [[ ${nb_pairs} -gt 0 ]]; then
    paired=$(samtools view -F 0x4 -F 0x100 -F 0x800 -f 0x1 -c $1)
    single=$(( $mapped - $paired ))
    uniq_paired=$(samtools view -q 1 -F 0x4 -F 0x100 -F 0x800 -f 0x1 $1 | grep -v XA:X | grep -v SA:Z | wc -l)
    uniq_single=$(( $uniq - $uniq_paired ))
    multi_paired=$(( $paired - $uniq_paired ))
    multi_single=$(( $single - $uniq_single ))
    unmapped=$(samtools view -f 12 $1)

    tot=$(( $tot / 2 ))
    mapped=$(( $mapped / 2))
    paired=$(( $paired / 2 ))
    uniq_paired=$(( $uniq_paired / 2 ))
    multi_paired=$(( $multi_paired / 2 ))
    unmapped=$(( $unmapped / 2 ))

    echo -e "Total\t${tot}" > $2
    echo -e "Mapped\t${mapped}" >> $2
    echo -e "PE neither mate aligned\t${unmapped}" >> $2
    echo -e "PE mapped uniquely\t${uniq_paired}" >> $2
    echo -e "PE one mate mapped uniquely\t${uniq_single}" >> $2
    echo -e "PE multi mapped\t${multi_paired}" >> $2
    echo -e "PE one mate multi\t${multi_single}" >> $2
else
    echo -e "Total\t${tot}" > $2
    echo -e "Mapped\t${mapped}" >> $2
    echo -e "Unmapped\t${unmapped}" >> $2
    echo -e "Uniquely mapped reads\t${uniq}" >> $2
    echo -e "Multi mapped reads\t${multi}" >> $2
fi
