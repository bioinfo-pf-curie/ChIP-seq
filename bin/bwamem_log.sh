filename=$1

awk 'NR == 1{s=$0} NR == 2 {s=s-$0} END {print s}' log.tmp >> log.tmp
awk 'NR==1{ printf "Total reads : "; print }' log.tmp > ${filename}.log
awk 'NR==2{ printf "Uniquely mapped reads : "; print }' log.tmp >> ${filename}.log
awk 'NR==3{ printf "Multiply mapped reads : "; print }' log.tmp >> ${filename}.log
