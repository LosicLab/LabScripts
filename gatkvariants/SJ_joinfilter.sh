#usage: SJ_joinfilter Nreads Nindivids
cat */SJ.out.tab | awk '{ if ($6=="0") print $0'} | awk '{ if ($1 !="M" && $1 !="MT" && $1 !="chrM" && $1 != "chrMT") print $0 '} \
| awk -v var=$1 '{ if ($7 >=var) print $1,$2,$3,$4 }' |sort |uniq -c |awk -v var=$2 '{if ($1 >=var) print $2,$3,$4,$5}' OFS="\t" 
#gather all SJ.out.tab, take the unnannotated junctions, skip mitochondrial, if read counts is >= 5 and and it's found in >= 2 individuals, run with it.  
