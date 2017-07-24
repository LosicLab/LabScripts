echo -e "ID\tchr\tpos\tVAF\tdepth\tclone" 
for f in scicloneOut_noCNV/* ; do id=`echo $f |sed 's/.*\///' |sed 's/.vaf.txt//' `; tail -n +2 $f |cut -f1,2,5,8,10 |while read line ;do  echo -e $id"\t"$line ; done ;done |tr " " "\t" 

