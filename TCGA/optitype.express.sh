echo -e "barcode\tlocus\tallele\tsharedReads\tuniqueReads" ; 
for f in data1k/* ; do id=`echo $f |sed 's/.*\///' `; tail -n +2 ${f}/HLA.counts |cut -f2,3,4 |tr "*" "\t" |while read line ;do echo $id $line  ; done ; done |tr " " "\t"

