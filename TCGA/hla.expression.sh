for f in seq2hla/*/*expression ;do id=`echo $f |sed 's/.*\///' |sed 's/-Class.*//' ` ; cat $f |sed 's/://' |sed 's/RPKM//' |while read line; do echo -e $id"\t"$line ; done ; done |tr " " "\t" 
