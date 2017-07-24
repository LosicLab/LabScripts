find . -name "*.hla" |while read f ;do  id=`echo $f |sed 's/.hla//' |sed 's/.*\///'` ; cat $f |while read line ; do echo -e $id"\t"$line ; done ; done  |tr "*" "\t"

