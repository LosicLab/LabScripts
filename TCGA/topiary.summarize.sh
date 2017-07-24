ls data/classI/*post.csv |head -1 |while read line; do head -1 $line |sed 's/#/ID/' |tr "," "\t" ; done 
for f in data/classI*/*post.csv ; do id=`echo $f |sed 's/.*\///' |sed 's/_[0-9]*.post.csv//' `; tail -n +2 $f |cut -f2- -d"," | while read line; do echo $id","$line ; done ; done |tr "," "\t"
