mkdir -p ${2}/HLA/ClassI
find $1 -name *result.tsv |while read line; do 
	id=`dirname $(dirname $line) | sed 's/.*\///'`; 
	echo $id ; 
	tail -n +2 $line |tr "\t" "\n" | tail -n +2 |head -n 6 > ${2}/HLA/ClassI/${id}.hla ; 
done
