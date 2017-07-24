echo -e "barcode\tmixcrClones\tmixcrReads"
for f in data_mixcr/*Rescued.txt ; do 
	id=`echo $f |sed 's/_clones.Rescued.txt//' |sed 's/.*\///' `; 
	lines=`tail -n +2 $f |wc -l` ;
	if [ "$lines" -eq "0" ] ; then 
		echo -e $id"\t"$lines"\t"0 ;
	else 
		reads=`tail -n +2 $f |cut -f2 | paste -sd+ |bc ` ; 
		echo -e $id"\t"$lines"\t"$reads  
	fi
done 
