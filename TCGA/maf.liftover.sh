for f in *.maf.txt ; do 
cut -f5,6,7 ${f} |sed 's/^/chr/' |awk '{ if ($2==$3) print $1,$2,($3+1),"MODDED" ; else print $1,$2,$3 }' OFS="\t" |tail -n +2 > ${f}.bed
liftOver ./${f}.bed ~/misc_files/hg18ToHg19.over.chain ${f}.bed.lifted ${f}.unmapped.txt >> liftover.log 2>>liftover.log
sed -i 's/chr//' ${f}.bed.lifted 
sed -i '1 i Chromosome\tStart_position\tEnd_position' ${f}.bed.lifted 
#check line numbers are equal
var1=`wc -l  ${f}.bed.lifted |awk '{ print $1}'`
var2=`wc -l  ${f}|awk '{ print $1}'`
if [ $var1 == $var2 ] 
then 
	cut -f1,2,3 ${f} |awk '{ print $1,$2,$3,"37"}' OFS="\t"  | sed 's/Center\t37/Center\tNCBI_Build/' > temp1
	paste temp1 ${f}.bed.lifted  |awk '{ if ($8=="MODDED") print $1,$2,$3,$4,$5,$6,$6 ; else print $1,$2,$3,$4,$5,$6,$7}' OFS="\t"> temp2
	cut -f8- ${f} > temp3
	#keep old file
	mv ${f} ${f}.orig
	paste temp2 temp3 > ${f}  
else 
	echo "$f some missed liftovers"
	grep -v '#' ${f}.unmapped.txt |cut -f1,2 |sed 's/chr//' > ${f}.remove
	fgrep -v -f ${f}.remove ${f} > ${f}.nomissing
	cut -f1,2,3 ${f}.nomissing |awk '{ print $1,$2,$3,"37"}' OFS="\t"  | sed 's/Center\t37/Center\tNCBI_Build/' > temp1
	paste temp1 ${f}.bed.lifted  |awk '{ if ($8=="MODDED") print $1,$2,$3,$4,$5,$6,$6 ; else print $1,$2,$3,$4,$5,$6,$7}' OFS="\t"> temp2
	cut -f8- ${f}.nomissing > temp3
	#keep old file
	mv ${f} ${f}.orig
	paste temp2 temp3 > ${f}  
fi 
done 
