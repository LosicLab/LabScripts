# $1 is a dir full of MAF files
# $2 is output dir.

for f in $1/*.maf.txt ; do 
	name=`basename $f |sed 's/.maf.txt/.vaf.txt/' `;
	echo $name 
	echo -e "Chromosome\tStart_Position\ttumor_ref_reads\ttumors_var_reads\ttumor_vaf" >${2}/${name}
	tail -n +2 $f |	cut -f5,6,43,44 | awk '{ print $1,$2,($3-$4),$4,(100*$4/$3) }' OFS="\t" >> ${2}/${name} ; 
done 



