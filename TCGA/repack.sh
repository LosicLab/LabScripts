mkdir -p log
mkdir -p runfiles
while read line; do
	parentdir="$(dirname "$line")"
	sample=`echo $line |sed 's/.*\///'`
	#extract archive to directory
	echo "
#BSUB -J "$sample"
#BSUB -e log/"${sample}".e
#BSUB -o log/"${sample}".o
#BSUB -q alloc
#BSUB -W 4:00
#BSUB -n 1
#BSUB -P acc_immunotherapy
tar -zxvf "$line" -C "$parentdir"
#recompress the individual files
gzip -r -f "${parentdir}" && rm "${parentdir}"/*tar.gz
" > runfiles/unzip_${sample}.lsf
done <$1
