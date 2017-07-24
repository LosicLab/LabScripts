
snpvcf=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/dbSNP_142.All.nocDNA.chr.vcf
myfiles=$1
mygenomefile=$2

mkdir -p log run data
while read line ; do
	ID=`echo $line |sed 's/.*\///' |sed 's/\..*//'`
echo "
#! /bin/bash
#BSUB -J intersect_"${ID}"
#BSUB -e log/"${ID}".e
#BSUB -o log/"${ID}".o
#BSUB -q low
#BSUB -W 0:30
#BSUB -n 2
#BSUB -P acc_PBG
###BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
module load bedtools
head -n 225 $line > data/${ID}.nosnp.vcf
fgrep -v chrUn $line |fgrep -v GL00 |fgrep -v KI27 > ${line}.std
bedtools intersect -sorted -v -g $mygenomefile -a ${line}.std -b $snpvcf >> data/${ID}.nosnp.vcf
" > run/${ID}.lsf

done < $myfiles
