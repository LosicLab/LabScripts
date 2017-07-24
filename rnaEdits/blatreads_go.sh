## $1 is file with list of vcfs to filter.  

bamdir=/sc/orga/scratch/akersn01/full_HB/star_align/mergedbams
mkdir -p log run data
while read line ; do
        ID=`echo $line |sed 's/.*\///' |sed 's/\..*//'`
echo "
#! /bin/bash
#BSUB -J intersect_"${ID}"
#BSUB -e log/"${ID}".e
#BSUB -o log/"${ID}".o
#BSUB -q low
#BSUB -W 4:00
#BSUB -n 2
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=7000]" 
#BSUB -R "span[hosts=1]"
module load blat samtools
~/scripts/gatkvariants/rnaEdits/convert_vcf.pl ${line} data/${ID}.var
~/scripts/gatkvariants/rnaEdits/blat_candidates.pl data/${ID}.var ${bamdir}/${ID}.md.realn.bam data/${ID}.out

" >run/${ID}.lsf
done < $1
