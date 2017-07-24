FilesToRun=./HB_merge.list
picard=/hpc/packages/minerva-common/picard/1.93/bin
gatk=/hpc/packages/minerva-common/gatk-mssm/3.2.0/target
#ref=/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/Homo_sapiens.GRCh37.75.dna.primary_assembly.sort.fa
ref=/sc/orga/scratch/akersn01/tempref/GRCh38_Gencode21/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
vcfmills=/sc/orga/work/akersn01/ref/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
vcf1000g=/sc/orga/work/akersn01/ref/gatk/1000G_phase1.indels.hg19.sites.vcf
vcfdbsnp=/sc/orga/work/akersn01/ref/gatk/dbsnp_138.hg19.sort.vcf
#sampleguide=./HB_merge.list
#this 'sampleguide' should be in the format (for all samples):
# $filename	samplename

mkdir -p run_star/mergestarbams
mkdir -p star_align/mergedstarbams

while read line ; do
	#echo $line 
        ID=`echo "${line}" | sed 's/\s.*//g'`
        sampleline=`echo "${line}" | sed 's/.*\s//g'`
	#create array of samples to merge. 
	IFS=',' read -a samples <<<"$sampleline"
	input=()
	for samp in "${samples[@]}" ;do
		input="I=star_align/"${samp}"/alignments.bam "${input}
	done
echo "
#! /bin/bash
#BSUB -J mergcall_"${ID}"
#BSUB -e log/mergecall"${ID}".e
#BSUB -o log/mergecall"${ID}".o
#BSUB -q alloc
#BSUB -W 12:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_PBG
# mem requirements: 
#BSUB -R "rusage[mem=2500]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools
##Merge BAMs
java -jar "${picard}"/MergeSamFiles.jar "${input}" O=star_align/mergedstarbams/"${ID}".bam ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=SILENT
##Dedup
java -jar "${picard}"/MarkDuplicates.jar I=star_align/mergedstarbams/"${ID}".bam O=star_align/mergedbams/"${ID}".md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/mergedbams/"${ID}".MDoutputmetrics
##Realign

" > run_star/mergestarbams/${ID}.lsf
echo $ID 
done < ${FilesToRun}
