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

mkdir -p run_star/mergecall
mkdir -p star_align/mergedbams

while read line ; do
	#echo $line 
        ID=`echo "${line}" | sed 's/\s.*//g'`
        sampleline=`echo "${line}" | sed 's/.*\s//g'`
	#create array of samples to merge. 
	IFS=',' read -a samples <<<"$sampleline"
	input=()
	for samp in "${samples[@]}" ;do
		input="I=star_align/"${samp}"/recalib.bam "${input}
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
java -jar "${picard}"/MergeSamFiles.jar "${input}" O=star_align/mergedbams/"${ID}".bam ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=SILENT
##Dedup
java -jar "${picard}"/MarkDuplicates.jar I=star_align/mergedbams/"${ID}".bam O=star_align/mergedbams/"${ID}".md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/mergedbams/"${ID}".MDoutputmetrics
##Realign
java -jar ${gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} -I star_align/mergedbams/${ID}.md.bam -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.target_intervals.list
java -jar ${gatk}/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref} -I star_align/mergedbams/${ID}.md.bam -targetIntervals star_align/mergedbams/${ID}.target_intervals.list -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.md.realn.bam
#call variants, filter variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${ref} -I star_align/mergedbams/${ID}.md.realn.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o star_align/mergedbams/${ID}.vcf
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V star_align/mergedbams/${ID}.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o star_align/mergedbams/${ID}.filter.vcf

" > run_star/mergecall/${ID}.lsf
echo $ID 
done < ${FilesToRun}
