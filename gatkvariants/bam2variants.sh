#USAGE ./run_bsub_star.sh
####requires a file: bams.list  in the format:
#bam1
#bam2
# ..etc.
##### UNTESTED!! ###
mkdir -p data
mkdir -p log
mkdir -p run

###CHANGE THESE PARAMS AS NEEDED:
FilesToRun=$1

#star parameters
overhang=99 #should be read length - 1.  
Genome=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/
stargtf=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/gencode.v21.annotation.gtf
star=/hpc/users/akersn01/software/STAR/bin/Linux_x86_64/STAR
fasta=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa 
SJminsubjects=5
SJminreads=5
#gatk parameters
picard=/hpc/packages/minerva-common/picard/1.93/bin
gatk=/hpc/packages/minerva-common/gatk-mssm/3.2.0/target
vcfmills=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.sort.vcf
vcf1000g=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/1000G_phase1.indels.hg38.sites.sort.vcf 
vcfdbsnp=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/dbsnp_138.hg38.sort.vcf


##create run files for each file  
while read full_line ; do
	line=`echo $full_line | sed 's/\s.*//' `
        #if your samples are named by the file:
	temp="$(basename "${fastqs[0]}")"
	ID=`echo $temp |sed 's/\..*//' `
	#if your samples are named by the directory:
	#temp="$(dirname "${fastqs[0]}")"
        #group=`echo $temp |sed 's/.*\///'`
        echo -n $group" "
# Run GATK
echo "
#! /bin/bash
#BSUB -J bam2var_"${ID}"
#BSUB -e log/bam2var_"${ID}".e
#BSUB -o log/bam2var_"${ID}".o
#BSUB -q low
#BSUB -W 12:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_PBG
# mem requirements: 
#BSUB -R "rusage[mem=2500]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools
##Dedup
java -jar "${picard}"/MarkDuplicates.jar I=${full_line} O=data/"${ID}".md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/mergedbams/"${ID}".MDoutputmetrics && rm star_align/mergedbams/"${ID}".bam
##Realign
java -jar ${gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${fasta} -I star_align/mergedbams/${ID}.md.bam -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.target_intervals.list 
java -jar ${gatk}/GenomeAnalysisTK.jar -T IndelRealigner -R ${fasta} -I star_align/mergedbams/${ID}.md.bam -targetIntervals star_align/mergedbams/${ID}.target_intervals.list -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.md.realn.bam && rm star_align/mergedbams/${ID}.md.bam
#call variants, filter variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${fasta} -I star_align/mergedbams/${ID}.md.realn.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o star_align/mergedbams/${ID}.vcf
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantFiltration -R ${fasta} -V star_align/mergedbams/${ID}.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o star_align/mergedbams/${ID}.filter.vcf
#final bam: ID.md.realn.bam
" > run_star/mergecall/${ID}.lsf
done < $FilesToRun
