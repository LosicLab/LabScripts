FilesToRun=./hb_fastq.list
picard=/hpc/packages/minerva-common/picard/1.93/bin
gatk=/hpc/packages/minerva-common/gatk-mssm/3.2.0/target
#ref=/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/Homo_sapiens.GRCh37.75.dna.primary_assembly.sort.fa
#ref=/sc/orga/scratch/akersn01/tempref/GRCh38_Gencode21/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
ref=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
vcfmills=/sc/orga/work/akersn01/ref/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
vcf1000g=/sc/orga/work/akersn01/ref/gatk/1000G_phase1.indels.hg19.sites.vcf
vcfdbsnp=/sc/orga/work/akersn01/ref/gatk/dbsnp_138.hg19.sort.vcf
sampleguide=./HB_Sample_Codes
#this 'sampleguide' should be in the format (for all samples):
# $filename	samplename

mkdir -p run_star/premerge

while read line ; do
	echo $line 
        R1=`echo "${line}" | sed 's/\s.*//g'`
        R2=`echo "${line}" | sed 's/.*\s//g'`
        filename=$(basename "${R1}") #remove path
        filename="${filename%.*}" #remove extension
        shortname="${filename}" # sometimes I use this line to shorten up filenames.
	shortername=`grep "${shortname}" $sampleguide | cut -f2 |head -n 1 ` 
	echo $shortname $shortername
echo "
#! /bin/bash
#BSUB -J picard_"${shortname}"
#BSUB -e log/pic"${shortname}".e
#BSUB -o log/pic"${shortname}".o
#BSUB -q alloc
#BSUB -W 12:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_PBG
# mem requirements: ~25Gb for markdups. 
#BSUB -R "rusage[mem=2500]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools R
##Most of the below is from the GATK best practices

#remove dups, ~45 mins
java -jar "${picard}"/MarkDuplicates.jar I=star_align/"${shortname}"/alignments.bam O=star_align/${shortname}/md.alignments.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/${shortname}/output.metrics
	#this next bit can be removed when I fix my star coordinates/read-groups.  Takes ~1.5hr
	#fix header, readgroups in bamfiles.
		#samtools view -H star_align/${shortname}/md.alignments.bam | sed -e 's/DUMMYRG/${shortname}\tLB:rna\tPL:ILLUMINA\tPU:machine\tSM:${shortername}/' | samtools reheader - star_align/${shortname}/md.alignments.bam > star_align/${shortname}/rh.md.alignments.bam
		#samtools view -h star_align/${shortname}/rh.md.alignments.bam | sed 's/DUMMYRG/${shortname}/g' |samtools view -bS -  >star_align/${shortname}/rg1.rh.md.alignments.bam
	#replace the readgroups
		#java -jar ${picard}/AddOrReplaceReadGroups.jar I=star_align/${shortname}/rg1.rh.md.alignments.bam O=star_align/${shortname}/rg.rh.md.alignments.bam SO=coordinate RGID=${shortname} RGLB=rna RGPL=illumina RGPU=machine RGSM=${shortername}
	#re-order the bam file. 
		#java -jar ${picard}/ReorderSam.jar I= star_align/${shortname}/rg.rh.md.alignments.bam O= star_align/${shortname}/sort.rg.rh.md.alignments.bam REFERENCE=${ref}
	#index the sorted file. 
		#samtools index star_align/${shortname}/sort.rg.rh.md.alignments.bam 

##Split reads. Use the top line if going through the mess above.  Use the bottom line if star read groups are fixed.  Takes ~1.5 hr
#java -jar ${gatk}/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref} -I star_align/${shortname}/sort.rg.rh.md.alignments.bam -o star_align/${shortname}/split.alignments.bam -U ALLOW_N_CIGAR_READS
java -jar ${gatk}/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref} -I star_align/${shortname}/md.alignments.bam -o star_align/${shortname}/split.alignments.bam -U ALLOW_N_CIGAR_READS

#Indel Realignment
java -jar ${gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} -I star_align/${shortname}/split.alignments.bam -known ${vcfmills} -known ${vcf1000g} -o star_align/${shortname}/target_intervals.list
java -jar ${gatk}/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref} -I star_align/${shortname}/split.alignments.bam -targetIntervals star_align/${shortname}/target_intervals.list -known ${vcfmills} -known ${vcf1000g} -o star_align/${shortname}/realigned_reads.bam

#Base recalibration
java -jar ${gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref} -I star_align/${shortname}/realigned_reads.bam -knownSites ${vcfmills} -knownSites ${vcf1000g} -knownSites ${vcfdbsnp} -o star_align/${shortname}/recalibration_report.grp
java -jar ${gatk}/GenomeAnalysisTK.jar -T PrintReads -R ${ref} -I star_align/${shortname}/realigned_reads.bam -BQSR star_align/${shortname}/recalibration_report.grp -o star_align/${shortname}/recalib.bam

#generate a before table, plot (this is new and untested)
java -jar ${gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${ref} -I star_align/${shortname}/realigned_reads.bam -knownSites ${vcfmills} -knownSites ${vcf1000g} -knownSites ${vcfdbsnp} -BQSR star_align/${shortname}/recalibration_report.grp -o star_align/${shortname}/post_recal_data.table
java -jar ${gatk}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${ref} -before star_align/${shortname}/recalibration_report.grp -after star_align/${shortname}/post_recal_data.table -plots star_align/${shortname}/recalibration_plots.pdf

" > run_star/premerge/${shortname}.lsf

echo $shortname 
done < ${FilesToRun}
