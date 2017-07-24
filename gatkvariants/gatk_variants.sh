#USAGE ./run_bsub_star.sh
####requires a file: fastq.list  in the format:
# /path/to/file1_R1.fastq /path/to/file1_R2.fastq
# /path/to/file2_R1.fastq.gz /path/to/file2_R2.fastq.gz
# fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
# ..etc.


mkdir -p run_star/pass1
mkdir -p run_star/pass2
mkdir -p run_star/gatk_premerge
mkdir -p log
mkdir -p star_align/custom_genome

###CHANGE THESE PARAMS AS NEEDED:
FilesToRun=$1
mergefile=$2 # if you want to merge multiple files, have this file present.  format: filename samplename  ie: Case1  123455_R1,123456_R1
#if no files to merge, just have no file.  
if [[ $mergefile == "" ]] ; then
	mergefile=NOMERGEFILE_DUMMYNAME12345609834203948034895
fi
if [[ -f $mergefile ]]
then
	mkdir -p star_align/mergedbams
	mkdir -p run_star/mergecall
fi

#star parameters
overhang=99 #should be read length - 1.  
Genome=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/2.4.2/99bp
stargtf=/sc/orga/projects/losicb01a/common_folder/ref/hg38/gtf/gencode.v23.primary_assembly.annotation.gtf
star=STAR
module="module load star/2.4.2a"
fasta=/sc/orga/projects/losicb01a/common_folder/ref/hg38/fasta/GRCh38.primary_assembly.genome.gencode.fa 
SJminsubjects=5
SJminreads=5
#gatk parameters
picard=/hpc/packages/minerva-common/picard/1.93/bin
gatk=/hpc/packages/minerva-common/gatk-mssm/3.2.0/target
vcfmills=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.sort.vcf
vcf1000g=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/1000G_phase1.indels.hg38.sites.sort.vcf 
vcfdbsnp=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/dbsnp_138.hg38.sort.vcf

### Create genome generation script:
echo "
#! /bin/bash
#BSUB -J create_custom_genome
#BSUB -e log/mkgenome.e
#BSUB -o log/mkgenome.o
#BSUB -q expressalloc
#BSUB -W 2:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_apollo
#BSUB -R "rusage[mem=6700]" 
#BSUB -R "span[hosts=1]"
$module
cd star_align/
/hpc/users/akersn01/scripts/gatkvariants/SJ_joinfilter.sh $SJminreads $SJminsubjects > custom_genome/SJ.all.out
cd ../
${star} --runMode genomeGenerate --genomeDir star_align/custom_genome --genomeFastaFiles ${fasta} --sjdbGTFfile ${stargtf} --limitGenomeGenerateRAM 90000000000 --sjdbFileChrStartEnd star_align/custom_genome/SJ.all.out --sjdbOverhang ${overhang} --runThreadN 11

" > run_star/custom_genome.lsf


##create run files for each file  
while read full_line ; do
	line=`echo $full_line | sed 's/\s.*//' `
	replicates=`echo $line | sed 's/,/\n/g' |wc -l ` #number of replicates for a given tissue/subject combo
        IFS=',' read -a fastqs <<< "$line" #put fastq files into an array ie ${fastqs[0]}
        #if your samples are named by the file:
	temp="$(basename "${fastqs[0]}")"
	group=`echo $temp |sed 's/\..*//' `
	#if your samples are named by the directory:
	#temp="$(dirname "${fastqs[0]}")"
        #group=`echo $temp |sed 's/.*\///'`
        echo -n $group" "
	if [[ ${fastqs[0]} == *.gz ]] ; then
                decompression='zcat'
        elif [[ ${fastqs[0]} == *.qp ]] ; then
                decompression='/sc/orga/projects/STARNET/Tissues.RNA-seq.raw/apps/quip/src/quip -d -c'
        else
                decompression='cat'
        fi
	seq=`$decompression ${fastqs[0]} |head -n 2 |tail -n 1 `
        readLength=`expr length $seq`
	echo "read length $readLength"
	readGroups=()
        for element in ${fastqs[@]} ; do
                # define read group:
                header=`$decompression $element |head -n 1 `
                IFS=':' read -a sampleinfo <<< "$header"
                #casava 1.8+
                #0 @<instrument> 1:<run number> 2:<flowcell ID> 3:<lane> 4:<tile> 5:<x-pos> 6:<y-pos> <read> 7:<is filtered> 8:<control number> 9:<index sequence>
                #or older:  0 @<instrument> 1:<lane> 2:<tile> 3:<x-pos> 4:<y-pos> 5:index? 6:pairID
                #read Groups using:     ID: LB: PL: PU: SM:
                # ID and PU should both be unique, can be the same.  sample_flowcell.lane
                # LB is the prep library.  unless I have better info, sample is fine (same as SM)
                # PL is Illumina
                machine=`echo ${sampleinfo[0]} | sed 's/@//' `
                #casava 1.8
                rglinestring="ID:"${group}"_"${sampleinfo[2]}.${sampleinfo[3]}" LB:"${group}" PL:ILLUMINA PU:"${group}"_"${sampleinfo[2]}.${sampleinfo[3]}" SM:"${group}" CN:"${machine}
                #older illumina
                #rglinestring="ID:"${group}"."${sampleinfo[1]}" LB:"${group}" PL:ILLUMINA PU:"${group}"."${sampleinfo[1]}" SM:"${group}" CN:"${machine}
                readGroups+=$rglinestring","
        done
        #echo ${group} ${line} >>starnet_ID_fastq.kipp
        #remove the trailing comma from the read group.
	RGline=`echo $readGroups | sed 's/,$//' `
##Begin Pass1 files:
echo "
#! /bin/bash
#BSUB -J star_"${group}"
#BSUB -e log/"${group}".e
#BSUB -o log/"${group}".o
#BSUB -q alloc
#BSUB -W 1:00
#BSUB -n 12
#BSUB -P acc_apollo
#BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
$module
cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir ${Genome} --sjdbGTFfile ${stargtf} --readFilesIn ${full_line} --runThreadN 11 --readFilesCommand ${decompression} --chimSegmentMin 15 \
--chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_Unsorted --outSAMtype BAM Unsorted >alignments.bam
" > run_star/pass1/"${group}".lsf

##Begin Pass2 files:
echo "
#! /bin/bash
#BSUB -J star_"${group}"
#BSUB -e log/"${group}".e
#BSUB -o log/"${group}".o
#BSUB -q alloc
#BSUB -W 1:00
#BSUB -n 12
#BSUB -P acc_apollo
#BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
$module samtools
cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir ../custom_genome/ --sjdbGTFfile ${stargtf} --readFilesIn ${full_line} --outSAMattrRGline $RGline \
--readFilesCommand ${decompression} --outSAMmapqUnique 60 --outSAMattributes NH HI AS nM NM MD --runThreadN 11 --outReadsUnmapped Fastx --chimSegmentMin 15 \
--chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate >alignments.bam
samtools index alignments.bam

" > run_star/pass2/"${group}".lsf

# Run GATK
echo "
#! /bin/bash
#BSUB -J gatkvars_"${group}"
#BSUB -e log/gatkvars_"${group}".e
#BSUB -o log/gatkvars_"${group}".o
#BSUB -q alloc
#BSUB -W 18:00
#BSUB -n 12
#BSUB -P acc_apollo
# mem requirements: ~25Gb for markdups. 
#BSUB -R "rusage[mem=2500]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools R
#rm dups
java -jar "${picard}"/MarkDuplicates.jar I=star_align/"${group}"/alignments.bam O=star_align/${group}/md.alignments.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/${group}/output.metrics && rm star_align/"${group}"/alignments.bam
#split reads
java -jar ${gatk}/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${fasta} -I star_align/${group}/md.alignments.bam -o star_align/${group}/split.alignments.bam -U ALLOW_N_CIGAR_READS && rm star_align/${group}/md.alignments.bam
#indel realignment
java -jar ${gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${fasta} -I star_align/${group}/split.alignments.bam -known ${vcfmills} -known ${vcf1000g} -o star_align/${group}/target_intervals.list
java -jar ${gatk}/GenomeAnalysisTK.jar -T IndelRealigner --maxReadsInMemory 1000000 -R ${fasta} -I star_align/${group}/split.alignments.bam -targetIntervals star_align/${group}/target_intervals.list -known ${vcfmills} -known ${vcf1000g} -o star_align/${group}/realigned_reads.bam && rm star_align/${group}/split.alignments.bam
#Base recalibration
java -jar ${gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${fasta} -I star_align/${group}/realigned_reads.bam -knownSites ${vcfmills} -knownSites ${vcf1000g} -knownSites ${vcfdbsnp} -o star_align/${group}/recalibration_report.grp
java -jar ${gatk}/GenomeAnalysisTK.jar -T PrintReads -R ${fasta} -I star_align/${group}/realigned_reads.bam -BQSR star_align/${group}/recalibration_report.grp -o star_align/${group}/recalib.bam
#generate a before table, plot
java -jar ${gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${fasta} -I star_align/${group}/realigned_reads.bam -knownSites ${vcfmills} -knownSites ${vcf1000g} -knownSites ${vcfdbsnp} -BQSR star_align/${group}/recalibration_report.grp -o star_align/${group}/post_recal_data.table
java -jar ${gatk}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${fasta} -before star_align/${group}/recalibration_report.grp -after star_align/${group}/post_recal_data.table -plots star_align/${group}/recalibration_plots.pdf && rm star_align/${group}/realigned_reads.bam
#final bam output: recalib.bam
" > run_star/gatk_premerge/"${group}".lsf
#if we don't need to merge bams:
if [ ! -f $mergefile ] ; then
echo "
#call variants, filter variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${fasta} -I star_align/${group}/recalib.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o star_align/${group}/variants.vcf
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantFiltration -R ${fasta} -V star_align/${group}/variants.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o star_align/${group}/variants.filter.vcf
" >> run_star/gatk_premerge/"${group}".lsf
fi
done < $FilesToRun

#here we fork.  If mergefile exists, we merge bams, and redo all the de-dupping, indel realignment, etc.  If not, we just call variants (above).
if [ -f $mergefile ]
then
	while read mergeline ; do
        ID=`echo "${mergeline}" | sed 's/\s.*//g'`
	echo $ID
        sampleline=`echo "${mergeline}" | sed 's/.*\s//g'`
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
#BSUB -P acc_apollo
# mem requirements: 
#BSUB -R "rusage[mem=2500]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools
##Merge BAMs
java -jar "${picard}"/MergeSamFiles.jar "${input}" O=star_align/mergedbams/"${ID}".bam ASSUME_SORTED=true USE_THREADING=true VALIDATION_STRINGENCY=SILENT
##Dedup
java -jar "${picard}"/MarkDuplicates.jar I=star_align/mergedbams/"${ID}".bam O=star_align/mergedbams/"${ID}".md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=star_align/mergedbams/"${ID}".MDoutputmetrics && rm star_align/mergedbams/"${ID}".bam
##Realign
java -jar ${gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${fasta} -I star_align/mergedbams/${ID}.md.bam -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.target_intervals.list 
java -jar ${gatk}/GenomeAnalysisTK.jar -T IndelRealigner --maxReadsInMemory 1000000 -R ${fasta} -I star_align/mergedbams/${ID}.md.bam -targetIntervals star_align/mergedbams/${ID}.target_intervals.list -known ${vcfmills} -known ${vcf1000g} -o star_align/mergedbams/${ID}.md.realn.bam && rm star_align/mergedbams/${ID}.md.bam
#call variants, filter variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${fasta} -I star_align/mergedbams/${ID}.md.realn.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o star_align/mergedbams/${ID}.vcf
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantFiltration -R ${fasta} -V star_align/mergedbams/${ID}.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o star_align/mergedbams/${ID}.filter.vcf
#final bam: ID.md.realn.bam
" > run_star/mergecall/${ID}.lsf
done < $mergefile
fi
