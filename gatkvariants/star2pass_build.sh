#USAGE ./run_bsub_star.sh
####requires a file: fastq.list  in the format:
# /path/to/file1_R1.fastq /path/to/file1_R2.fastq
# /path/to/file2_R1.fastq /path/to/file2_R2.fastq
# fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
# ..etc.


mkdir -p run_star/pass1
mkdir -p run_star/pass2
mkdir -p log
mkdir -p star_align/custom_genome
mkdir -p counts

###CHANGE THESE PARAMS AS NEEDED:
FilesToRun=$1
overhang=99 #should be read length - 1.  
fasta=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa 
Genome=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/
stargtf=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/gencode.v21.annotation.gtf
star=/hpc/users/akersn01/software/STAR/bin/Linux_x86_64/STAR
####

### Create genome generation script:
echo "
#! /bin/bash
#BSUB -J create_custom_genome
#BSUB -e log/mkgenome.e
#BSUB -o log/mkgenome.o
#BSUB -q alloc
#BSUB -W 7:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=9000]" 
#BSUB -R "span[hosts=1]"

${star} --runMode genomeGenerate --genomeDir star_align/custom_genome --genomeFastaFiles ${fasta} --sjdbGTFfile ${stargtf} --limitGenomeGenerateRAM 90000000000 --sjdbFileChrStartEnd star_align/custom_genome/SJ.all.out --sjdbOverhang ${overhang} --runThreadN 24

" > run_star/custom_genome.lsf

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
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir ${Genome} --sjdbGTFfile ${stargtf} --readFilesIn ${full_line} --runThreadN 11 --readFilesCommand ${decompression} --chimSegmentMin 15 \
--chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_Unsorted --outSAMtype BAM Unsorted >alignments.bam
" > run_star/pass1/"${group}".lsf""

##Begin Pass2 files:
echo "
#! /bin/bash
#BSUB -J star_"${group}"
#BSUB -e log/"${group}".e
#BSUB -o log/"${group}".o
#BSUB -q alloc
#BSUB -W 1:00
#BSUB -n 12
#BSUB -m manda
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
module load samtools
cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir ../custom_genome/ --sjdbGTFfile ${stargtf} --readFilesIn ${full_line} --outSAMattrRGline $RGline \
--readFilesCommand ${decompression} --outSAMmapqUnique 60 --outSAMattributes NH HI AS nM NM MD --runThreadN 11 --outReadsUnmapped Fastx --chimSegmentMin 15 \
--chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate >alignments.bam
samtools index alignments.bam

" > run_star/pass2/"${group}".lsf""

done < $FilesToRun
