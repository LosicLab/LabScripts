#USAGE ./run_bsub_star.sh
####requires a file: fastq.list  in the format:
# /path/to/file1_R1.fastq /path/to/file1_R2.fastq
# /path/to/file2_R1.fastq /path/to/file2_R2.fastq
# ..etc.


mkdir -p run_star/pass1
mkdir -p run_star/pass2
mkdir -p log
mkdir -p star_align/custom_genome
mkdir -p counts

###CHANGE THESE PARAMS AS NEEDED:
FilesToRun=./hb_fastq.list
sampleguide=./HB_Sample_Codes
readgroups=./readgroups
overhang=99 #should be read length - 1.  
###trying out this sorted fasta for genome generation, hopefully will save time later in the process.  
#fasta=/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/Homo_sapiens.GRCh37.75.dna.primary_assembly.sort.fa
fasta=/sc/orga/scratch/akersn01/tempref/GRCh38_Gencode21/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
#Genome=/sc/orga/work/akersn01/ref/star/star/overhang49/
#Genome=/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/
Genome=/sc/orga/scratch/akersn01/tempref/GRCh38_Gencode21/
stargtf=/sc/orga/scratch/akersn01/tempref/GRCh38_Gencode21/gencode.v21.annotation.gtf
#gtf=/sc/orga/scratch/akersn01/ensembl/Homo_sapiens.GRCh37.75.gtf
#gtf=/sc/orga/projects/losicb01a/reference_seq/ensembl/Homo_sapiens.GRCh37.74.gtf
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
###BSUB -m manda
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=9000]" 
#BSUB -R "span[hosts=1]"


${star} --runMode genomeGenerate --genomeDir star_align/custom_genome --genomeFastaFiles ${fasta} --sjdbGTFfile ${stargtf} --limitGenomeGenerateRAM 90000000000 --sjdbFileChrStartEnd star_align/custom_genome/SJ.all.out --sjdbOverhang ${overhang} --runThreadN 11

" > run_star/custom_genome.lsf





while read line ; do
	R1=`echo "${line}" | sed 's/\s.*//g'`  
	R2=`echo "${line}" | sed 's/.*\s//g'`  
	filename=$(basename "${R1}") #remove path
	filename="${filename%.*}" #remove extension
	shortname=${filename} # sometimes I use this line to shorten up filenames.
	sample=`grep "${shortname}" $sampleguide | cut -f2 |head -n 1 `
	flowcell=`grep ${shortname} ${readgroups} |cut -f4`
	lane=`grep ${shortname} ${readgroups} |cut -f5 `
	rgID=${sample}_${flowcell}.${lane}
	rgPL=ILLUMINA
	rgLB=${sample}
	rgPU=${sample}_${flowcell}.${lane} #PU should be flowcell.lane.barcode, but i don't have barcodes so sample names will suffice here.

##Begin Pass1 files:
echo "
#! /bin/bash
#BSUB -J star_"${shortname}"
#BSUB -e log/"${shortname}".e
#BSUB -o log/"${shortname}".o
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
cd star_align
mkdir -p ${shortname}
cd ${shortname}
${star} --genomeDir ${Genome} --sjdbGTFfile ${stargtf} --readFilesIn ${R1} ${R2} --runThreadN 11 chimSegmentMin --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_Unsorted --outSAMtype BAM Unsorted >alignments.bam
" > run_star/pass1/"${shortname}".lsf""

##Begin Pass2 files:
echo "
#! /bin/bash
#BSUB -J star_"${shortname}"
#BSUB -e log/"${shortname}".e
#BSUB -o log/"${shortname}".o
#BSUB -q alloc
#BSUB -W 1:00
#BSUB -n 12
###BSUB -m manda
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=3500]" 
#BSUB -R "span[hosts=1]"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
module load samtools
cd star_align
mkdir -p ${shortname}
cd ${shortname}
${star} --genomeDir ../custom_genome/ --readFilesIn ${R1} ${R2} --outSAMattrRGline ID:${rgID} LB:${rgLB} PL:${rgPL} PU:${rgPU} SM:${sample} --outSAMmapqUnique 60 --outSAMattributes NH HI AS nM NM MD --runThreadN 11 outReadsUnmapped Fastx --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate >alignments.bam
samtools index alignments.bam

" > run_star/pass2/"${shortname}".lsf""

done < $FilesToRun
