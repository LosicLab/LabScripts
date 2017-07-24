#USAGE ./run_bsub_star.sh
#### Wrote this script for STARNET.  Goes through a file with all the fastqs for each subject/tissue combo and does one STAR run for each.
### modified here for tcga. 
#each different fastq file is separated by a comma. 

mkdir -p run_star
mkdir -p log
mkdir -p star_align

###CHANGE THESE PARAMS AS NEEDED:
#FilesToRun=/sc/orga/projects/PBG/TCGA/RNA/BRCA/manifests/fastq.gz.pairs  #format : fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
FilesToRun=$1  #format : fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
star=/hpc/users/akersn01/software/STAR/bin/Linux_x86_64/STAR
####

while read full_line ; do
        line=`echo $full_line | sed 's/\s.*//' `
	#line=$full_line
        #group=`echo $full_line | sed 's/\s.*//' `
	
	## Figure out if single or paired end ##
	pairs=`echo $full_line | sed 's/\s/\n/g' |wc -l ` #number of replicates for a given tissue/subject combo
	if [[ ${pairs} == 1 ]] ; then trinityline='--single Unmapped.out.mate1'
	elif [[ ${pairs} == 2 ]] ; then trinityline='--left Unmapped.out.mate1 --right Unmapped.out.mate2'
	else
		echo "Error, wrong number of pairs (${pairs}) for ${full_line}, trinity will fail"; 
	fi
	IFS=',' read -a fastqs <<< "$line" #put fastq files into an array ie ${fastqs[0]}
	###Group/individual identier:
	##Option 1:
	group="$(basename "${fastqs[0]}")"
	##Option 2:
	#temp="$(dirname "${fastqs[0]}")"
	#group=`echo $temp |sed 's/.*\///'`
	echo $group
	if [[ ${fastqs[0]} == *.gz ]] ; then
		decompression='zcat'
	elif [[ ${fastqs[0]} == *.qp ]] ; then
		decompression='/sc/orga/projects/STARNET/Tissues.RNA-seq.raw/apps/quip/src/quip -d -c'
	else
		decompression='cat'
	fi
	#get the read length (I know that all read lengths are the same within a given tissue/subject combo, with 16 exceptions.  These will be aligned to a 50bp overhang)
	seq=`$decompression ${fastqs[0]} |head -n 2 |tail -n 1 `
        readLength=`expr length $seq`
	if (( $readLength < 60 )) ; then
		Genome=/sc/orga/projects/losicb01a/common_folder/ref/hg19/star/overhang49
	else
		Genome=/sc/orga/projects/losicb01a/common_folder/ref/hg19/star/overhang100
	fi
	#get the tissue
	#IFS='/' read -a path <<< "${fastqs[0]}"
	#pathlength=${#path[@]}
	#tissue=${path[(${pathlength}-2)]} #assigns the parent directory value to tissue eg VAF or AOR
	#pull read groups from fastq files
	readGroups=()
	for element in ${fastqs[@]} ; do
		# define read group:
		header=`$decompression $element |head -n 1 `
		IFS=':' read -a sampleinfo <<< "$header"
		#casava 1.8+
		#0 @<instrument> 1:<run number> 2:<flowcell ID> 3:<lane> 4:<tile> 5:<x-pos> 6:<y-pos> <read> 7:<is filtered> 8:<control number> 9:<index sequence>
	#or older:  0 @<instrument> 1:<lane> 2:<tile> 3:<x-pos> 4:<y-pos> 5:index? 6:pairID
		#read Groups using:	ID: LB: PL: PU: SM:
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
##Begin LSF files:
echo "
#! /bin/bash
#BSUB -J star_"${group}"
#BSUB -e log/"${group}".e
#BSUB -o log/"${group}".o
#BSUB -q alloc
#BSUB -W 2:30
#BSUB -n 12
#BSUB -P acc_PBG
#BSUB -R \"rusage[mem=3500]\" 
#BSUB -R \"span[hosts=1]\"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
#### manda can be changed to mothra for faster runs but generally will queue slower
module load samtools bedtools
module load python py_packages
module load trinity

cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir $Genome --readFilesIn $full_line --outSAMattrRGline $RGline --outSAMmapqUnique 60 --outSAMattributes NH HI AS nM NM MD --runThreadN 11 --outReadsUnmapped Fastx \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --readFilesCommand ${decompression} --outSAMtype BAM Unsorted

#--chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --outStd BAM_SortedByCoordinate --readFilesCommand ${decompression} --outSAMtype BAM Unsorted SortedByCoordinate >alignments.bam
#samtools index alignments.bam

Trinity --seqType fq --JM 35G --CPU 11 --no_run_chrysalis --no_run_butterfly ${trinityline} --output trinity_inchworm

python /hpc/users/akersn01/scripts/sj2psi_x.py SJ.out.tab 5 10

" > run_star/"${group}".lsf""
done < $FilesToRun
