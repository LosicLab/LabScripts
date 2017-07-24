mkdir -p run/star1
mkdir -p log
mkdir -p star_align/trinity


###CHANGE THESE PARAMS AS NEEDED:
FilesToRun=$1  #format : fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
IDSource=$2

## New STAR > 2.4.1
star="STAR"
Genome50=/sc/orga/projects/losicb01a/common_folder/ref/hg19/star/2.4.2/49bp
Genome100=/sc/orga/projects/losicb01a/common_folder/ref/hg19/star/2.4.2/99bp
module="module load "$3" "$4" "samtools
GTF=/sc/orga/projects/losicb01a/common_folder/ref/hg19/gtf/gencode.v19.annotation.gtf

###
while read full_line ; do
        line=`echo $full_line | sed 's/\s.*//' `
	replicates=`echo $line | sed 's/,/\n/g' |wc -l ` #number of replicates for a given tissue/subject combo
	IFS=',' read -a fastqs <<< "$line" #put fastq files into an array ie ${fastqs[0]}
	## Figure out if single or paired end ##
        pairs=`echo $full_line | sed 's/\s/\n/g' |wc -l ` #number of replicates for a given tissue/subject combo
        if [[ ${pairs} == 1 ]] ; then trinityline='--single Unmapped.out.mate1'
        elif [[ ${pairs} == 2 ]] ; then trinityline='--left Unmapped.out.mate1 --right Unmapped.out.mate2'
        else
                echo "Error, wrong number of pairs (${pairs}) for ${full_line}, trinity will fail"; 
        fi
	#group/individual identier:
	if [[ ${IDSource} == 'FileName' ]] ; then
		temp="$(basename "${fastqs[0]}")"
		group=`echo $temp |sed 's/\.fastq.gz//'`
	else
		temp="$(dirname "${fastqs[0]}")"
		group=`echo $temp |sed 's/.*\///'`
	fi
	echo $group
	#Determine Decompression
	if [[ ${fastqs[0]} == *.gz ]] ; then
		decompression='zcat'
	elif [[ ${fastqs[0]} == *.qp ]] ; then
		decompression='/sc/orga/projects/STARNET/Tissues.RNA-seq.raw/apps/quip/src/quip -d -c'
	elif [[ ${fastqs[0]} == *.bz2 ]] ; then
		decompression='bzcat'
	else
		decompression='cat'
	fi
	#get the read length.  Simplifying assumption to use the first read length uncountered.
	seq=`$decompression ${fastqs[0]} |head -n 2 |tail -n 1 `
        readLength=`expr length $seq`
	if (( $readLength < 60 )) ; then
		Genome=$Genome50
	else
		Genome=$Genome100
	fi
	#pull read groups from fastq files
	readGroups=()
	for element in ${fastqs[@]} ; do
		# define read group:
		header=`$decompression $element |head -n 1 `
		IFS=':' read -a sampleinfo <<< "$header"
		#casava 1.8+
			#0 @<instrument> 1:<run number> 2:<flowcell ID> 3:<lane> 4:<tile> 5:<x-pos> 6:<y-pos> <read> 7:<is filtered> 8:<control number> 9:<index sequence>
		#or older:
			# 0 @<instrument> 1:<lane> 2:<tile> 3:<x-pos> 4:<y-pos> 5:index? 6:pairID
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
	#remove the trailing comma from the read group.
	RGline=`echo $readGroups | sed 's/,$//' `
##Begin LSF files:
echo "
#! /bin/bash
#BSUB -J star_"${group}"
#BSUB -e log/"${group}".e
#BSUB -o log/"${group}".o
#BSUB -q low
#BSUB -W 2:00
#BSUB -n 12
#BSUB -P acc_PBG
#BSUB -R \"rusage[mem=3500]\" 
#BSUB -R \"span[hosts=1]\"
##### STAR needs >32GB RAM, ~40 to be safe, and spanned across just one host.
#### this sets us up with 12*3.5=42GB
${module}
cd star_align
mkdir -p ${group}
cd ${group}
${star} --genomeDir $Genome --readFilesIn $full_line --outSAMattrRGline $RGline --outSAMmapqUnique 60 --outSAMattributes NH HI AS nM NM MD --runThreadN 11 --outReadsUnmapped Fastx \
--quantMode GeneCounts \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --readFilesCommand ${decompression} --outSAMtype BAM Unsorted

Trinity --seqType fq --JM 35G --CPU 11 --no_run_chrysalis --no_run_butterfly ${trinityline} --output trinity_inchworm
	#pulling out fasta reads that are >= 150bp and have >=1000 reads of support.  
	fgrep '>' trinity_inchworm/inchworm.K25.L25.DS.fa | sed 's/;/\t/g' |awk '{ if ($2 > 1000 && $10 >=150 ) print $1,$2 }' OFS="\t" |sed 's/\t/;/g' |sed 's/^>//' > inchworm.150bp.1000reads.ids
	samtools faidx trinity_inchworm/inchworm.K25.L25.DS.fa 
	while read ids ; do samtools faidx trinity_inchworm/inchworm.K25.L25.DS.fa $ids ; done < inchworm.150bp.1000reads.ids > inchworm.150bp.1000reads.fa
	

	


" > run/star1/"${group}".lsf""
done < $FilesToRun
