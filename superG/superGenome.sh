##The goal of this pipeline is to create a more complete assignment of reads, with less dependence on a reference.
##I'm writing this for a batch of sequencing files. 

FilesToRun=$1	#format : fastq1_R1,fastq2_R1 fastq1_R2,fastq2_R2
IDSource="ParentDir" ##alternatively use "FileName" if unique IDs come from the Filename.  


#important pathways:
SCRIPTPATH=$( cd "$(dirname "$0")" ; pwd -P )
starload="star/2.4.2a"
trinityload="trinity/2014-04-23"
blastload="blast/2.2.26+"
subreadload="subread/1.4.4"
cdhitload="cdhit/3.1"
#check input fastqs:
#STEP 1: STAR+Trinity
	#this will create submit scripts that run star and trinity.
	${SCRIPTPATH}/star_complete_pass1.sh $FilesToRun $IDsource $starload $trinityload
	##concatenate the trinity output files
	#Option 1: Concat all unmapped reads, run Trinity---think it will be prohibitively slow, depending on cohort size. 
	#Option 2: run trinity on everything, merge results (via blast?)
	find star_align -name "inchworm.150bp.1000reads.fa" -exec cat {} \; > star_align/trinity/all.150bp.1000reads.fa
	#Reduce these hits to representative sequences.
	module load $cdhitload
	cd-hit-est -n 8 -c 0.90 -g 1 -r 1 -i star_align/trinity/all.150bp.1000reads.fa -o star_align/trinity/reps.150bp.1000reads.fa 1> log/cdest.o 2> log/cdest.e

#STEP 2: BLAST
##BLAST my contigs 
	#first blast against human, then bacteria/virus? or just nr 
	${SCRIPTPATH}/blast_trinity_out.sh

#STEP 3: Create a new STAR index
#submit a script to minerva
${SCRIPTPATH}/star_genome.sh $

#wait for it to finish:
while [ ! -f /tmp/list.txt ]
do
  sleep 2
done

#STEP 4: Align with STAR

#STEP 5: Featurecounts
