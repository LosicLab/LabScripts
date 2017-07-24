#input #1 is a file in the format [FASTQ1][whitespace][FASTQ2]
#input #2 is my tcga summary file

mkdir -p log_mixcr run_mixcr data_mixcr

cat $1 | while read full_line ; do
        line=`echo $full_line | sed 's/\s.*//' ` #gather just R1
	#use the parent Dir + summary file to get a name.  
	temp="$(dirname "${line}")"
	group=`echo $temp |sed 's/.*\///'`
	nom=$group
	#nom=`grep $group $2 | cut -f2` #this is pretty situation specific.  In tcga the handy name is in column 2.  
	echo $group $nom
echo "
#! /bin/bash
#BSUB -J "${nom}"
#BSUB -e log_mixcr/"${nom}".e
#BSUB -o log_mixcr/"${nom}".o
#BSUB -q alloc
#BSUB -W 5:00
#BSUB -n 12
#BSUB -P acc_apollo
#BSUB -R "rusage[mem=4500]" 
#BSUB -R "span[hosts=1]"
module load mixcr/1.8.2

mixcr align --parameters rna-seq -OallowPartialAlignments=true ${full_line} data_mixcr/${nom}_alignments.partial.vdjca
mixcr assemblePartial data_mixcr/${nom}_alignments.partial.vdjca data_mixcr/${nom}_alignmentsRescued.vdjca
#specify index to assemble step to write index (and index_p) file to preserve clone-readID mapping
mixcr assemble data_mixcr/${nom}_alignmentsRescued.vdjca data_mixcr/${nom}_clones.Rescued.clns
mixcr exportClones -s data_mixcr/${nom}_clones.Rescued.clns data_mixcr/${nom}_clones.Rescued.txt

" > run_mixcr/${nom}_mixcr.lsf
done ;
