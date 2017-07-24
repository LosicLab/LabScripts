mkdir -p postlog 
mkdir -p postrun

#class="classI"
class="classII"
myfile="$(basename "${1}")"
newfile=`echo $1 |sed 's/csv/post.csv/' `
echo "$myfile"
echo "
#! /bin/bash
#BSUB -J "${myfile}"
#BSUB -e postlog/"${myfile}".e
#BSUB -o postlog/"${myfile}".o
#BSUB -q alloc
#BSUB -W 30:00
#BSUB -n 1
#BSUB -P acc_apollo
#BSUB -R "rusage[mem=5000]" 
##BSUB -R "span[hosts=1]"
module load samtools python py_packages

tail -n +2 $1 | tr \",\" \"\\t\" | sort -k2,2 -k4,4 | /hpc/users/akersn01/scripts/self-ligandome-compare.pl $class | /hpc/users/akersn01/scripts/hla/netmhc/topiary_original_peptide.pl $class |tr \"\\t\" \",\" > ${newfile}.temp
head -n 1 $1 |sed 's/$/,LevenshteinDist,NearestSelfEpitope,NearestSelfIC50,WTEpitope,WTEpIC50,WTEpPercentile/' > ${newfile}.header
cat ${newfile}.header ${newfile}.temp > $newfile
rm ${newfile}.header ${newfile}.temp

" > postrun/${myfile}.lsf 
