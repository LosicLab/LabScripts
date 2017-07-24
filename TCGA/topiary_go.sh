# topiary_go.sh #FilesList
# file list should be a two column file with MAF file in the first column and an HLA column in the second column.

#directory containing your peptides
peptides=/sc/orga/projects/immunotherapy/peptide_prediction/human/classI/alllengths/
class=classI
#class=classII
lengths=(8 9 10 11 12)
predictor=netmhccons  #other options
if [ $class = "classII" ]; then
	predictor=netmhciipan
	peptides=/sc/orga/projects/immunotherapy/peptide_prediction/human/classII/lt500/
	lengths=(15)
fi
#setup
mkdir -p run/${class}
mkdir -p log/${class}
mkdir -p data/${class}

while read line ; do
	MAF=`echo $line |sed 's/\s.*//'`
	MHC=`echo $line |sed 's/.*\s//'`
	ID=`echo $MAF |sed 's/.*\///' |sed 's/.maf.txt//'`
	echo $ID
	for peplength in ${lengths[@]} ; do 
echo "
#! /bin/bash
#BSUB -J Topiary_"${ID}"_"${peplength}"
#BSUB -e log/"${class}"/"${ID}"_"${peplength}".e
#BSUB -o log/"${class}"/"${ID}"_"${peplength}".o
#BSUB -q alloc
#BSUB -W 2:00
#BSUB -n 3
#BSUB -P acc_PBG
#BSUB -R "rusage[mem=12500]" 
#BSUB -R "span[hosts=1]"
module load python py_packages R
PYTHONPATH=/hpc/users/akersn01/.local/lib/python2.7/site-packages/:\${PYTHONPATH}
PATH=/hpc/users/akersn01/.local/bin/:/hpc/users/akersn01/bin/:\${PATH}
topiary \
--mhc-predictor "${predictor}" \
--mhc-alleles-file "${MHC}" \
--vcf "${MAF}" \
--mhc-epitope-lengths "${peplength}" \
--wildtype-ligandome-directory "${peptides}" \
--ic50-cutoff 500 \
--percentile-cutoff 100 \
--skip-variant-errors \
--output-csv data/"${class}"/"${ID}"_"${peplength}".csv

" > run/${class}/${ID}_${peplength}.lsf
done 
done < $1
