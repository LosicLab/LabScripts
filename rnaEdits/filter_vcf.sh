FilesToRun=$1 #one vcf per line, full path. 
picard=/hpc/packages/minerva-common/picard/1.93/bin
gatk=/hpc/packages/minerva-common/gatk-mssm/3.2.0/target
#ref=/sc/orga/work/akersn01/ref/star/ENSEMBL.homo_sapiens.release-75_overhang100/Homo_sapiens.GRCh37.75.dna.primary_assembly.sort.fa
#ref=/sc/orga/projects/losicb01a/common_folder/ref/hg38/fasta/hg38.fa
ref=/sc/orga/projects/losicb01a/common_folder/ref/hg38/star/GRCh38_Gencode21_overhang100/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
alubed=/sc/orga/work/akersn01/ref/Alu_elements_hg38.bed
custom_genome=/sc/orga/scratch/akersn01/hb/gatk/star_align/custom_genome
mkdir -p run log data

while read line ; do
	ID=`echo $line |sed 's/.*\///' |sed 's/\..*//'`
echo "
#! /bin/bash
#BSUB -J filter_"${ID}"
#BSUB -e log/"${ID}".e
#BSUB -o log/"${ID}".o
#BSUB -q low
#BSUB -W 4:00
#BSUB -n 4
#BSUB -m manda
#BSUB -P acc_PBG
# mem requirements: 
#BSUB -R "rusage[mem=10000]" 
#BSUB -R "span[hosts=1]"
module purge
module load picard/1.93 gatk-mssm/3.2.0 samtools bedtools vcftools/0.1.12b
#annotate variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantAnnotator -R ${ref} -V $line \
-A HomopolymerRun -A StrandOddsRatio -A TandemRepeatAnnotator -o data/${ID}.anno.vcf

#filter variants
java -jar ${gatk}/GenomeAnalysisTK.jar -T VariantFiltration -R ${ref} -V data/${ID}.anno.vcf \
-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -filterName MQ -filter \"MQ >= 25 \" \
-o data/${ID}.filter.vcf
~/scripts/gatkvariants/rnaEdits/vcf_parser_rnaedits.pl data/${ID}.filter.vcf > data/${ID}.filter2.vcf

#separate into Alu/nonAlu
bedtools intersect -header -wa -v -a data/${ID}.filter2.vcf -b ${alubed} > data/${ID}.NotAlu.vcf
bedtools intersect -header -wa -u -a data/${ID}.filter2.vcf -b ${alubed} > data/${ID}.Alu.vcf

#further filter nonAlu
grep -v STR data/${ID}.NotAlu.vcf | ~/scripts/gatkvariants/rnaEdits/hrun_filter.pl > data/${ID}.NotAlu.filter1.vcf
	###How to generate the splice junctions db:
	# awk '{ print $1,($2-5),($2+5)"\n"$1,($3-5),($3+5) }' OFS="\t" sjdbList.out.tab > splice_junctions_abparts.bed
	# cat Abparts.bed >> splice_junctions_abparts.bed
bedtools intersect -header -wa -v -a data/${ID}.NotAlu.filter1.vcf -b ${custom_genome}/splice_junctions_abparts.bed > data/${ID}.NotAlu.filtered.vcf

##remove indels:
vcftools --vcf data/${ID}.Alu.vcf --remove-indels --recode --recode-INFO-all --out data/${ID}.Alu.noindel
vcftools --vcf data/${ID}.NotAlu.filtered.vcf --remove-indels --recode --recode-INFO-all --out data/${ID}.NotAlu.filtered.nodindel

cut -f1-5 data/${ID}.NotAlu.filtered.nodindel.recode.vcf |grep -v -P "^#" >>data/all.NotAlu.sites
cut -f1-5 data/${ID}.Alu.noindel.recode.vcf |grep -v -P "^#" >>data/all.Alu.sites

" > run/${ID}.lsf
echo $ID 
done < ${FilesToRun}
