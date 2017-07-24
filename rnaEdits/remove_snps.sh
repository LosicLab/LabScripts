myvcf=$1
vcfdbsnp=/sc/orga/projects/losicb01a/common_folder/ref/hg38/vcf/dbsnp_138.hg38.sort.vcf

~/software/vcflib/vcfintersect -v -w 1 \
--intersect-vcf $vcfdbsnp \
-r /sc/orga/projects/losicb01a/common_folder/ref/hg38/fasta/hg38.fa $myvcf > ${myvcf}.nosnps.vcf


