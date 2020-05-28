%%bash
# Number of SNPs without for each individual
KG_DIR=/datasets/cs284s-sp20-public/1000Genomes
HOME_DIR=~/cse284/datasets/1KGenomes
for chrom in $(seq 1 22)
do
    NUM_SNPS=$(zcat ${KG_DIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | wc -l)
    echo -e "${chrom}\t${NUM_SNPS}" >> ${HOME_DIR}/num_snps.tab
done