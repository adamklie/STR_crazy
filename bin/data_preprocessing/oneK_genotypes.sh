"""
Adam Klie
05/22/2020
Code to extract all SNPs found in openSNP from 1000Genomes vcf
"""

KG_DIR=/datasets/cs284s-sp20-public/1000Genomes  # Path to 1000Genomes vcf (per chromosome)
SNP_LIST=~/project/datasets/openSNP/data/openSNP_filtered_rsids.txt  # Path to SNP list (one per line from openSNP)
DATA=~/project/datasets/oneKGenomes/data  # Path to directory to house output

for chrom in $(seq 1 22)
do
    echo -e "Adding SNPs from chromosome $chrom"
    vcftools \
        --snps ${SNP_LIST} \
        --gzvcf ${KG_DIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --recode \
        --out ${chrom}
    bgzip ${chrom}.recode.vcf
done

bcftools concat *recode.vcf.gz > oneK_genotypes.vcf
rm *recode.vcf.gz
bgzip oneK_genotypes.vcf
tabix -p vcf oneK_genotypes.vcf.gz

bcftools query -H -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" oneK_genotypes.vcf.gz | \
    sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' |  \
    grep -v "|" | sed -r 's/\[.{1,4}\]//g' | sed 's/:GT//g' \
    > oneK_genotypes.tsv
    
mv oneK_genotypes.tsv $DATA
mv oneK_genotypes.vcf* $DATA
