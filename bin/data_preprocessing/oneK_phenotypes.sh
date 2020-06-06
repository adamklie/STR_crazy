"""
Adam Klie
05/22/2020
Code to extract an IrisPlex SNPs from 1000Genomes vcf files
"""


KG_DIR=/datasets/cs284s-sp20-public/1000Genomes  # Path to 1000Genomes vcf (per chromosome)
SNP_LIST=~/project/datasets/iris_rsids.txt  # Path to IrisPlex SNP list (one per line)
DATA=~/project/datasets/oneKGenomes/data  # Path to directory to house output

declare -a arr=(15 14 11 6 5)

## now loop through the above array
for chrom in "${arr[@]}"
do
    echo -e "Adding iris SNPs from chromosome $chrom"
    vcftools \
        --snps ${SNP_LIST} \
        --gzvcf ${KG_DIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --recode \
        --out ${chrom}.iris
    bgzip ${chrom}.iris.recode.vcf
done

bcftools concat *iris.recode.vcf.gz > iris_oneK_genotypes.vcf
rm *iris.recode.vcf.gz
bgzip iris_oneK_genotypes.vcf
tabix -p vcf iris_oneK_genotypes.vcf.gz

bcftools query -H -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n" iris_oneK_genotypes.vcf.gz | \
    sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' |  \
    grep -v "|" | sed -r 's/\[.{1,4}\]//g' | sed 's/:GT//g' \
    > iris_oneK_genotypes.tsv
    
mv iris_oneK_genotypes.tsv $DATA
mv iris_oneK_genotypes.vcf* $DATA