"""
Adam Klie
05/27/2020
Convert an openSNP DataFrame with allele characters as genotypes (e.g., CC, AT, etc.)
to a similar DataFrame with numbers representing the number of minor alleles from 1000Genomes (e.g. 0, 1, 2)
"""

import os
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
tqdm.pandas()


# Function to convert genotype from nucleotides to numbers 
# 0 (homozygous reference), 1(heteroozygous), 2 (homozygous minor allele)
def get_gt(x):
    for col in x.index[3:-2]:       
        if isinstance(x[col], str):
            ref_allele = x["REF"]
            gt = 2 - x[col].upper().count(ref_allele)
            x[col] = gt
    return x


#Main function
def main(args):
    
    # Open log file for writing
    log = open("{}/alleleToNum.log".format(args.log_dir), 'w')
    
    # Read in filtered genotypes
    filtered_genotypes = pd.read_csv(args.genotype_file, sep='\t')
    log.writelines("Number of SNPs in {}: {}\n".format(args.genotype_file, len(filtered_genotypes)))
    
    # Replace -- with NaN
    filtered_genotypes = filtered_genotypes.replace('--', np.nan)

    # Read in the final SNP the final SNP set
    oneK_ref = pd.read_csv(args.rsids, sep='\t')
    log.writelines("Number of SNPs in {}: {}\n".format(args.rsids, len(oneK_ref)))
    
    # Keep only SNPs found in 1000Genomes
    gt_char = filtered_genotypes.merge(oneK_ref, left_on="rsid", right_on="ID")

    # WARNING: Getting the numbers takes about an hour and half
    gt_nums = gt_char.progress_apply(get_gt, axis=1)
    gt_nums = gt_nums.drop("ID", axis=1)
    
    # Save a final SNP vs indiviudal table with other SNP info in it
    log.writelines("Number of SNPs in output: {}\n".format(len(gt_nums)))
    gt_nums.to_csv("{}/openSNP_final_genotypes.tsv".format(args.output_dir), sep='\t', index=False)
    
               
if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_filtered_genotypes.tsv")
    ONEK_RSIDS = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "oneK_rsids.tsv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to openSNP filtered genotype tsv, default is {}'.format(GENOTYPES))
    parser.add_argument('-id', '--rsids', type=str, default=ONEK_RSIDS, help='filepath to 1000Genomes rsids with ref allele default is {}'.format(ONEK_RSIDS))
    parser.add_argument('-p', '--pickle', type=bool, default=False, help='Boolean for whether or not to save output as a pickle file')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)