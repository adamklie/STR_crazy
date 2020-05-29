"""
Adam Klie
05/27/2020

Input: 
    - 
Output: 
    - 
    - 
Notes:
    - 
"""
import os
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
tqdm.pandas()


# Function to convert genotype from nucleotides to numbers 0 (homozygous reference), 1(heteroozygous), 2 (homozygous minor allele)
def get_gt(x):
    for col in x.index[3:-2]:       
        if isinstance(x[col], str):
            ref_allele = x["REF"]
            gt = 2 - x[col].upper().count(ref_allele)
            x[col] = gt
    return x


#Function to z-score SNPs
def z_score(df): 
    return (df-df.mean())/df.std()


# Map label to numeric
def get_label(x, mapping):
    return mapping[x["eye_color"]]


# Main function
def main(args):
    
    # Read in filtered genotypes
    filtered_genotypes = pd.read_csv(args.genotype_file, sep='\t')
    
    # Replace -- with NaN
    filtered_genotypes = filtered_genotypes.replace('--', np.nan)

    # Read in the final SNP the final SNP set
    oneK_ref = pd.read_csv(args.rsids, sep='\t')

    # Keep only SNPs found in 1000Genomes
    gt_char = filtered_genotypes.merge(oneK_ref, left_on="rsid", right_on="ID")

    # WARNING: Getting the numbers takes about an hour and half
    gt_nums = gt_char.progress_apply(get_gt, axis=1)
    
    # Save a final SNP vs indiviudal table with other SNP info in it
    gt_nums = gt_nums.drop("ID", axis=1)
    gt_nums.to_csv("{}/openSNP_final_genotypes.tsv".format(args.output_dir), sep='\t', index=False)

    # Z-score SNPs
    z_scored_snps = gt_nums.set_index("rsid").loc[:, "6":"6131"].progress_apply(z_score, axis=1)
    z_scored_snps = z_scored_snps.replace(np.nan, 0)
    
    # Save SNPs included
    snps = z_scored_snps.index
    with open('{}/openSNP_final_rsids.txt'.format(args.output_dir), 'w') as f:
        f.writelines("%s\n" % id for id in snps)
    
    # Format data
    test_set = z_scored_snps.T
    
    # Reorder to match
    if args.snp_order != None:
        with open(args.snp_order, 'r') as filename:
            snp_order = [line.rstrip() for line in filename.readlines()]
        test_set = test_set[snp_order]
    
    # Save test set and test labels
    test_set.to_csv("{}/test_set.csv".format(args.output_dir), index=True)
    
    # Read in phenotypes
    phenotypes = pd.read_csv(args.phenotype_file, sep='\t')
    final_phenotypes = phenotypes[phenotypes["user_id"].isin(test_set.index)]
    
    # Use numeric labels for phenotypes
    label_mapping = {"brown":0, "blue":1, "green":2}
    final_phenotypes["label"] = final_phenotypes.apply(get_label, mapping=label_mapping, axis=1)
    final_phenotypes["label"].to_csv('{}/test_labels.csv'.format(args.output_dir), index=True)
                   
    # Write log info
    with open("{}/extractGenotypes.log".format(args.log_dir), 'w') as filename:
        filename.writelines("Test matrix: {} X {}\n".format(test_set.shape[1], test_set.shape[0]))
        

if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_filtered_genotypes.tsv")
    PHENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_initial_phenotypes.tsv")
    ONE_K_RSIDS = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "oneK_rsids.tsv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to openSNP filtered genotype tsv')
    parser.add_argument('-p', '--phenotype_file', type=str, default=PHENOTYPES, help='filepath to openSNP phenotypes tsv')
    parser.add_argument('-id', '--rsids', type=str, default=ONE_K_RSIDS, help='filepath to 1000Genomes rsids')
    parser.add_argument('-so', '--snp_order', type=str, default=None, help='filepath to ordering of SNPs to match')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)