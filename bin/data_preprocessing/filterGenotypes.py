"""
Adam Klie
05/27/2020
Filtering for SNPs in a large fraction individuals
Input: 
    - openSNP_initial_genotypes.tsv: genotype data for downloaded openSNP genotype files
Output: 
    - openSNP_filtered_genotypes.tab: a tabular version of genotypes filtered
    - openSNP_filtered_rsids.txt: list of rsids to pass to 1000Genomes to find intersection with
Notes:
    - Filter this dataframe to bring the number of SNPs down (currently at least 80% of samples have rsid as default)
"""

import pandas as pd
import os
import numpy as np
import argparse


# Main function
def main(args):

    # Read in genotypes
    initial_genotypes = pd.read_csv(args.genotype_file, sep='\t')
    
    # For reporting later
    num_SNPs = len(initial_genotypes)
    num_samples = len(initial_genotypes.loc[:, "6":"6131"].columns)

    # Filter for SNPs that appear in at least 80% of samples (> than min iris SNP fraction)
    filtered_genotypes = initial_genotypes[initial_genotypes.loc[:,'6':].notna().apply(np.mean, 1) > args.percentage]

    # Write this filtered set to a smaller dataframe
    filtered_genotypes.to_csv('{}/openSNP_filtered_genotypes.tsv'.format(args.output_dir), sep='\t', index=False)

    # Save this list of SNP rsids for oneK_genotpyes.sh
    filtered_rsids = filtered_genotypes["rsid"].values.tolist()
    with open("{}/openSNP_filtered_rsids.txt".format(args.output_dir), 'w') as filehandle:
        filehandle.writelines("%s\n" % id for id in filtered_rsids)
    
    with open("{}/filterGenotypes.log".format(args.log_dir), 'w') as log:
        log.writelines("After extracting genotypes we have {} SNP x {} sample matrix\n".format(num_SNPs, num_samples))
        log.writelines("After filtering, %d SNPs and %d samples left" % (len(filtered_rsids), len(user_id_list)))

        
if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_initial_genotypes.tsv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to openSNP genotype tsv')
    parser.add_argument('-p', '--percentage', type=float, default=0.8, help='percentage of individuals necessary to keep SNP')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)