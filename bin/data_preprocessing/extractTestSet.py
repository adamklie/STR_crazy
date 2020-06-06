"""
Adam Klie
05/27/2020
Z-score a SNP x individual matrix and return a usable test set for model training
"""

import os
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
tqdm.pandas()


#Function to z-score SNPs based on input dataframe with standard deviations and means
def z_score(x, z_df):
    mean = z_df["means"].loc[x.name]
    std =  z_df["std"].loc[x.name]
    return ((x-mean)/std)


# Map label to numeric
def get_label(x, mapping):
    return mapping[x["eye_color"]]


# Main function
def main(args):
    
    # Read in final genotypes
    print("Reading in genotypes")
    raw_genotypes = pd.read_csv(args.genotype_file, sep='\t').drop_duplicates("rsid").set_index("rsid").loc[:, "6":"6131"]
    print(len(raw_genotypes))

    # Read in SNP means and stds
    print("Reading in means and stds from 1000K genomes")
    mean_std = pd.read_csv(args.stats, sep='\t', index_col=0).loc[raw_genotypes.index]
    print(len(mean_std))
          
    # Z-score SNPs
    print("Z-scoring SNPs")
    z_genotypes = raw_genotypes.progress_apply(z_score, z_df=mean_std, axis=1)
    z_genotypes = z_genotypes.replace(np.nan, 0)
    
    # Save SNPs included
    print("Saving SNP list")
    snps = z_genotypes.index
    with open('{}/openSNP_final_rsids.txt'.format(args.output_dir), 'w') as f:
        f.writelines("%s\n" % id for id in snps)
    
    # Format data as SNPs in rows and individuals in columns
    test_set = z_genotypes.T
    
    # Save test set and test labels
    print("Saving test set as pickle object")
    test_set.to_pickle("{}/test_set.pickle".format(args.output_dir))
    
    # Read in phenotypes
    print("Reading in phenotpyes")
    phenotypes = pd.read_csv(args.phenotype_file, sep='\t')
    final_phenotypes = phenotypes[phenotypes["user_id"].isin(test_set.index)]
    
    # Use numeric labels for phenotypes
    print("Saving test labels as csv")
    label_mapping = {"brown":0, "blue":1, "green":2}
    final_phenotypes["label"] = final_phenotypes.apply(get_label, mapping=label_mapping, axis=1)
    final_labels = final_phenotypes.set_index("user_id")["label"]
    final_labels.to_csv('{}/test_labels.csv'.format(args.output_dir), index=True)
                   
    # Write log info
    with open("{}/extractTestSet.log".format(args.log_dir), 'w') as filename:
        filename.writelines("Initial SNPs in {}: {}\n".format(args.genotype_file, len(raw_genotypes)))
        filename.writelines("Test dataset size: {} X {}\n".format(test_set.shape[1], test_set.shape[0]))
        filename.writelines("Test labels size: {}\n".format(len(final_labels)))
        

if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_final_genotypes.tsv")
    PHENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_initial_phenotypes.tsv")
    STATS = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "train_set_stats.tsv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to openSNP filtered genotype tsv, default is {}'.format(GENOTYPES))
    parser.add_argument('-p', '--phenotype_file', type=str, default=PHENOTYPES, help='filepath to openSNP phenotypes tsv, default is {}'.format(PHENOTYPES))
    parser.add_argument('-s', '--stats', type=str, default=STATS, help='means and standard devs of SNPs from training set, default is {}'.format(STATS))
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)