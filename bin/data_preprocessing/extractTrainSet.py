"""
Adam Klie
05/27/2020
Extract a test set and labels from openSNP for input into neural network
"""

import os
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
tqdm.pandas()


#Function to z-score SNPs
def z_score(x): 
    return (x-x.mean())/x.std()


# Main function
def main(args):
    
    # Format dataset and z-score
    print("Reading in genotype file and Z-scoring SNPs")
    data = pd.read_csv(args.genotype_file, sep='\t').drop_duplicates("ID").set_index("ID").loc[:, "HG00096":"NA21144"].progress_apply(z_score, axis=1)
    data = data.replace(np.nan, 0)
    
    # Reorder to match
    print("Subsetting SNPs from openSNP")
    if args.snp_order != None:
        with open(args.snp_order, 'r') as filename:
            snp_order = [line.rstrip() for line in filename.readlines()]
        data = data.loc[snp_order]
    data = data.T
             
    # Split data into train and val
    print("Train/val split")
    with open(args.train_ids, 'r') as filename:
        train_id = [line.rstrip() for line in filename.readlines()]
    with open(args.val_ids, 'r') as filename:
        val_id = [line.rstrip() for line in filename.readlines()]      
    train = data.loc[train_id]
    val = data.loc[val_id]
    
    print("Writing to sets")    
    # Write to a file
    train.to_pickle('{}/train_set.pickle'.format(args.output_dir))
    val.to_pickle('{}/val_set.pickle'.format(args.output_dir))
    
    # Write log info
    with open("{}/extractTrainSet.log".format(args.log_dir), 'w') as filename:
        filename.writelines("Train dimensions: {} X {}\n".format(train.shape[1], train.shape[0]))
        filename.writelines("Val dimensions: {} X {}\n".format(val.shape[1], val.shape[0]))
        
    
if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "oneK_genotypes.tsv")
    TRAIN_IDS = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "train_ids.txt")
    VAL_IDS = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "val_ids.txt")
    SNP_ORDERING = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_final_rsids.txt")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-t', '--train_ids', type=str, default=TRAIN_IDS, help='filepath to training id labels')
    parser.add_argument('-v', '--val_ids', type=str, default=VAL_IDS, help='filepath to val id labels')
    parser.add_argument('-so', '--snp_order', type=str, default=SNP_ORDERING, help='filepath to ordering of SNPs to match OpenSNP test set')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)
