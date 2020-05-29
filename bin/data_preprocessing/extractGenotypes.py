"""
Adam Klie
05/22/2020
Code to extract an initial SNP (rsID) x individuals genotype table from OpenSNP data
Input: 
    - directory with OpenSNP genotypes data
    - list of userids
Output: 
    - openSNP_initial_genotypes.tsv: tab deliminated table of SNPs vs individuals
"""

import pandas as pd
import os
import glob
import numpy as np
import tqdm
import re
import argparse


# Main function
def main(args):
    
    # Extract sample ids as list
    with open(args.id_file) as f:
        ids = [line.rstrip() for line in f.readlines()]
    
    # Log file for reporting
    log = open("{}/extractGenotypes.log".format(args.log_dir), 'w')
    
    # Open directory with genotype files
    os.chdir(args.genotypes_dir)

    # Create a directory to hold genotypes
    merged = pd.DataFrame()

    # Useful counters
    num_bad_files = 0
    num_valid_files = 0
    
    # Optional subset
    if args.subset != None:
        ids = ids[:args.subset]
        
    # Loop through the files
    for i, id in enumerate(tqdm.tqdm(ids)):
        filename = glob.glob("user{}_*".format(id))[0]

        # Check for valid 23andMe filetype
        if "23andme.txt" in filename:
            try:
                current = pd.read_csv(filename, sep='\t', names=['rsid', 'chromosome', 'position', 'gt'], comment='#')
            except:
                log.writelines("Error when reading {} as csv\n".format(filename))
                num_bad_files += 1
                continue
            else:
                current[str(id)] = current["gt"]

        # Check for valid AncestryDNA filetype
        elif "ancestry.txt" in filename:
            try:
                current = pd.read_csv(filename, sep='\t', comment='#')
            except:
                log.writelines("Error when reading {} as csv\n".format(filename))
                num_bad_files += 1
                continue
            else:
                current[str(id)] = current["allele1"] + current["allele2"]

        # Report unused filetype
        else:
            log.writelines("Error: filetype not currently used: {}\n".format(filename))
            num_bad_files += 1
            continue

        # If this is the first file set-up the dataframe
        if i == 0:
            merged = current[["rsid", "chromosome", "position", str(id)]]

        # Else add the column
        else:
            current = current[["rsid", str(id)]]
            merged = merged.merge(current, how="left", on="rsid")
        num_valid_files += 1

    # Get only valid rsid SNPs in list
    rsid_re = re.compile('rs[0-9]+')
    merged = merged[merged["rsid"].str.match(rsid_re)]

    # Write the dataframe to a file
    merged.to_csv('{}/openSNP_initial_genotypes.tsv'.format(args.output_dir), sep='\t', index=False, na_rep=np.nan)
    
    # Report number of invalid files
    log.writelines("Number of bad file formats: {}\n".format(num_bad_files))
    log.writelines("Number of valid files: {}\n".format(num_valid_files))
    log.writelines("Number of SNPs captured: {}\n".format(len(merged)))
    log.close()

if __name__ == '__main__':
    USERIDS = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_initial_userids.txt")
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/genotypes/")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--id_file', type=str, default=USERIDS, help='filepath to openSNP userid list')
    parser.add_argument('-g', '--genotypes_dir', type=str, default=GENOTYPES, help='filepath to openSNP genotype files')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    parser.add_argument('-s', '--subset', type=int, default=None, help='number of ids to subset (optional)')
    args = parser.parse_args()
    main(args)