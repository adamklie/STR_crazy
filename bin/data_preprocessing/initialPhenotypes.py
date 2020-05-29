"""
Adam Klie
05/27/2020
Retrieve user IDs with eye color phenotype data from OpenSNP
Input: 
    - phenotype data for downloaded openSNP genotype files
Output: 
    - openSNP_initial_phenotypes.tab: a tabular version of initial set of users with phenotype data (638)
    - openSNP_initial_userid.txt: list of user IDs (one per line) for those users with eye color phenotype data
Notes:
    - Only looking at individuals with the following eye colors: 'Brown', 'brown', 'Blue', 'blue', 'Green', 'green'
    - In the following columns: "Eye color", "Eye Color", "Eye pigmentation "
    - 638 such individuals come out of this
"""

import os
import pandas as pd
import argparse


# Collapes to a single column
def collapseEyeColors(x, values):
    if x["Eye color"] not in values:
        if x["Eye Color"] in values:
            return x["Eye Color"].lower()
        elif x["Eye pigmentation "] in values:
            return x["Eye pigmentation "].lower()
    else:
        return x["Eye color"].lower()

    
# Main function
def main(args):
    
    # Read in OpenSNP phenotypes
    phenotypes = pd.read_csv(args.phenotype_file, sep=';')
    phenotypes = phenotypes[['user_id', 'genotype_filename', 'Eye color', 'Eye Color', 'Eye pigmentation ']]

    # Remove non-reporters
    phenotypes = phenotypes[(phenotypes["Eye color"] != '-') | (phenotypes["Eye Color"] != '-') | (phenotypes["Eye pigmentation "] != '-')]

    # Valid phenotype labels to keep
    valid_phenotypes = ['Brown', 'brown', 'Blue', 'blue', 'Green', 'green']

    # Only take those samples that have valid phenotypes in one of these three columns
    phenotypes = phenotypes[(phenotypes["Eye color"].isin(valid_phenotypes)) | 
                            (phenotypes["Eye Color"].isin(valid_phenotypes)) | 
                            (phenotypes["Eye pigmentation "].isin(valid_phenotypes))]


    # Use function to collapse eye color metadata to single columns
    phenotypes["eye_color"] = phenotypes.apply(lambda x: collapseEyeColors(x, valid_phenotypes), 1)

    # Remove duplicate individuals
    phenotypes = phenotypes.drop_duplicates(subset="user_id")

    # Only take Ancestry and 23andMe data to simplify
    phenotypes = phenotypes[phenotypes["genotype_filename"].str.contains('ancestry') | phenotypes["genotype_filename"].str.contains('23andme')]
    
    # Save these initial phenotypes to a file
    phenotypes.to_csv('{}/openSNP_initial_phenotypes.tsv'.format(args.output_dir), sep='\t', index=False)

    # Save this id set as a list
    ids = set(phenotypes["user_id"].values)
    with open("{}/openSNP_initial_userids.txt".format(args.output_dir), 'w') as filehandle:
        filehandle.writelines("%s\n" % id for id in ids)
    
    # Save potentially useful log information
    with open("{}/initialPhenotypes.log".format(args.log_dir), 'w') as log:
        log.writelines("Counts for each phenotype:\n{}\n".format(str(phenotypes["eye_color"].value_counts())))
        log.writelines("Total number of userids: {}\n".format(phenotypes["eye_color"].value_counts().sum()))

        
if __name__ == '__main__':
    PHENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP", "phenotypes_202004220659.csv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--phenotype_file', type=str, default=PHENOTYPES, help='filepath to openSNP phenotype csv')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)