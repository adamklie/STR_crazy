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

# Function to switch genotype label for SNPs whose ref allele matches iris
def fix_genotypes(x):
    if x["rsid"] in ['rs12896399', 'rs12913832', 'rs16891982']:
        for col in x.index[3:-8]:
            x[col] = 2 - x[col]
    return x


# Functions for model probability calculations
def p_blue(m1_sum, m2_sum):
    return np.exp(m1_sum)/(1+np.exp(m1_sum)+np.exp(m2_sum))


def p_other(m1_sum, m2_sum):
    return np.exp(m2_sum)/(1+np.exp(m1_sum)+np.exp(m2_sum))


def get_color_pred(x):
    if x["blue"] >= x["brown"] and x["blue"] >= x["other"]:
        return "blue"
    elif x["brown"] >= x["other"]:
        return "brown"
    else:
        return "other"
    

def main():
    
    # Open log for reporting
    log = open("{}/testIris.log".format(args.log_dir), 'w')
    
    # Load irisplex parameters for each SNP, sort on ID
    iris_header = ["chr", "pos1", "pos2", "id", "minor_allele", "b1", "b2"]
    iris = pd.read_csv(args.iris_plex_file, sep='\t', header=None, names=iris_header)

    # Load the genotypes of individuals at each SNP, sort on ID and reconfigure genotypes for 3 alleles
    gt = pd.read_csv(args.genotype_file, sep='\t')
    merged_gt = gt.merge(iris, left_on="ID", right_on="id")
    fixed_gt = merged_gt.apply(fix_genotypes, 1)
    
    # Perform predictions for each individual based on genotype and parameters
    a1 = 3.94
    a2 = 0.65

    predictions = {}
    for col in fixed_gt.columns:
        if ("HG" in col) or ("NA" in col):
            pred1 = np.dot(fixed_gt[col], fixed_gt["b1"]) + a1
            pred2 = np.dot(fixed_gt[col], fixed_gt["b2"]) + a2
            blue = p_blue(pred1, pred2)
            other = p_other(pred1, pred2)
            brown = 1 - blue - other
            predictions[col] = [pred1, pred2, blue, other, brown]

    predict = pd.DataFrame.from_dict(predictions, orient='index', columns=["pred1", "pred2", "blue", "other", "brown"])
    
    # Get the predicted phenotype
    predict["predicted_eye_color"] = predict.apply(get_color_pred, 1)
    
    # Read in actual phenotypes
    phenotypes = pd.read_csv(args.phenotype_file, sep='\t')
    
    # Merge actual and predicted phenotypes and determine correctness
    result = phenotypes.merge(predictions, on="user_id")
    result["pred_correct"] = result_nona["eye_color"].str.lower().values == result_nona["predicted_eye_color"].str.lower().values
    
    result.groupby("eye_color").mean()
    np.mean(result_nona["eye_color"].str.lower().values == result_nona["predicted_eye_color"].str.lower().values)

    # Drop those individuals that did not have full iris info and repeat
    result_nona = result.dropna()
    result_nona["pred_correct"] = result_nona["eye_color"].str.lower().values == result_nona["predicted_eye_color"].str.lower().values

    result_nona.groupby("eye_color").mean()
    np.mean(result_nona["eye_color"].str.lower().values == result_nona["predicted_eye_color"].str.lower().values)

    
if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_genotypes.tsv")
    PHENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "iris_oneK_genotypes.tsv")
    IRISPLEX = os.path.join(os.environ["HOME"], "project/datasets", "irisplex.bed")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-p', '--phenotype_file', type=str, default=PHENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-ip', '--iris_plex_file', type=str, default=IRISPLEX, help='filepath to IrisPlex bed file')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)