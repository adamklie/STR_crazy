"""
Adam Klie
05/27/2020
Test IrisPlex on openSNP data
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
    

def main(args):
    
    # Open log for reporting
    log = open("{}/testIris.log".format(args.log_dir), 'w')
    
    """
    # Load irisplex parameters for each SNP, sort on ID
    iris_header = ["chr", "pos1", "pos2", "id", "minor_allele", "b1", "b2"]
    iris = pd.read_csv(args.iris_plex_file, sep='\t', header=None, names=iris_header)

    # Load the genotypes of individuals at each SNP, sort on ID and reconfigure genotypes for 3 alleles
    gt = pd.read_csv(args.genotype_file, sep='\t')
    merged_gt = gt.merge(iris, left_on="ID", right_on="id")
    fixed_gt = merged_gt.apply(fix_genotypes, 1)
    """
    
    # Read in already merged table
    fixed_gt = pd.read_csv(args.genotype_file, sep='\t')
    
    # Perform predictions for each individual based on genotype and parameters
    a1 = 3.94
    a2 = 0.65

    predictions = {}
    for col in fixed_gt.columns[3:-8]:
            pred1 = np.nansum(fixed_gt[col]*fixed_gt["b1"]) + a1
            pred2 = np.nansum(fixed_gt[col]*fixed_gt["b2"]) + a2
            blue = p_blue(pred1, pred2)
            other = p_other(pred1, pred2)
            brown = 1 - blue - other
            predictions[col] = [pred1, pred2, blue, other, brown]

    predict = pd.DataFrame.from_dict(predictions, orient='index', columns=["pred1", "pred2", "blue", "other", "brown"])
    
    # Get the predicted phenotype
    predict["predicted_eye_color"] = predict.apply(get_color_pred, 1)
    predict = predict.reset_index()
    print(predict.head())
    predict['user_id'] = predict['index'].astype(int)
    
    # Read in actual phenotypes
    phenotypes = pd.read_csv(args.phenotype_file, sep='\t')

    # Merge actual and predicted phenotypes and determine correctness
    results = phenotypes.merge(predict, on="user_id")
    results_nona = results.dropna()
    results["pred_correct"] = results["eye_color"].str.lower().values == results["predicted_eye_color"].str.lower().values
    results_nona["pred_correct"] = results_nona["eye_color"].str.lower().values == results_nona["predicted_eye_color"].str.lower().values
    results.to_csv('{}/openSNP_iris_prediction.tsv'.format(args.output_dir), sep='\t', index=False)
    
    log.writelines("Overally accuracy of IrisPlex on openSNP: {}\n".format(np.mean(results["pred_correct"])))
    log.writelines("By label accuracy of IrisPlex on openSNP: {}\n".format(results.groupby("eye_color").mean()["pred_correct"]))
    log.writelines("Accuracy on individuals that contain all 6 SNPs: {}\n".format(np.mean(results_nona["pred_correct"])))
    log.writelines("By label accuracy of of indiviudals that contain all 6 SNPs: {}\n".format(results_nona.groupby("eye_color").mean()["pred_correct"]))
    log.close()
    
    
if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_iris_genotypes.tsv")
    PHENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/openSNP/data", "openSNP_initial_phenotypes.tsv")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-p', '--phenotype_file', type=str, default=PHENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)