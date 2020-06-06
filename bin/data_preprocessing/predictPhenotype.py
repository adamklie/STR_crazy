"""
Adam Klie
05/27/2020
Predicting phenotypes of 1000Genomes individuals using IrisPlex model
"""

import os
import pandas as pd
import argparse
import numpy as np
from sklearn.model_selection import train_test_split

# Functions for model probability calculations
def p_blue(m1_sum, m2_sum):
    return np.exp(m1_sum)/(1+np.exp(m1_sum)+np.exp(m2_sum))


def p_other(m1_sum, m2_sum):
    return np.exp(m2_sum)/(1+np.exp(m1_sum)+np.exp(m2_sum))


# Function to switch genotype label for SNPs whose ref allele matches iris
def fix_genotypes(x):
    if x["ID"] in ['rs12896399', 'rs12913832', 'rs16891982']:
        for col in x.index:
            if ("HG" in col) or ("NA" in col):
                x[col] = 2 - x[col]
    return x


# Function for getting the color prediction out of probabilities
def get_color_pred(x):
    if x["blue"] >= x["brown"] and x["blue"] >= x["other"]:
        return "blue"
    elif x["brown"] >= x["other"]:
        return "brown"
    else:
        return "other"
    
    
# Function to get a numeric label for predicted eye color
def get_label(x, mapping):
    return mapping[x["predicted_eye_color"]]
  
    
# Main function
def main(args):
    
    # Open log for reporting
    log = open("{}/predictPhenotypes.log".format(args.log_dir), 'w')
    
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
    
    # Get a numeric and string representation of predicted eye_color
    predict["predicted_eye_color"] = predict.apply(get_color_pred, 1)
    
    # Use numeric labels for phenotypes
    label_mapping = {"brown":0, "blue":1, "other":2}
    predict["label"] = predict.apply(get_label, mapping=label_mapping, axis=1)
   
    # Sanity check to see if matches homework
    color = predict[["blue", "other", "brown"]]
    for ind in ["NA12249", "NA20509", "NA12750"]:
        color_pred = color.loc[ind].idxmax(axis=1)
        log.writelines(ind + "\t" + color_pred + "\n")
    
    # Get labels and count
    labels = predict["label"]
    train_labels, val_labels = train_test_split(labels, test_size=0.2, stratify=labels.values)
    train_counts = train_labels.value_counts()
    val_counts = val_labels.value_counts()
    label_counts = labels.value_counts()
    
    # Save to files
    train_labels.to_csv('{}/train_labels.csv'.format(args.output_dir), index=True)
    val_labels.to_csv('{}/val_labels.csv'.format(args.output_dir), index=True)
    
    # Save train and val ids to txt files
    with open('{}/train_ids.txt'.format(args.output_dir), 'w') as f:
        f.writelines("%s\n" % id for id in train_labels.index)
              
    with open('{}/val_ids.txt'.format(args.output_dir), 'w') as f:
        f.writelines("%s\n" % id for id in val_labels.index)
              
    # Output to log
    log.writelines("{}\n".format(str(predict["predicted_eye_color"].value_counts())))
    log.writelines("{}\n".format(str(predict["predicted_eye_color"].value_counts().sum())))
    log.writelines("{}\n".format(str(predict.loc["NA12249"])))
    log.writelines("{}\n".format(str(predict.loc["NA20509"])))
    log.writelines("{}\n".format(str(predict.loc["NA12750"])))
    log.writelines("{}\n".format(str(label_counts)))
    log.writelines("{}\n".format(str(train_counts)))  
    log.writelines("{}\n".format(str(val_counts)))
    log.close()

if __name__ == '__main__':
    GENOTYPES = os.path.join(os.environ["HOME"], "project/datasets/oneKGenomes/data", "iris_oneK_genotypes.tsv")
    IRISPLEX = os.path.join(os.environ["HOME"], "project/datasets", "irisplex.bed")
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genotype_file', type=str, default=GENOTYPES, help='filepath to 1000Genomes genotype tsv')
    parser.add_argument('-ip', '--iris_plex_file', type=str, default=IRISPLEX, help='filepath to IrisPlex bed file')
    parser.add_argument('-o', '--output_dir', type=str, default='.', help='path to output dir')
    parser.add_argument('-l', '--log_dir', type=str, default='.', help='path to output log file to')
    args = parser.parse_args()
    main(args)