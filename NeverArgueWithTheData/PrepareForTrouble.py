'''
@author = james

Functions used to generate ROC plots for openSNP test set. If want to save instead of outputting to nb simply pass in save=True to TeamROCket.

Example to run (note first arg==Loaded trained model for first line):

whatAreMyScores = np.array(Meowth(whoeverMadeQualExamsAThingHasEarnedMyEternalIre, test_loader, pleaseBeGpu))
phenotypeLabels = list(pd.read_csv("../Sheen/test_labels.csv").loc[:, "label"])
irisPlexScores = np.array(pd.read_csv("../Sheen/openSNP_final_iris_preds.tsv", sep="\t").loc[:, "brown":"other"])

TeamROCket(whatAreMyScores, phenotypeLabels, "Get Gymwrecked")
TeamROCket(irisPlexScores, phenotypeLabels, "IrisPlex")


'''

import matplotlib.pyplot as plt
from itertools import cycle

from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score

def Meowth(model, loader, pleaseBeGpu):
    model.eval()
    predictedScores = list()
    #orderOfPheno = will just be pulled from the CSV of labels here: phenotype_file="../Sheen/test_labels.csv", 
    for i, (snpBatch, phenotypeBatch) in enumerate(loader): 
        #print(*phenotypeBatch.tolist())
        snpBatch = snpBatch.to(pleaseBeGpu)
        phenotypeBatch = phenotypeBatch.to(pleaseBeGpu)
        output = model(snpBatch)
        for el in output[0].to('cpu').detach().numpy():
            predictedScores.append(el)
     
    
    return predictedScores 

def TeamROCket(scores, labels, name, save=False):
    y = label_binarize(labels, classes=[0, 1, 2])
    n_classes = y.shape[1]

    y_score = scores

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    
    plt.figure(figsize=(7,7))
    lw=2.3
    
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='indigo', linestyle=':', linewidth=4)

    colors = cycle(['cyan', 'green', 'red'])
    for i, color in zip(range(n_classes), colors):   
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(name + ' Test Set ROC')
    plt.legend(loc="lower right")
    plt.show()
    if save:
        plt.savefig(name + "_ROC.png")
        
    return None
