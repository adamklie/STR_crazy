import pandas as pd
import torch
from torch.utils import data
import numpy as np

class SNPDataset(data.Dataset):
    """Characterizes a dataset for PyTorch."""

    def __init__(self, geno_file, pheno_file, pickle=False):
        """Initialization.
        Args:
            geno_file
            pheno_file
        """
        if pickle == False:
            self.genotypes = pd.read_csv(geno_file, index_col=0)
        else:
            self.genotypes = pd.read_pickle(geno_file)
        
        self.phenotypes = pd.read_csv(pheno_file, index_col=0)
        self.list_ids = self.phenotypes.index

    def __len__(self):
        return len(self.list_ids)

    def __getitem__(self, index):
        """Returns one data pair (genotypes and label)."""
        ID=self.list_ids[index]
        
        #load data
        X = torch.from_numpy(np.array(self.genotypes.loc[ID].values)).float()
        y = torch.from_numpy(np.array(self.phenotypes.loc[ID].values))
        return X,y


def get_loader(genotype_file, phenotype_file, batch_size, shuffle, num_workers, pickle=False):
    """Returns torch.utils.data.DataLoader for geno dataset."""

    geno = SNPDataset(geno_file=genotype_file, pheno_file=phenotype_file, pickle=pickle)
    params = {'batch_size': batch_size, 'shuffle': shuffle,'num_workers': num_workers}
    data_loader = torch.utils.data.DataLoader(dataset=geno,**params)
    return data_loader