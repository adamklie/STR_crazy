import pandas as pd
import torch
from torch.utils import data
import numpy as np

class SNPDataset(data.Dataset):
    """Characterizes a dataset for PyTorch."""

    def __init__(self, geno_file, pheno_file):
        """Initialization.
        Args:
            geno_file
            pheno_file
        """
        self.genotypes = pd.read_csv("../../datasets")
        self.file = geno_file
        self.labels = labels
        self.labels2 = labels2

    def __len__(self):
        return len(self.list_ids)

    def __getitem__(self, index):
        """Returns one data pair (genotypes and label)."""
        ID=self.list_ids[index]
        #load data
        with open(self.file) as csvfile:
            file = csv.reader(csvfile)  # read in genotypes
            for row in file:
                if row[1]==ID:
                    X=torch.from_numpy(np.array([int(x) for x in row[6:]])).float()
                    X=torch.unsqueeze(X,0)
            y=torch.from_numpy(np.array(self.labels[ID]))
            z=torch.from_numpy(np.array(self.labels2[ID]))
            return X,y,z


def get_loader(geno_file,ids,labels,batch_size, shuffle, num_workers):
    """Returns torch.utils.data.DataLoader for geno dataset."""

    geno = Dataset(geno_file=geno_file,list_ids=ids,labels=labels,labels2=labels)
    params = {'batch_size': batch_size, 'shuffle': shuffle,'num_workers': num_workers}
    data_loader = torch.utils.data.DataLoader(dataset=geno,**params)
    return data_loader