import numpy as np


#QPE_distribution = np.loadtxt("QPE_dist.txt")
#TDE_distribution = np.loadtxt("TDE_dist.txt")
sample_distribution = np.loadtxt("referenceCatalog_modif.txt")
sample_distribution = sample_distribution[:,[12,60,68]]

def getStats(dist):
    n_data = dist.shape[0] # number of data points
    means = np.mean(dist, axis=0)
    print(means)

getStats(sample_distribution)