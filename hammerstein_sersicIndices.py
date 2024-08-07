import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import seaborn as sns
from utils import myCombinedFinalPlot
import pandas as pd

hammerstein_TDE_data = np.loadtxt("hammerstein_TDE_distribution.txt")
QPE_data = np.loadtxt("QPE_distribution.txt")
refCat = np.loadtxt("referenceCatalog_final.txt")
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
refCat = pd.read_csv("referenceCatalog_final.txt", delimiter=" ", header=None, names=fieldnames)

hammerstein_Sersics = hammerstein_TDE_data[:,0]
QPE_Sersics = QPE_data[:,1]

x_min = np.min(QPE_Sersics) if np.min(hammerstein_Sersics) > np.min(QPE_Sersics) else np.min(hammerstein_Sersics)
x_max = np.max(QPE_Sersics) if np.max(hammerstein_Sersics) < np.max(QPE_Sersics) else np.max(hammerstein_Sersics)
X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
bandwidth = np.abs(x_max-x_min)/6

ax1 = plt.subplot(111)
kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/6).fit(np.array(refCat[f"col_60"])[:,np.newaxis])
log_dens = kde.score_samples(X_plot)
ax1.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=2)
kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_Sersics[:,np.newaxis])
log_dens = kde.score_samples(X_plot)
ax1.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(hammerstein_Sersics[:,np.newaxis])
log_dens = kde.score_samples(X_plot)
ax1.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

ax1.set_xlabel("Sérsic Index $n$")
plt.show()


# n vs mBH

plt.plot(QPE_Sersics, QPE_data[:,4], "o", color="blue", label="QPEs", markersize=8)

plt.plot(hammerstein_Sersics, hammerstein_TDE_data[:,1], "*", color="red", label="hammerstein TDEs", markersize=7)
plt.legend()
plt.ylabel("$M_\mathrm{BH}$", fontsize=16)
plt.xlabel("Sérsic index $n$", fontsize=16)
plt.show()