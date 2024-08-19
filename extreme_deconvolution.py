import numpy as np
from matplotlib import pyplot as plt

from astroML.density_estimation import XDGMM
from astroML.plotting.tools import draw_ellipse


# example
if False:
    # Sample the dataset. 
    # Here we use sample size = 400 in the example, 
    # which converges in shorter time, and gives reasonable result.
    N = 400
    np.random.seed(0)

    # generate the true data
    x_true = (1.4 + 2 * np.random.random(N)) ** 2
    y_true = 0.1 * x_true ** 2

    # add scatter to "true" distribution
    dx = 0.1 + 4. / x_true ** 2
    dy = 0.1 + 10. / x_true ** 2

    x_true += np.random.normal(0, dx, N)
    y_true += np.random.normal(0, dy, N)

    # define a function to plot all distributions in the same format
    def plot_distribution(text, sample_x, sample_y):
        plt.figure(figsize=(5, 3.75))
        plt.scatter(sample_x, sample_y, s=4,lw=0,c='k')
        plt.xlim(-1, 13)
        plt.ylim(-6, 16)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title(text,fontsize=10)

    # plot true distribution
    plot_distribution('True Distribution', x_true, y_true)

    plt.show()



    # add noise to get the "observed" distribution
    dx = 0.2 + 0.5 * np.random.random(N)
    dy = 0.2 + 0.5 * np.random.random(N)

    x = x_true + np.random.normal(0, dx)
    y = y_true + np.random.normal(0, dy)

    # plot noisy distribution
    plot_distribution('Noisy Distribution', x, y)

    plt.show()

    # stack the results for computation
    X = np.vstack([x, y]).T
    Xerr = np.zeros(X.shape + X.shape[-1:])
    diag = np.arange(X.shape[-1])
    Xerr[:, diag, diag] = np.vstack([dx ** 2, dy ** 2]).T

    clf = XDGMM(n_components=10, max_iter=200)

    clf.fit(X, Xerr)
    sample = clf.sample(N)

    # plot noisy distribution
    plot_distribution('Extreme Deconvolution Resampling', sample[:, 0], sample[:, 1])
    plt.show()

    # Plot the results
    fig = plt.figure(figsize=(5, 3.75))
    fig.subplots_adjust(left=0.1, right=0.95,
                        bottom=0.1, top=0.95,
                        wspace=0.02, hspace=0.02)

    ax1 = fig.add_subplot(221)
    ax1.scatter(x_true, y_true, s=4, lw=0, c='k')

    ax2 = fig.add_subplot(222)
    ax2.scatter(x, y, s=4, lw=0, c='k')

    ax3 = fig.add_subplot(223)
    ax3.scatter(sample[:, 0], sample[:, 1], s=4, lw=0, c='k')



    titles = ["True Distribution", "Noisy Distribution",
            "Extreme Deconvolution\n  resampling",
            "Extreme Deconvolution\n  cluster locations"]

    ax = [ax1, ax2, ax3]

    for i in range(len(ax)):
        ax[i].set_xlim(-1, 13)
        ax[i].set_ylim(-6, 16)

        ax[i].xaxis.set_major_locator(plt.MultipleLocator(4))
        ax[i].yaxis.set_major_locator(plt.MultipleLocator(5))

        ax[i].text(0.05, 0.95, titles[i],
                ha='left', va='top', transform=ax[i].transAxes)

        if i in (0, 1):
            ax[i].xaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax[i].set_xlabel('$x$')

        if i in (1, 3):
            ax[i].yaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax[i].set_ylabel('$y$')

    plt.show()






from utils import recombine_arrays
import pandas as pd
from paper_data import *
import random
from sklearn.neighbors import KernelDensity
import seaborn as sns

# get black hole masses:
# load QPE and TDE data
qpe_load_0 = np.loadtxt("QPE_allRelevantData_0.txt")
qpe_load_1 = np.loadtxt("QPE_allRelevantData_1.txt")
qpe_load_2 = np.loadtxt("QPE_allRelevantData_2.txt")
QPE_fullData = recombine_arrays(qpe_load_0,qpe_load_1,qpe_load_2)

tde_load_0 = np.loadtxt("TDE_allRelevantData_0.txt")
tde_load_1 = np.loadtxt("TDE_allRelevantData_1.txt")
tde_load_2 = np.loadtxt("TDE_allRelevantData_2.txt")
TDE_fullData = recombine_arrays(tde_load_0,tde_load_1,tde_load_2)
TDE_fullData = TDE_fullData[:10,:,:]

french_tde_load_0 = np.loadtxt("french_TDE_allRelevantData_0.txt")
french_tde_load_1 = np.loadtxt("french_TDE_allRelevantData_1.txt")
french_tde_load_2 = np.loadtxt("french_TDE_allRelevantData_2.txt")
french_TDE_fullData = recombine_arrays(french_tde_load_0,french_tde_load_1,french_tde_load_2)
french_TDE_fullData = french_TDE_fullData[:,:,:]


# load reference catalog
refCat = np.loadtxt(f"referenceCatalog_with_uncertainties_final_9bins.txt")
indices = list(range(len(refCat)))
random.shuffle(indices)
print(len(refCat))
print(indices[:1000])
refCat = refCat[indices[:1000]]
#fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
#refCat = pd.read_csv(f"referenceCatalog_with_uncertainties_final_9bins.txt", delimiter=" ", header=None, names=fieldnames)

QPE_stellar_masses, TDE_stellar_masses, french_TDE_stellar_masses = QPE_fullData[:,3,:], TDE_fullData[:,3,:], french_TDE_fullData[:,3,:]
QPE_mBH, TDE_mBH, french_TDE_mBH = QPE_fullData[:,4,:], TDE_fullData[:,4,:], french_TDE_fullData[:,4,:]
TDE_stellar_masses = np.concatenate((TDE_stellar_masses, french_TDE_stellar_masses))
TDE_redshifts = np.concatenate((TDE_redshifts, french_TDE_redshifts))
TDE_mBH = np.concatenate((TDE_mBH, french_TDE_mBH))
ref_mBH = refCat[:,100:102]
ref_redshifts = refCat[:,4]


# define a function to plot all distributions in the same format
def plot_distribution(sample_z, sample_mbh, err_mbh, host="QPE"):
    index = ["QPE","TDE","ref"].index(host)
    plt.errorbar(x=sample_z, y=sample_mbh, yerr=err_mbh, fmt=".", label=["QPE","TDE","ref"][index], color=["blue","red","grey"][index])

    

# plot true distribution
plt.figure(figsize=(5, 3.75))
plot_distribution(ref_redshifts, ref_mBH[:,0], ref_mBH[:,1], host="ref")
plot_distribution(QPE_redshifts, QPE_mBH[:,0], [QPE_mBH[:,1],QPE_mBH[:,2]], host="QPE")
plot_distribution(TDE_redshifts, TDE_mBH[:,0], [TDE_mBH[:,1],TDE_mBH[:,2]], host="TDE")

plt.xlabel('$z$')
plt.ylabel('$M_\mathrm{BH}$')
plt.legend()
plt.title("Distributions",fontsize=10)
plt.show()



def extreme_deconvolve(x, y, dx, dy, n_components=10, max_iter=200):
    X = np.vstack([x, y]).T
    Xerr = np.zeros(X.shape + X.shape[-1:])
    diag = np.arange(X.shape[-1])
    Xerr[:, diag, diag] = np.vstack([dx ** 2, dy ** 2]).T
    clf = XDGMM(n_components=n_components, max_iter=max_iter)
    clf.fit(X, Xerr)
    sample = clf.sample(len(X))
    return sample

QPE_mBH_dc = extreme_deconvolve(x=QPE_redshifts, y=QPE_mBH[:,0], dx=np.zeros(np.array(QPE_redshifts).shape), dy=QPE_mBH[:,1], n_components=1)[:,1]
TDE_mBH_dc = extreme_deconvolve(x=TDE_redshifts, y=TDE_mBH[:,0], dx=np.zeros(np.array(TDE_redshifts).shape), dy=TDE_mBH[:,1], n_components=1)[:,1]
ref_mBH_dc = extreme_deconvolve(x=ref_redshifts, y=ref_mBH[:,0], dx=np.zeros(np.array(ref_redshifts).shape), dy=ref_mBH[:,1], n_components=5)[:,1]

# plot noisy distribution
plt.plot(ref_redshifts, ref_mBH_dc, "o", color="grey", label="ref")
plt.plot(QPE_redshifts, QPE_mBH_dc, "o", color="blue", label="QPE")
plt.plot(TDE_redshifts, TDE_mBH_dc, "o", color="red", label="TDE")
plt.xlabel('$z$')
plt.ylabel('$M_\mathrm{BH}$')
plt.legend()
plt.title("Deconvolved distributions",fontsize=10)
plt.show()

import sys
sys.exit()

# Compare before/after histograms
ax_before = plt.subplot(211)
ax_after = plt.subplot(212)
# m_bh
x_min, x_max = (4.5,9)
X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
smoothness=12
referenceSmoothness=12
bandwidth = np.abs(x_max-x_min)/smoothness

kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[2]}"])[:,np.newaxis])
log_dens = kde.score_samples(X_plot)
#mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
ax_before.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=2)
kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[4,:,0][:,np.newaxis])
log_dens = kde.score_samples(X_plot)
ax_before.fill_betweenx(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[4,:,0][:,np.newaxis])
log_dens = kde.score_samples(X_plot)
ax_before.fill_betweenx(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)