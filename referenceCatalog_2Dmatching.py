"""
Trying to sample from a 2D distribution in order to match a second 2D distribution, will be fun!
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from paper_data import QPE_redshifts
import scipy as sp

# Let's load the background galaxies
refCat = np.loadtxt("referenceCatalog_modif2.txt")
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
fieldnames[1] = "redshift"
fieldnames[63] = "m_star"
refCat = pd.read_csv("referenceCatalog_modif2.txt", delimiter=" ", header=None, names=fieldnames)

# Let's load the QPE hosts:
fieldnames = ["bt_ratio", "n_sersic", "smsd", "m_star", "m_bh", "redshift"]
QPEGalaxies = pd.read_csv("QPE_distribution.txt", header=None, names=fieldnames, delimiter=" ")


# ***************************************************************************
# All parameters
# Choose bounds for that z-M_star grid
bounds = {"redshift": (0.01,0.1),
          "m_star": (9.5,10.4)}
resolution = 50
smoothing_param = 8
# ***************************************************************************



redshift_bins = np.linspace(bounds["redshift"][0], bounds["redshift"][1], resolution)
m_star_bins = np.linspace(bounds["m_star"][0], bounds["m_star"][1], resolution)

QPEGalaxies["redshift_bin"] = pd.cut(QPEGalaxies["redshift"], bins=redshift_bins)
QPEGalaxies["m_star_bin"] = pd.cut(QPEGalaxies["m_star"], bins=m_star_bins)

bin_counts = QPEGalaxies.groupby(['redshift_bin', 'm_star_bin'], observed=False).size()

#print(QPEGalaxies["redshift_bin"])
#print(QPEGalaxies["m_star_bin"])
#print(bin_counts)
QPE_grid = np.zeros((resolution,resolution))
for i in range(len(QPEGalaxies)):
    m_star, z = QPEGalaxies["m_star"][i], QPEGalaxies["redshift"][i]
    for mbin in list(m_star_bins)[::-1]:
        if m_star > mbin:
            m_star_index = list(m_star_bins).index(mbin)
            break
    for zbin in list(redshift_bins)[::-1]:
        if z > zbin:
            z_index = list(redshift_bins).index(zbin)
            break
    #print(m_star_index, z_index)
    QPE_grid[m_star_index,z_index] += 1

QPE_grid = QPE_grid/np.sum(QPE_grid) # normalize
print(np.sum(QPE_grid))
plt.imshow(QPE_grid, origin="lower", cmap="Greys")
plt.xlabel("$z$ bin", fontsize=15)
plt.ylabel("$M_\star$ bin", fontsize=15)
plt.title("Normalized probability density", fontsize=16)
plt.show()

# smooth it


QPE_grid = sp.ndimage.gaussian_filter(QPE_grid, sigma=smoothing_param, mode='constant')
QPE_grid = QPE_grid/np.sum(QPE_grid) # normalize
print(np.sum(QPE_grid))
plt.imshow(QPE_grid, origin="lower", cmap="Greys")
plt.xlabel("$z$ bin", fontsize=15)
plt.ylabel("$M_\star$ bin", fontsize=15)
plt.title("Normalized probability density", fontsize=16)
plt.show()


# Make that weird grid

# Sort the reference galaxies in their respective bins: