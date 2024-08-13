"""
Trying to sample from a 2D distribution in order to match a second 2D distribution, will be fun!
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from paper_data import QPE_redshifts
import scipy as sp
import sys
from tqdm import tqdm
from utils import cut_from_catalog

# ***************************************************************************
# All parameters
# Choose bounds for that z-M_star grid
bounds = {"redshift": (0.01,0.1),
          "m_star": (9.2,10.4),
          "m_bh": (4.5,8)}
resolution = 20
smoothing_param = resolution/6
sample_size=48000
show_plots = True
matching = "m_star" # choose wether to match on black hole mass or stellar mass
# ***************************************************************************





# Let's load the background galaxies
refCat = np.loadtxt("referenceCatalog_modif2.txt")
# Cut reference catalog to only have galaxies within selected bounds:
refCat = cut_from_catalog(refCat, 1, bounds["redshift"], verbose=True)
refCat = cut_from_catalog(refCat, (63 if matching == "m_star" else 67), bounds[matching], verbose=True)
np.savetxt("referenceCatalog_matching.txt", refCat)
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
fieldnames[1] = "redshift"
fieldnames[63] = "m_star"
fieldnames[67] = "m_bh"
refCat = pd.read_csv("referenceCatalog_matching.txt", delimiter=" ", header=None, names=fieldnames)

# Let's load the QPE hosts:
fieldnames = ["bt_ratio", "n_sersic", "smsd", "m_star", "m_bh", "redshift"]
QPEGalaxies = pd.read_csv("QPE_distribution.txt", header=None, names=fieldnames, delimiter=" ")






redshift_bins = np.linspace(bounds["redshift"][0], bounds["redshift"][1], resolution)
m_star_bins = np.linspace(bounds[matching][0], bounds[matching][1], resolution)

QPEGalaxies["redshift_bin"] = pd.cut(QPEGalaxies["redshift"], bins=redshift_bins)
QPEGalaxies["m_star_bin"] = pd.cut(QPEGalaxies[matching], bins=m_star_bins)

# Sort QPEs into their bins
QPE_grid = np.zeros((resolution,resolution))
QPE_bins = []
for i in range(len(QPEGalaxies)):
    m_star, z = QPEGalaxies[matching][i], QPEGalaxies["redshift"][i]
    for mbin in list(m_star_bins)[::-1]:
        if m_star > mbin:
            m_star_index = list(m_star_bins).index(mbin)
            print("QPE mass:", m_star, mbin)
            break
    for zbin in list(redshift_bins)[::-1]:
        if z > zbin:
            z_index = list(redshift_bins).index(zbin)
            print("QPE z:", z, zbin)
            break
    #print(m_star_index, z_index)
    QPE_grid[m_star_index,z_index] += 1
    QPE_bins.append((m_star_index,z_index))

QPE_grid = QPE_grid/np.sum(QPE_grid) # normalize
assert np.sum(QPE_grid) > 0.999999 and np.sum(QPE_grid) < 1.000001


# Sort the reference galaxies in their respective bins:
ref_grid = np.zeros((resolution,resolution))
ref_grid_list = [[[] for i in range(resolution)] for i in range(resolution)]
for i in tqdm(range(len(refCat))):
    m_star, z = refCat[matching][i], refCat["redshift"][i]
    for mbin in list(m_star_bins)[::-1]:
        if m_star > mbin:
            m_star_index = list(m_star_bins).index(mbin)
            break
    for zbin in list(redshift_bins)[::-1]:
        if z > zbin:
            z_index = list(redshift_bins).index(zbin)
            break
    ref_grid[m_star_index,z_index] += 1
    ref_grid_list[m_star_index][z_index].append(i)

# Normalize
ref_grid = ref_grid/np.sum(ref_grid)

# Compare the reference catalog bin to the 2D probability density
fig = plt.gcf()
fig.set_size_inches((9.4,4.5))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
ax1.imshow(QPE_grid, origin="lower", cmap="Greys")
ax1.set_xlabel("$z$ bin", fontsize=15)
ax1.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax1.set_title("QPE sample probability density", fontsize=16)

ax2.imshow(ref_grid, origin="lower", cmap="Greys")
ax2.set_xlabel("$z$ bin", fontsize=15)
ax2.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax2.set_title("Comparison sample probability density", fontsize=16)
plt.subplots_adjust(wspace=0.35)
plt.show()



refCat = np.loadtxt("referenceCatalog_matching.txt")
# Calculate how many of the QPE bins contain at least one galaxy and make a list of them:
non_empty_bins = []
for bin in QPE_bins:
    if ref_grid[bin] > 0:
        non_empty_bins.append(bin)
    else:
        print(bin, "empty!")
print(len(non_empty_bins))

# Resample the reference catalog galaxies in the QPE bins:
goodGalaxies = []
for i in tqdm(range(sample_size)):
    bin_to_pick_from = non_empty_bins[int((i/sample_size)*len(non_empty_bins))]
    random_index_in_bin = np.random.randint(len(ref_grid_list[bin_to_pick_from[0]][bin_to_pick_from[1]]))
    random_galaxy = ref_grid_list[bin_to_pick_from[0]][bin_to_pick_from[1]][random_index_in_bin]
    goodGalaxies.append(refCat[random_galaxy,:])


goodGalaxies = np.array(goodGalaxies)[:48000,:]
print(goodGalaxies.shape)
np.savetxt("referenceCatalog_final_9bins.txt", goodGalaxies)
refCat = np.loadtxt("referenceCatalog_final_9bins.txt")
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
fieldnames[1] = "redshift"
fieldnames[63] = "m_star"
fieldnames[67] = "m_bh"
refCat = pd.read_csv("referenceCatalog_final_9bins.txt", delimiter=" ", header=None, names=fieldnames)

ref_grid = np.zeros((resolution,resolution))
for i in tqdm(range(len(refCat))):
    m_star, z = refCat[matching][i], refCat["redshift"][i]
    for mbin in list(m_star_bins)[::-1]:
        if m_star > mbin:
            m_star_index = list(m_star_bins).index(mbin)
            break
    for zbin in list(redshift_bins)[::-1]:
        if z > zbin:
            z_index = list(redshift_bins).index(zbin)
            break
    ref_grid[m_star_index,z_index] += 1
    ref_grid_list[m_star_index][z_index].append(i)

# Normalize
ref_grid = ref_grid/np.sum(ref_grid)

# Compare the 2D probability density to the resampled reference catalogue
fig = plt.gcf()
fig.set_size_inches((10,4.5))
ax1, ax2 = plt.subplot(121), plt.subplot(122)
ax1.imshow(QPE_grid, origin="lower", cmap="Greys")
ax1.set_xlabel("$z$ bin", fontsize=15)
ax1.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax1.set_title("QPE sample probability density", fontsize=16)

ax2.imshow(ref_grid, origin="lower", cmap="Greys")
ax2.set_xlabel("$z$ bin", fontsize=15)
ax2.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax2.set_title("Resampled comparison sample probability density", fontsize=16)
plt.subplots_adjust(wspace=0.35)
plt.show()