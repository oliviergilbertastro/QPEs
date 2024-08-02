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
resolution = 100
smoothing_param = resolution/6
sample_size=1000000
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

QPE_grid = QPE_grid/np.sum(QPE_grid) # normalize
assert np.sum(QPE_grid) > 0.999999 and np.sum(QPE_grid) < 1.000001
plt.imshow(QPE_grid, origin="lower", cmap="Greys")
plt.xlabel("$z$ bin", fontsize=15)
plt.ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
plt.title("Normalized probability density", fontsize=16)
plt.show()

# smooth it
QPE_grid = sp.ndimage.gaussian_filter(QPE_grid, sigma=smoothing_param, mode='constant')
QPE_grid = QPE_grid/np.sum(QPE_grid) # normalize
assert np.sum(QPE_grid) > 0.999999 and np.sum(QPE_grid) < 1.000001
plt.imshow(QPE_grid, origin="lower", cmap="Greys")
plt.xlabel("$z$ bin", fontsize=15)
plt.ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
plt.title("Normalized probability density", fontsize=16)
plt.show()

# Marginalize the probability so it's only dependent on the redshift z
marginalized_prob = np.sum(QPE_grid, axis=0)
assert np.sum(marginalized_prob) > 0.999999 and np.sum(marginalized_prob) < 1.000001

norm_cdf = np.cumsum(marginalized_prob)
ax1, ax2 = plt.subplot(121), plt.subplot(122)
ax1.plot(redshift_bins, marginalized_prob)
ax1.set_xlabel("$z$", fontsize=15)
ax1.set_ylabel("Marginalized probability density", fontsize=15)
ax2.plot(redshift_bins, norm_cdf)
ax2.set_xlabel("$z$", fontsize=15)
ax2.set_ylabel("CDF", fontsize=15)
plt.show()

def get_sampled_element(cdf, sample_size=1):
    """Returns index"""
    a_list = np.random.uniform(0, 1, sample_size)
    indices = []
    for a in a_list:
        indices.append(np.argmax(cdf>=a))
    return indices

sample_indices_z = get_sampled_element(norm_cdf, sample_size=sample_size)
if show_plots:
    # This is just to make sure I am sampling correctly
    #print(redshift_bins[sample_index])
    ax1 = plt.subplot(111)
    ax1.plot(redshift_bins, marginalized_prob, linewidth=2, label="PDF")
    ax1.set_xlabel("$z$", fontsize=15)
    ax1.set_ylabel("Probability density", fontsize=15)
    counts2, bins2 = np.histogram(redshift_bins[sample_indices_z], bins=len(redshift_bins))
    ax1.stairs(counts2/len(sample_indices_z), bins2, fill=True, label="Sampled from CDF")
    plt.legend()
    plt.show()

# Get CDF of m_star at the specific redshift bin:
def get_m_star_cdf(z_bin):
    return np.cumsum(QPE_grid[:,z_bin]/np.sum(QPE_grid[:,z_bin]))

# For each sampled redshift, sample a stellar mass
sample_indices_m_star = []
sampled_grid = np.zeros((resolution,resolution))
for i in tqdm(range(len(sample_indices_z))):
    if show_plots and i == 0:
        # This is just to make sure I am sampling correctly
        #print(redshift_bins[sample_index])

        m_star_CDF = get_m_star_cdf(sample_indices_z[i])
        sample_indices_star = get_sampled_element(m_star_CDF, sample_size=1000000)
        ax1 = plt.subplot(111)
        ax1.plot(m_star_bins, QPE_grid[:,sample_indices_z[i]]/np.sum(QPE_grid[:,sample_indices_z[i]]), linewidth=2, label="PDF")
        #ax1.plot(m_star_bins, m_star_CDF, linewidth=2, label="CDF")
        ax1.set_xlabel("$M_\star$" if matching == "m_star" else "$M_\mathrm{BH}$", fontsize=15)
        ax1.set_ylabel("Probability density", fontsize=15)
        counts2, bins2 = np.histogram(m_star_bins[sample_indices_star], bins=len(redshift_bins))
        ax1.stairs(counts2/len(sample_indices_star), bins2, fill=True, label="Sampled from CDF")
        plt.legend()
        plt.show()
    m_star_CDF = get_m_star_cdf(sample_indices_z[i])
    sample_index_m_star = get_sampled_element(m_star_CDF, sample_size=1)[0]
    sample_indices_m_star.append(sample_index_m_star)
    sampled_grid[sample_indices_m_star[i], sample_indices_z[i]] += 1

# Normalize the sampled grid:
sampled_grid = sampled_grid/np.sum(sampled_grid)


# Compare the sampled bins to the 2D probability density
fig = plt.gcf()
fig.set_size_inches((8,4.5))
ax1, ax2 = plt.subplot(121), plt.subplot(122)
ax1.imshow(QPE_grid, origin="lower", cmap="Greys")
ax1.set_xlabel("$z$ bin", fontsize=15)
ax1.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax1.set_title("QPE sample probability density", fontsize=16)

ax2.imshow(sampled_grid, origin="lower", cmap="Greys")
ax2.set_xlabel("$z$ bin", fontsize=15)
ax2.set_ylabel("$M_\star$ bin" if matching == "m_star" else "$M_\mathrm{BH}$ bin", fontsize=15)
ax2.set_title("Resampled probability density", fontsize=16)
plt.subplots_adjust(wspace=0.35)
plt.show()




# Sort the reference galaxies in their respective bins:
# Sort QPEs into their bins
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
ax1, ax2 = plt.subplot(121), plt.subplot(122)
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

old_skips = 0
new_skips = 0
goodGalaxies = []
refCat = np.loadtxt("referenceCatalog_matching.txt")
# Resample the reference catalog galaxies so they look gut.
for z_index, m_index in tqdm(zip(sample_indices_z, sample_indices_m_star)):
    try:
        random_index_in_bin = np.random.randint(len(ref_grid_list[m_index][z_index]))
        random_galaxy = ref_grid_list[m_index][z_index][random_index_in_bin]
        goodGalaxies.append(refCat[random_galaxy,:])
    except:
        done = False
        old_skips += 1
        for x in [-1,1]:
            for y in [-1,1]:
                if not done:
                    try:
                        random_index_in_bin = np.random.randint(len(ref_grid_list[m_index][z_index]))
                        random_galaxy = ref_grid_list[m_index+x][z_index+y][random_index_in_bin]
                        goodGalaxies.append(refCat[random_galaxy,:])
                        done = True
                    except:
                        pass
        new_skips += 1
        pass
print("Old skip percentage:", old_skips/sample_size)
print("New skip percentage:", new_skips/sample_size)
    #input("...")
goodGalaxies = np.array(goodGalaxies)[:48000,:]
print(goodGalaxies.shape)
np.savetxt("referenceCatalog_final.txt", goodGalaxies)
refCat = np.loadtxt("referenceCatalog_final.txt")
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
fieldnames[1] = "redshift"
fieldnames[63] = "m_star"
fieldnames[67] = "m_bh"
refCat = pd.read_csv("referenceCatalog_final.txt", delimiter=" ", header=None, names=fieldnames)

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