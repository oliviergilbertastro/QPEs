"""
Trying to sample from a 2D distribution in order to match a second 2D distribution, will be fun!
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from paper_data import QPE_redshifts

# Let's load the background galaxies
refCat = np.loadtxt("referenceCatalog_modif2.txt")
fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
fieldnames[1] = "redshift"
fieldnames[63] = "m_star"
refCat = pd.read_csv("referenceCatalog_modif2.txt", delimiter=" ", header=None, names=fieldnames)

print(refCat)
# Let's load the QPE hosts:
fieldnames = ["bt_ratio", "n_sersic", "smsd", "m_star", "m_bh", "redshift"]
QPEGalaxies = pd.read_csv("QPE_distribution.txt", header=None, names=fieldnames, delimiter=" ")
print(QPEGalaxies)


# Choose bounds for that z-M_star grid
bounds = {"redshift": (0.01,0.2),
          "m_star": (9,11)}
resolution = 30

redshift_bins = np.linspace(bounds["redshift"][0], bounds["redshift"][1], resolution)
m_star_bins = np.linspace(bounds["m_star"][0], bounds["m_star"][1], resolution)

QPEGalaxies["redshift_bin"] = pd.cut(QPEGalaxies["redshift"], bins=redshift_bins)
QPEGalaxies["m_star_bin"] = pd.cut(QPEGalaxies["m_star"], bins=m_star_bins)

bin_counts = QPEGalaxies.groupby(['redshift_bin', 'm_star_bin']).size()

print(QPEGalaxies["redshift_bin"])
print(QPEGalaxies["m_star_bin"])
print(bin_counts)
# Make that weird grid