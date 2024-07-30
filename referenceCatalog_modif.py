"""
This code is QUICK to run and serves to quickly test different cuts from a larger reference catalog.
Ideally, a complete uncut catalog would be made beforehand and cuts would only be made in here, but it
would take days to run referenceCatalog.py, so cuts are pre-made, but I am trying to make less invasive
cuts to continually add more and more data to the reference catalog.
"""


import numpy as np

# Just delete unphysical BH masses and add a stellar surface density column

from utils import cut_from_catalog

reference_catalog = np.loadtxt("referenceCatalog.txt")
print("Physical black hole mass cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(0, 20), verbose=True)

ssmd = np.log10((10**reference_catalog[:,63])/reference_catalog[:,59]**2)
reference_catalog = np.vstack((reference_catalog.T, ssmd)).T
print("Physical SSMD cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=68, bounds=(0, 20), verbose=True)

print("Physical (B/T)g cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=12, bounds=(0, 1), verbose=True)

# Stricter bulge g-r cut (changed it from 0.51 to 1 in referenceCatalog.py)
print("Bulge g-r cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=61, bounds=(-50, 1), verbose=True)

#print("Black hole mass cut...")
#reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(5.5, 7), verbose=True)
#reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(None, 8), verbose=True)

print("Stellar mass cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=63, bounds=(9.5, 10.4), verbose=True)

print("r50 cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=59, bounds=(0, 100), verbose=True)


# redshift cut
print("Redshift cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=1, bounds=(0.01, 0.1), verbose=True)


#print("Sigma_hl,g cut...")
#reference_catalog = cut_from_catalog(reference_catalog, index=62, bounds=(2.05, None), verbose=True)

if False:
    # add Stellar Masses from Mendel2014
    mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table4.dat")
    # Add a stellar mass column by fetching values from the Mendel2014 catalog
    sorter = np.argsort(mendel2014[:,0])
    indices = np.searchsorted(mendel2014[:,0], reference_catalog[:,0], sorter=sorter)
    stellarMasses = mendel2014[:,2][indices]
    reference_catalog = np.vstack((reference_catalog.T, np.array(stellarMasses))).T


np.savetxt("referenceCatalog_modif2.txt", reference_catalog[::2])
