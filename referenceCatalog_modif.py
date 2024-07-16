import numpy as np

# Just delete unphysical BH masses and add a stellar surface density column


def cut_from_catalog(catalog, index, bounds, verbose=False):
    """
    catalog: numpy array
    index: int of index of parameter column which is under study here
    bounds: tuple of bound which we want to keep

    returns: numpy array with only remaining objects
    """
    catalog = np.array(catalog)
    if bounds[0] == None:
        good_indices_lo = catalog[:,index] == catalog[:,index]
    else:
        good_indices_lo = catalog[:,index] >= bounds[0]
    if bounds[1] == None:
        good_indices_hi = catalog[:,index] == catalog[:,index]
    else:
        good_indices_hi = catalog[:,index] <= bounds[1]
    good_indices = []
    for i in range(len(good_indices_lo)):
        good_indices.append(good_indices_lo[i] and good_indices_hi[i])
    cut_catalog = catalog[good_indices]
    if verbose:
        print(f"\x1b[31m{catalog.shape[0]-cut_catalog.shape[0]} objects cut\x1b[0m")
        print(f"\x1b[32m{cut_catalog.shape[0]} objects remaining\x1b[0m")
    return cut_catalog

reference_catalog = np.loadtxt("referenceCatalog.txt")
print("Physical black hole mass cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(0, 20), verbose=True)

smsd = np.log10((10**reference_catalog[:,63])/reference_catalog[:,59]**2)
reference_catalog = np.vstack((reference_catalog.T, smsd)).T
print("Physical SMSD cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=68, bounds=(0, 20), verbose=True)

print("Physical (B/T)g cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=12, bounds=(0, 1), verbose=True)

# redshift cut
print("Redshift cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=1, bounds=(0.01, 0.09), verbose=True)



# add Stellar Masses from Mendel2014
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table4.dat")
# Add a stellar mass column by fetching values from the Mendel2014 catalog
sorter = np.argsort(mendel2014[:,0])
indices = np.searchsorted(mendel2014[:,0], reference_catalog[:,0], sorter=sorter)
stellarMasses = mendel2014[:,2][indices]
reference_catalog = np.vstack((reference_catalog.T, np.array(stellarMasses))).T


np.savetxt("referenceCatalog_modif2.txt", reference_catalog)
