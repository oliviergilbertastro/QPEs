import numpy as np
from astropy.io import fits
from tqdm import tqdm

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
mpajhu = fits.open("data/catalogs/MPA_JHU/galSpecInfo-dr8.fits")
simard2011a = np.loadtxt("data/catalogs/Simard2011/table1.dat")
#simard2011b = np.loadtxt("data/catalogs/Simard2011/table2.dat")
simard2011c = np.loadtxt("data/catalogs/Simard2011/table3.dat")

# Extract the data in an astropy table
from astropy.table import Table, Column
my_table = Table(mpajhu[1].data)

# Extract the names of the columns
colnames = my_table.colnames


mpajhu = []
for i in tqdm(range(len(my_table[colnames[0]]))):
    try:
        SpecObjID = int(my_table["SPECOBJID"][i])
        mpajhu.append([SpecObjID, my_table["RA"][i], my_table["DEC"][i], my_table["V_DISP"][i], my_table["V_DISP_ERR"][i]])
    except:
        pass
mpajhu = np.array(mpajhu)







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




# Create a reference catalog, which is just the simard catalog with sersic indices and r50s
reference_catalog = np.vstack((simard2011a.T, simard2011c[:,15])).T #r50
reference_catalog = np.vstack((reference_catalog.T, simard2011c[:,34])).T #n_sersic

# Add a bulge g-r column to the reference catalog
# Add a Sigma_hl,g column to the reference catalog
bulge_gr = []
sigma_hlg = []
for i in range(len(simard2011a[:,0])):
    bulge_gr.append(simard2011a[i,46]-simard2011a[i,52])
    sigma_hlg.append((simard2011c[i,5]+0.75254)/(np.pi*(simard2011c[i,14]/simard2011c[i,3])**2)) # this formula comes from section 4.6 of https://arxiv.org/pdf/1707.01559
reference_catalog = np.vstack((reference_catalog.T, np.array(bulge_gr))).T
reference_catalog = np.vstack((reference_catalog.T, np.array(sigma_hlg))).T

# Start cutting the reference catalogue:

# redshift cut
print("Redshift cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=1, bounds=(0.01, 0.09), verbose=True)

# bulge g-r cut
print("Bulge g-r cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=61, bounds=(-50, 0.51), verbose=True)

# Sigma_hl,g cut
print("Sigma_hl,g cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=62, bounds=(2.05, None), verbose=True)


# Add a stellar mass column by fetching values from the Mendel2014 catalogue
if False:
    def getSM(objID):
        """Returns the stellar mass from the Mendel2014 catalogue by inputting a PhotObjID"""
        diff_in_objID = np.abs(mendel2014[:,0] - np.ones_like(mendel2014[:,0])*objID)
        index = list(diff_in_objID).index(np.min(diff_in_objID))
        good_line = mendel2014[index]
        return good_line[2]
    stellarMasses = []
    for objID in tqdm(reference_catalog[:,0]):
        stellarMasses.append(getSM(objID))
sorter = np.argsort(mendel2014[:,0])
indices = np.searchsorted(mendel2014[:,0], reference_catalog[:,0], sorter=sorter)
print(indices)
stellarMasses = mendel2014[:,0][indices]
reference_catalog = np.vstack((reference_catalog.T, np.array(stellarMasses))).T

import sys
sys.exit()

example_catalog = np.array([[10, 13, 14],
                            [11, 12, 14],
                            [18, 11, 10],
                            [8, 9, 7]])

print(cut_from_catalog(example_catalog, 2, bounds=(3,6)))
print()