import numpy as np
from astropy.io import fits
from tqdm import tqdm

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
mpajhu = fits.open("data/catalogs/MPA_JHU/galSpecInfo-dr8.fits")
ra_decs = fits.open("data/catalogs/asu.fit")
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

from astropy.table import Table, Column
my_table = Table(ra_decs[1].data)
# Extract the names of the columns
colnames = my_table.colnames
print(colnames)
simard_ra_decs = []
for i in tqdm(range(len(my_table[colnames[0]]))):
    try:
        SpecObjID = int(my_table["objID"][i])
        simard_ra_decs.append([SpecObjID, my_table["_RAJ2000"][i], my_table["_DEJ2000"][i]])
    except:
        simard_ra_decs.append([SpecObjID, my_table["_RA"][i], my_table["_DE"][i]])
simard_ra_decs = np.array(simard_ra_decs)





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
Sigma_hlg = []
for i in range(len(simard2011a[:,0])):
    bulge_gr.append(simard2011a[i,46]-simard2011a[i,52])
    Sigma_hlg.append((simard2011c[i,5]+0.75254)/(np.pi*(simard2011c[i,14]/simard2011c[i,3])**2)) # this formula comes from section 4.6 of https://arxiv.org/pdf/1707.01559
reference_catalog = np.vstack((reference_catalog.T, np.array(bulge_gr))).T
reference_catalog = np.vstack((reference_catalog.T, np.array(Sigma_hlg))).T

# Start cutting the reference catalog:

# redshift cut
print("Redshift cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=1, bounds=(0.01, None), verbose=True)

# bulge g-r cut
print("Bulge g-r cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=61, bounds=(-50, 0.51), verbose=True)

# Sigma_hl,g cut
print("Sigma_hl,g cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=62, bounds=(2.05, None), verbose=True)


# Add a stellar mass column by fetching values from the Mendel2014 catalog
sorter = np.argsort(mendel2014[:,0])
indices = np.searchsorted(mendel2014[:,0], reference_catalog[:,0], sorter=sorter)
stellarMasses = mendel2014[:,2][indices]
reference_catalog = np.vstack((reference_catalog.T, np.array(stellarMasses))).T


# Add RA and DECs by fetching values from the online Simard Vizier tool:

sorter = np.argsort(simard_ra_decs[:,0])
indices = np.searchsorted(simard_ra_decs[:,0], reference_catalog[:,0], sorter=sorter)
ra_s, dec_s = simard_ra_decs[:,1][indices], simard_ra_decs[:,2][indices]
reference_catalog = np.vstack((reference_catalog.T, np.array(ra_s))).T
reference_catalog = np.vstack((reference_catalog.T, np.array(dec_s))).T
# Add the velocity dispersions in the reference catalog:
from utils import get_smallest_sep, get_smallest_sep_v2

# Latitude cut to be physical
print("Latitude cut in the MPA-JHU catalog...")
mpajhu = cut_from_catalog(mpajhu, index=2, bounds=(-90, 90), verbose=False)

sigma_a = []
seps = []
#reference_catalog = reference_catalog[:100,:] # make it small so it doesn't take an hour to test
for i in tqdm(range(len(reference_catalog[:,0]))):
    index, smallest_sep = get_smallest_sep_v2(reference_catalog[i,64:66], mpajhu[:,1], mpajhu[:,2])
    sigma_a.append(mpajhu[index,3])
    seps.append(smallest_sep)
sigma_a = np.array(sigma_a)
seps = np.array(seps)
reference_catalog = np.vstack((reference_catalog.T, np.array(sigma_a))).T
# Add the separation between the catalogs in the reference catalogs, make a cut if the sep is > 3", then remove this column
reference_catalog = np.vstack((reference_catalog.T, np.array(seps))).T
print("Separation cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(None, 3), verbose=True)

reference_catalog = reference_catalog[:,:-1]


# Calculate the black hole mass
log_ratio_sigmas = -0.065*np.log10(1.5/(reference_catalog[:,21]/reference_catalog[:,3]))-0.013*np.log10(1.5/(reference_catalog[:,21]/reference_catalog[:,3]))**2
sigma_e = reference_catalog[:,66]/(10**(log_ratio_sigmas)) # in km/s
mBH = np.log10(10**9*(0.309)*(sigma_e/200)) # in solar masses
reference_catalog = np.vstack((reference_catalog.T, mBH)).T


if False: # let's not do a black hole mass cut for the moment
    # Black hole mass cut
    print("Black hole mass cut...")
    reference_catalog = cut_from_catalog(reference_catalog, index=67, bounds=(5.5, 7), verbose=True)

np.savetxt(f"referenceCatalog.txt", reference_catalog)

import sys
sys.exit()

# Now we calculate black hole masses for the MPA_JHU catalog from its velocity dispersions
sigma_e = []
mBH = []
for i in range(len(mpajhu[:,3])):
    pass




example_catalog = np.array([[10, 13, 14],
                            [11, 12, 14],
                            [18, 11, 10],
                            [8, 9, 7]])

print(cut_from_catalog(example_catalog, 2, bounds=(3,6)))
print()