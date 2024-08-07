"""
Find velocity dispersions of the 30 TDE hosts in the MPA-JHU catalog to calculate mBH
"""

from download_data import hammerstein_TDE_coords, hammerstein_TDE_names
from astropy.io import fits
import numpy as np
from tqdm import tqdm
from utils import cut_from_catalog, get_smallest_sep_v2
from paper_data import hammerstein_TDE_redshifts
from ned_wright_cosmology import calculate_cosmo

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


hammerstein_properties = np.loadtxt("hammerstein_TDE_distribution.txt")

our_galaxies = np.array(hammerstein_TDE_coords) # coordinates ra: 1, dec: 2
our_galaxies = np.vstack((range(len(hammerstein_TDE_coords)), our_galaxies.T)).T # index: 0
our_galaxies = np.vstack((our_galaxies.T, hammerstein_TDE_redshifts)).T # redshifts: 3
for i in [5,6,11,23]:
    our_galaxies = np.vstack((our_galaxies[:i],our_galaxies[i+1:]))
our_galaxies = np.vstack((our_galaxies.T, hammerstein_properties[:,0])).T # sersic indices: 4
scales = []
for i in range(len(our_galaxies)):
    scales.append(calculate_cosmo(our_galaxies[i,3])["kpc_DA"])
our_galaxies = np.vstack((our_galaxies.T, scales)).T # scale: 5

print(our_galaxies)

# Latitude cut to be physical
print("Latitude cut in the MPA-JHU catalog...")
mpajhu = cut_from_catalog(mpajhu, index=2, bounds=(-90, 90), verbose=True)

sigma_a = []
seps = []
#reference_catalog = reference_catalog[:100,:] # make it small so it doesn't take an hour to test
print("Cross-matching RA&DEC...")
for i in tqdm(range(len(our_galaxies))):
    index, smallest_sep = get_smallest_sep_v2(our_galaxies[i,1:3], mpajhu[:,1], mpajhu[:,2])
    print(f"{our_galaxies[i,1:3]} vs {mpajhu[index,1:3]} -> {smallest_sep < 3}")
    sigma_a.append(mpajhu[index,3])
    seps.append(smallest_sep)
sigma_a = np.array(sigma_a)
seps = np.array(seps)
our_galaxies = np.vstack((our_galaxies.T, np.array(sigma_a))).T # velocity dispersion: 6
our_galaxies = np.vstack((our_galaxies.T, np.array(seps))).T # separation between catalog and my coords in arcsec: 7
print("Separation cut...")
our_galaxies = cut_from_catalog(our_galaxies, index=-1, bounds=(None, 3), verbose=True)

# Calculate the black hole mass
R_e = (our_galaxies[:,4]/our_galaxies[:,5])
log_ratio_sigmas = -0.065*np.log10(1.5/R_e)-0.013*(np.log10(1.5/R_e))**2
sigma_e = our_galaxies[:,3]/(10**(log_ratio_sigmas)) # in km/s
mBH = np.log10(10**9*(0.309)*(sigma_e/200)**4.38) # in solar masses
our_galaxies = np.vstack((our_galaxies.T, mBH)).T


print(our_galaxies)
for i in our_galaxies[:,0]:
    print(hammerstein_TDE_names[int(i)], our_galaxies[int(i)])