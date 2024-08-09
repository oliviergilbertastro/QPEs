"""
Find velocity dispersions of the 30 TDE hosts in the MPA-JHU catalog to calculate mBH
"""

from download_data import hammerstein_TDE_coords, hammerstein_TDE_names
from astropy.io import fits
import numpy as np
from tqdm import tqdm
from utils import cut_from_catalog, get_smallest_sep_v2
from paper_data import hammerstein_TDE_redshifts, hammerstein_TDE_mBH
from ned_wright_cosmology import calculate_cosmo

mpajhu = fits.open("data/catalogs/MPA_JHU/galSpecInfo-dr8.fits")
ra_decs = fits.open("data/catalogs/asu.fit")
simard2011a = np.loadtxt("data/catalogs/Simard2011/table1.dat")
#simard2011b = np.loadtxt("data/catalogs/Simard2011/table2.dat")
simard2011c = np.loadtxt("data/catalogs/Simard2011/table3.dat")


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

superSimard = np.vstack((simard2011a.T, simard2011c.T)).T
print("This shape:", superSimard.shape)
superSimard = np.vstack((superSimard.T, simard_ra_decs.T)).T
print(superSimard.shape)
#1123718 lines by 98 columns, the last two are ra and dec, some columns are repeated (e.g. SPECOBJID), but I don't care

# Find the sersic indices of the hammerstein galaxies
import copy
remaining_hammerstein_TDE_coords = copy.deepcopy(hammerstein_TDE_coords)
remaining_hammerstein_TDE_coords.pop(23)
remaining_hammerstein_TDE_coords.pop(11)
remaining_hammerstein_TDE_coords.pop(6)
remaining_hammerstein_TDE_coords.pop(5)
hammerstein_in_Simard = np.array(remaining_hammerstein_TDE_coords)
from download_data import remaining_hammerstein_TDE_names
sim_remaining_hammerstein_TDE_names = []
sim_seps = []
sim_sersics = []
for i in tqdm(range(len(remaining_hammerstein_TDE_coords))):
    index, smallest_sep = get_smallest_sep_v2(remaining_hammerstein_TDE_coords[i], superSimard[:,96], superSimard[:,97])
    print(f"{remaining_hammerstein_TDE_coords[i]} vs {superSimard[index,96:98]} -> {smallest_sep < 3}")
    sim_seps.append(smallest_sep)
    sim_sersics.append(superSimard[index,93])
    if smallest_sep < 3:
        sim_remaining_hammerstein_TDE_names.append(remaining_hammerstein_TDE_names[i])
sim_seps = np.array(sim_seps)
sim_sersics = np.array(sim_sersics)
hammerstein_in_Simard = np.vstack((hammerstein_in_Simard.T, sim_sersics)).T
hammerstein_in_Simard = np.vstack((hammerstein_in_Simard.T, sim_seps)).T # separation between catalog and my coords in arcsec
print("Separation cut...")
hammerstein_in_Simard = cut_from_catalog(hammerstein_in_Simard, index=-1, bounds=(None, 3), verbose=True)

print(hammerstein_in_Simard)
print(sim_remaining_hammerstein_TDE_names)
import sys
sys.exit()

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
for i in [23,11,6,5]:
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
sigma_e = our_galaxies[:,6]/(10**(log_ratio_sigmas)) # in km/s
mBH = np.log10(10**9*(0.309)*(sigma_e/200)**4.38) # in solar masses
our_galaxies = np.vstack((our_galaxies.T, mBH)).T # velocity dispersion: 8


print(our_galaxies)
for i, index in enumerate(our_galaxies[:,0]):
    print(hammerstein_TDE_names[int(index)], our_galaxies[i,8], hammerstein_TDE_mBH[int(index)])