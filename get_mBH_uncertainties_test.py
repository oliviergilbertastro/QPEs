import numpy as np
from astropy.io import fits
from tqdm import tqdm
import sys
from utils import cut_from_catalog, mergeCatalogs_withObjIDs, get_smallest_sep_v2

# load all catalog text files
mpajhu = fits.open("data/catalogs/MPA_JHU/galSpecInfo-dr8.fits")
ra_decs = fits.open("data/catalogs/asu.fit")
simard2011a = np.loadtxt("data/catalogs/Simard2011/table1.dat")
#simard2011b = np.loadtxt("data/catalogs/Simard2011/table2.dat")
simard2011c = np.loadtxt("data/catalogs/Simard2011/table3.dat")
print(simard2011a.shape)
print(simard2011c.shape)
# Extract the data in an astropy table
from astropy.table import Table, Column
my_table = Table(mpajhu[1].data)
# Extract the names of the columns
colnames = my_table.colnames
mpajhu = []
print(colnames)
for i in tqdm(range(len(my_table[colnames[0]]))):
    try:
        SpecObjID = int(my_table["SPECOBJID"][i])
        mpajhu.append([SpecObjID, my_table["RA"][i], my_table["DEC"][i], my_table["V_DISP"][i], my_table["V_DISP_ERR"][i]])
    except:
        pass
mpajhu = np.array(mpajhu)
# Latitude cut to be physical
print("Latitude cut in the MPA-JHU catalog...")
mpajhu = cut_from_catalog(mpajhu, index=2, bounds=(-90, 90), verbose=False)

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








# Create a reference catalog, which is just the simard catalogs with ra and decs
print(simard_ra_decs.shape)
reference_catalog = np.vstack((simard_ra_decs.T, simard2011a.T)).T
print(reference_catalog.shape)
reference_catalog = np.vstack((reference_catalog.T, simard2011c.T)).T
print(reference_catalog.shape)


from download_data import *
gal_names = np.concatenate((QPE_names, TDE_names, french_TDE_names))
gal_coords = np.concatenate((QPE_coords, TDE_coords, french_TDE_coords))

#reference_catalog = reference_catalog[:100,:] # make it small so it doesn't take an hour to test
print("Cross-matching RA&DEC...")
simard_specObjIDs = []
seps = []
simard_indices = []
simard_gal_names = []
for i in tqdm(range(len(gal_names))):
    index, smallest_sep = get_smallest_sep_v2(gal_coords[i], reference_catalog[:,1], reference_catalog[:,2])
    seps.append(smallest_sep)
    if smallest_sep < 3:
        simard_specObjIDs.append(reference_catalog[index,0])
        simard_indices.append(index)
        simard_gal_names.append(gal_names[i])
seps = np.array(seps)
# Add the separation between the catalogs in the reference catalogs, make a cut if the sep is > 3", then remove this column

for i in range(len(simard_gal_names)):
    print(i, simard_specObjIDs[i], "vs", reference_catalog[simard_indices[i],0])

print("Simard Separation cut...")
reference_catalog = reference_catalog[simard_indices]
print(reference_catalog.shape)

for i in range(len(simard_gal_names)):
    print(i, simard_specObjIDs[i], "vs", reference_catalog[i,0])


# Add the mpa-jhu catalog
print("Cross-matching RA&DEC...")
mpa_specobjid = []
mpa_vdisp = []
mpa_vdisp_sigma = []
seps = []
mpa_indices = []
mpa_gal_names = []
for i in tqdm(range(len(simard_gal_names))):
    index, smallest_sep = get_smallest_sep_v2(reference_catalog[i,1:3], mpajhu[:,1], mpajhu[:,2])
    seps.append(smallest_sep)
    if smallest_sep < 3:
        mpa_specobjid.append(mpajhu[index,0])
        mpa_vdisp.append(mpajhu[index,3])
        mpa_vdisp_sigma.append(mpajhu[index,4])
        mpa_indices.append(index)
        mpa_gal_names.append(simard_gal_names[i])

print(mpa_specobjid)

print("MPA-JHU Separation cut...")
mpajhu = mpajhu[mpa_indices]
print(mpajhu.shape)

simard_and_mpa = []
for i in range(len(simard_gal_names)):
    if simard_gal_names[i] in mpa_gal_names:
        simard_and_mpa.append(True)
    else:
        simard_and_mpa.append(False)

reference_catalog = reference_catalog[simard_and_mpa]
for i in range(len(mpa_gal_names)):
    print(mpa_gal_names[i], end=" ")
    print(reference_catalog[i,1:3]-mpajhu[i,1:3]) # distance in ra and dec (degrees)

# Now we have the 13 remaining galaxies with all the information necessary to calculate m_BH uncertainties.
# Calculate the black hole mass
R_ap = 1.5 # Aperture radius [arcsec]
spheroid_radius = reference_catalog[:,24] # Rchl_r [kpc]
scale = reference_catalog[:,6] # Scale [kpc/arcsec]
v_disp = np.array(mpa_vdisp)
v_disp_sig = np.array(mpa_vdisp_sigma)

def get_mBH(spheroid_radius,scale,v_disp,v_disp_sig,R_ap=1.5):
    '''Kormendy relation'''
    R_e = (spheroid_radius/scale)
    log_ratio_sigmas = -0.065*np.log10(R_ap/R_e)-0.013*(np.log10(R_ap/R_e))**2
    # Still no uncertainty yet (none listed in Simard catalog at least)
    sigma_e = v_disp/(10**(log_ratio_sigmas)) # in km/s
    sigma_e_sig = v_disp_sig/(10**(log_ratio_sigmas)) # in km/s
    f1 = sigma_e/200
    sig_f1 = sigma_e_sig/200
    f2 = f1**4.38
    sig_f2 = (f2*4.38*sig_f1)/f1
    f3 = 10**9*(0.309)*f2
    sig_f3 = 10**9*(0.309)*sig_f2
    f4 = np.log10(f3)
    sig_f4 = np.abs(sig_f3/(f3*np.log(10)))
    #mBH = np.log10(10**9*(0.309)*(sigma_e/200)**4.38) # in solar masses
    return f4, sig_f4

mBH, mBH_sigma = get_mBH(spheroid_radius,scale,v_disp,v_disp_sig)

for i in range(len(mpa_gal_names)):
    print(f"{mpa_gal_names[i]}  mBH: {mBH[i]} +/- {mBH_sigma[i]}")

sys.exit()

# Add a bulge g-r column to the reference catalog
# Add a Sigma_hl,g column to the reference catalog
bulge_gr = []
Sigma_hlg = []
for i in range(len(reference_catalog[:,0])):
    bulge_gr.append(reference_catalog[i,46]-reference_catalog[i,52])
    Sigma_hlg.append((simard2011c[i,5]+0.75254)/(np.pi*(simard2011c[i,14]/simard2011c[i,3])**2)) # this formula comes from section 4.6 of https://arxiv.org/pdf/1707.01559
reference_catalog = np.vstack((reference_catalog.T, np.array(bulge_gr))).T
reference_catalog = np.vstack((reference_catalog.T, np.array(Sigma_hlg))).T

# Start cutting the reference catalog:


# Add RA and DECs by fetching values from the online Simard Vizier tool:

print(reference_catalog.shape)
reference_catalog = mergeCatalogs_withObjIDs(reference_catalog, simard_ra_decs, columnsToAdd=[1,2])
print(reference_catalog.shape)
# Add the velocity dispersions in the reference catalog:

# Latitude cut to be physical
print("Latitude cut in the MPA-JHU catalog...")
mpajhu = cut_from_catalog(mpajhu, index=2, bounds=(-90, 90), verbose=False)

sigma_a = []
seps = []
#reference_catalog = reference_catalog[:100,:] # make it small so it doesn't take an hour to test
print("Cross-matching RA&DEC...")
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
R_e = (reference_catalog[:,21]/reference_catalog[:,3])
log_ratio_sigmas = -0.065*np.log10(1.5/R_e)-0.013*(np.log10(1.5/R_e))**2
sigma_e = reference_catalog[:,66]/(10**(log_ratio_sigmas)) # in km/s
mBH = np.log10(10**9*(0.309)*(sigma_e/200)**4.38) # in solar masses
reference_catalog = np.vstack((reference_catalog.T, mBH)).T
