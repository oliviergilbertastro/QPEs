import numpy as np
from astropy.io import fits
from tqdm import tqdm
import sys
from utils import cut_from_catalog, mergeCatalogs_withObjIDs, get_smallest_sep_v2

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
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

# redshift cut
print("Redshift cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=4, bounds=(0.01, 0.1), verbose=True)

#reference_catalog = reference_catalog[:1000] # uncomment to make tests
# Add the mpa-jhu catalog
print("Cross-matching RA&DEC...")
mpa_vdisp = []
mpa_vdisp_sigma = []
seps = []
for i in tqdm(range(reference_catalog.shape[0])):
    index, smallest_sep = get_smallest_sep_v2(reference_catalog[i,1:3], mpajhu[:,1], mpajhu[:,2])
    seps.append(smallest_sep)
    mpa_vdisp.append(mpajhu[index,3])
    mpa_vdisp_sigma.append(mpajhu[index,4])

print("Adding velocity dispersion...")
reference_catalog = np.vstack((reference_catalog.T, mpa_vdisp)).T
print(reference_catalog.shape)
print("Adding velocity dispersion error...")
reference_catalog = np.vstack((reference_catalog.T, mpa_vdisp_sigma)).T
print(reference_catalog.shape)
print("Adding separation...")
reference_catalog = np.vstack((reference_catalog.T, seps)).T
print(reference_catalog.shape)

print("Separation cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=-1, bounds=(None, 3), verbose=True)

print("Removing separation column...")
reference_catalog = reference_catalog[:,:-1]
print(reference_catalog.shape)



# Now we have the 13 remaining galaxies with all the information necessary to calculate m_BH uncertainties.
# Calculate the black hole mass
R_ap = 1.5 # Aperture radius [arcsec]
spheroid_radius = reference_catalog[:,24] # Rchl_r [kpc]
scale = reference_catalog[:,6] # Scale [kpc/arcsec]
v_disp = reference_catalog[:,98] # Velocity dispersion [km/s]
v_disp_sig = reference_catalog[:,99] # Velocity dispersion error [km/s]

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
print("Adding mBH...")
reference_catalog = np.vstack((reference_catalog.T, mBH)).T
print(reference_catalog.shape)
print("Adding mBH error...")
reference_catalog = np.vstack((reference_catalog.T, mBH_sigma)).T
print(reference_catalog.shape)

print("Physical black hole mass cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=100, bounds=(0, 20), verbose=True)

# Add a stellar mass column by fetching values from the Mendel2014 catalog
print("Adding stellar mass and its uncertainties columns")
print(reference_catalog.shape)
reference_catalog = mergeCatalogs_withObjIDs(reference_catalog, mendel2014, columnsToAdd=[2,3,4])
print(reference_catalog.shape)

reference_catalog[:,-2] = reference_catalog[:,-3]-reference_catalog[:,-2] # convert 16th percentiles stellar mass log to uncertainties
reference_catalog[:,-1] = reference_catalog[:,-1]-reference_catalog[:,-3] # convert 84th percentiles stellar mass log to uncertainties

print("Adding stellar mass surface density")
ssmd = np.log10((10**reference_catalog[:,102])/reference_catalog[:,77]**2)
reference_catalog = np.vstack((reference_catalog.T, ssmd)).T
print(reference_catalog.shape)
print("Physical SSMD cut...")
reference_catalog = cut_from_catalog(reference_catalog, index=105, bounds=(0, 20), verbose=True)


np.savetxt("referenceCatalog_with_uncertainties.txt", reference_catalog)

sys.exit()
