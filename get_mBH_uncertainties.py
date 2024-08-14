"""
Compare my calculated stellar masses to the stellar masses in the Mendel catalogue to see if there is a systematic offset.

Procedure outline:
1. Find the SPECOBJID of the galaxy by cross-matching RA and DEC in the Simard catalog ra and decs.
2. Cross-match the SPECOBJID to the Mendel2014 catalog.
3. ???
4. Profit
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm
from utils import mergeCatalogs_withObjIDs, get_smallest_sep_v2, cut_from_array, cut_from_catalog
from download_data import TDE_names, TDE_coords, french_TDE_names, french_TDE_coords, QPE_names, QPE_coords

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


gal_names = np.concatenate((QPE_names, TDE_names, french_TDE_names))
gal_coords = np.concatenate((QPE_coords, TDE_coords, french_TDE_coords))

assert len(TDE_names) == len(TDE_coords)

good_gal_names = []
gal_mBH = []
for i in (range(len(gal_names))):
    index, smallest_sep = get_smallest_sep_v2(gal_coords[i], mpajhu[:,1], mpajhu[:,2])
    print(f"{gal_names[i]}: {gal_coords[i]} vs {mpajhu[index,1:3]} -> {smallest_sep < 3}")
    if smallest_sep < 3:
        good_gal_names.append(gal_names[i])
        gal_mBH.append((mpajhu[i,3],mpajhu[i,4],mpajhu[i,4]))
gal_mBH = np.array(gal_mBH)



for i in range(len(gal_names)):
    if gal_names[i] in good_gal_names:
        print(f"\x1b[32m{gal_names[i]}\x1b[0m", gal_mBH[good_gal_names.index(gal_names[i])])
    else:
        print(f"\x1b[31m{gal_names[i]}\x1b[0m")







