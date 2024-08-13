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
from utils import mergeCatalogs_withObjIDs, get_smallest_sep_v2, cut_from_array
from download_data import TDE_names, TDE_coords, french_TDE_names, french_TDE_coords, QPE_names, QPE_coords

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
ra_decs = fits.open("data/catalogs/asu.fit")

QPE_m_star = np.loadtxt("QPE_stellarMasses.txt")
TDE_m_star = np.loadtxt("TDE_stellarMasses.txt")
french_TDE_m_star = np.loadtxt("french_TDE_stellarMasses.txt")

from astropy.table import Table, Column
my_table = Table(ra_decs[1].data)
# Extract the names of the columns
colnames = my_table.colnames
simard_ra_decs = []
for i in tqdm(range(len(my_table[colnames[0]]))):
    try:
        SpecObjID = int(my_table["objID"][i])
        simard_ra_decs.append([SpecObjID, my_table["_RAJ2000"][i], my_table["_DEJ2000"][i]])
    except:
        simard_ra_decs.append([SpecObjID, my_table["_RA"][i], my_table["_DE"][i]])
simard_ra_decs = np.array(simard_ra_decs)

TDE_names = np.concatenate((QPE_names, TDE_names, french_TDE_names))
TDE_coords = np.concatenate((QPE_coords, TDE_coords, french_TDE_coords))
TDE_m_star = np.concatenate((QPE_m_star, TDE_m_star, french_TDE_m_star))[:,0]

assert len(TDE_names) == len(TDE_coords)



bad_indices = [] # indices to cut
TDE_names = cut_from_array(TDE_names, bad_indices)
TDE_coords = cut_from_array(TDE_coords, bad_indices)
TDE_m_star = cut_from_array(TDE_m_star, bad_indices)



sim_TDE_names = []
sim_seps = []
sim_objIDs = []
m_star_in_sim = []
for i in (range(len(TDE_names))):
    index, smallest_sep = get_smallest_sep_v2(TDE_coords[i], simard_ra_decs[:,1], simard_ra_decs[:,2])
    print(f"{TDE_names[i]}: {TDE_coords[i]} vs {simard_ra_decs[index,1:]} -> {smallest_sep < 3}")
    if smallest_sep < 3:
        sim_TDE_names.append(TDE_names[i])
        sim_seps.append(smallest_sep)
        sim_objIDs.append(simard_ra_decs[index,0])
        m_star_in_sim.append(TDE_m_star[i])
sim_seps = np.array(sim_seps)
sim_objIDs = np.array(sim_objIDs)

print(len(sim_objIDs))
print(len(m_star_in_sim))
myCat = np.array([sim_objIDs, m_star_in_sim]).T


print(myCat.shape)

myCat = mergeCatalogs_withObjIDs(myCat, mendel2014, columnsToAdd=[2,])
mendel_TDE_names = []
mendel_indices = []
for i in range(len(sim_objIDs)):
    if sim_objIDs[i] in myCat[:,0]:
        mendel_TDE_names.append(sim_TDE_names[i])
        print(f"\x1b[32m{sim_TDE_names[i]}\x1b[0m")#, int(sim_objIDs[i])-int(myCat[list(myCat[:,0]).index(sim_objIDs[i]),0]))
    else:
        print(f"\x1b[31m{sim_TDE_names[i]}\x1b[0m")
print(myCat.shape)
myCat[:,1] = np.log10(myCat[:,1])


def slope45(x, offset):
    return x+offset

def slope(x, slo, offset):
    return x*slo+offset

def slope_noOffset(x, slo):
    return x*slo

from scipy.optimize import curve_fit
res = curve_fit(slope45, myCat[:,1], myCat[:,2])[0]
res1 = curve_fit(slope, myCat[:,1], myCat[:,2])[0]
res2 = curve_fit(slope_noOffset, myCat[:,1], myCat[:,2])[0]

print("Fixed slope free offset fit:", res)
print("Free slope free offset fit:", res1)
print("Free slope fixed offset fit:", res2)

plt.plot(myCat[:,1], myCat[:,2], "o", color="blue", markersize=8)
min, max = np.min(np.concatenate((myCat[:,1], myCat[:,2]))), np.max(np.concatenate((myCat[:,1], myCat[:,2])))
plt.plot(np.linspace(min, max, 1000), np.linspace(min, max, 1000), "--", color="red", linewidth=2, markersize=7)
#plt.plot(np.linspace(min, max, 1000), np.linspace(min, max, 1000)+res[0], "--", color="blue", linewidth=2, markersize=7)
#plt.plot(np.linspace(min, max, 1000), np.linspace(min, max, 1000)*res1[0]+res1[1], "--", color="orange", linewidth=2, markersize=7)
#plt.plot(np.linspace(min, max, 1000), np.linspace(min, max, 1000)*res2[0], "--", color="purple", linewidth=2, markersize=7)
plt.ylabel("My prospector $\log(M_\star)$", fontsize=16)
plt.xlabel("Mendel 2014 $\log(M_\star)$", fontsize=16)
plt.show()