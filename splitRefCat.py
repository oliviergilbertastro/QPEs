"""
Code to split referenceCatalog.txt into smaller files so github doesn't whine.
We put referenceCatalog.txt into gitignore and simply rebuild it from the smaller files
"""

import numpy as np
from tqdm import tqdm

# The reference catalog is around 200mB while the largest file size is 100mB (and git recommends <50mB).

def split_catalog(catalog, number_of_children=4, savename="referenceCatalog"):
    length_cat = catalog.shape[0]
    for i in range(number_of_children):
        refCat = catalog[int(i*length_cat/number_of_children):int((i+1)*length_cat/number_of_children),:]
        print(int(i*length_cat/number_of_children),int((i+1)*length_cat/number_of_children))
        np.savetxt(f"{savename}_{i}.txt", refCat)
    pass

def recombine_catalog(catalogs, savename="referenceCatalog"):
    refCat = catalogs[0]
    for cat in catalogs[1:]:
        refCat = np.vstack((refCat,cat))
    print(refCat.shape)
    if input("Type 'YES' if you are absolutely sure you want to erase referenceCatalog.txt. This cannot be undone.") == "YES":
        np.savetxt(f"{savename}.txt", refCat)
    pass


if __name__ == "__main__":
    if input("Split reference catalog UNCERTAINTIES? [y/n]") == "y":
        refCat = np.loadtxt("referenceCatalog_with_uncertainties.txt")
        split_catalog(refCat, number_of_children=9, savename="referenceCatalog_with_uncertainties")

    if input("Recombine reference catalogs UNCERTAINTIES? [y/n]") == "y":
        refCats = []
        for i in range(9):
            refCats.append(np.loadtxt(f"referenceCatalog_with_uncertainties_{i}.txt"))
        recombine_catalog(refCats, savename="referenceCatalog_with_uncertainties")


    if input("Split reference catalog? [y/n]") == "y":
        refCat = np.loadtxt("referenceCatalog.txt")
        split_catalog(refCat, number_of_children=9)

    if input("Recombine reference catalogs? [y/n]") == "y":
        refCats = []
        for i in range(9):
            refCats.append(np.loadtxt(f"referenceCatalog_{i}.txt"))
        recombine_catalog(refCats, savename="referenceCatalog")

    #refCat1 = np.loadtxt("referenceCatalog.txt")
    #refCat2 = np.loadtxt("REFCAAT.txt")
    #for col in tqdm(range(refCat1.shape[1])):
    #    for line in range(refCat1.shape[0]):
    #        if refCat1[line,col] != refCat2[line,col]:
    #            print(line, col, ":", refCat1[line,col], refCat2[line,col])