import numpy as np
from astropy.io import fits

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
mpajhu = fits.open("data/catalogs/MPA_JHU/galSpecInfo-dr8.fits")
simard2011a = np.loadtxt("data/catalogs/Simard2011/table1.dat")
#simard2011b = np.loadtxt("data/catalogs/Simard2011/table2.dat")
simard2011c = np.loadtxt("data/catalogs/Simard2011/table3.dat")
