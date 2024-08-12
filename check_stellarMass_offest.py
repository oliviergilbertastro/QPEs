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

# load all catalog text files
mendel2014 = np.loadtxt("data/catalogs/Mendel2014/table6.dat")
ra_decs = fits.open("data/catalogs/asu.fit")