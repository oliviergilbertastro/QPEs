import numpy as np
import astropy.io.fits as pyfits
import galight.tools.astro_tools as astro_tools
from galight_modif.data_process import DataProcess
from galight_modif.fitting_specify import FittingSpecify
from galight_modif.fitting_process import FittingProcess
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mpl_interactions.ipyplot as iplt
from matplotlib.widgets import Slider
from download_data import objects

#THIS IS FOR OBJECT0 IN THE i-BAND ONLY
img_path = "/Users/oliviergilbert/Downloads/image-decam-498255-N24-i.fits.gz"
oow_path = "/Users/oliviergilbert/Downloads/iv-decam-498255-N24-i.fits.gz"

img = pyfits.open(img_path)
wht_img = pyfits.open(oow_path)

if False:
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    ax1.imshow(wht_img[1].data)
    ax2.imshow(img[1].data)
    plt.show()

fov_image = img[1].data
header = img[0].header
exp =  astro_tools.read_fits_exp(header)  #Read the exposure time 
wht = wht_img[1].data
mean_wht = exp * (0.0642/0.135)**2  #The drizzle information is used to derive the mean WHT value.
exp_map = exp * wht/mean_wht  #Derive the exposure time map for each pixel


data_process = DataProcess(fov_image = fov_image, target_pos = [1432., 966.], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=True, zp = 24.897)  #zp use 27.0 for convinence.


data_process.generate_target_materials(radius=60, create_mask = True, nsigma=10,
                                      exp_sz= 0.5, npixels = 20, if_plot=True)