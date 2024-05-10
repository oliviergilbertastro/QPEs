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

objects = [
    (19.7861250,-34.1916944),
    (195.500588,27.7827129),
    (37.9469167,-10.3361972),
    (38.7040417,-44.3254583),
    (189.734917,+33.165911),
    (42.322161985,-4.214494783),
    (210.2222,-28.7665),
    (71.3921,-10.2005),
    ]

def makeFalseColorImage(img):
    """
    Make an RGB false-color image from griz filter bands.
    We choose a visible color to represent each filter.
    """
    g = np.array([0,0,255])
    r = np.array([0,127.5,63.75])
    i = np.array([63.75,127.5,0])
    z = np.array([255,0,0])
    filters = [g,r,i,z]
    rgb_img = np.zeros((img.shape[0],img.shape[1],3))
    for i, filter in enumerate(filters):
        rgb_img[:,:,0] += img[:,:,i]*filter[0]
        rgb_img[:,:,1] += img[:,:,i]*filter[1]
        rgb_img[:,:,2] += img[:,:,i]*filter[2]
    #Now we need to renormalize:
    for i in range(3):
        rgb_img[:,:,i] = rgb_img[:,:,i]/np.max(rgb_img[:,:,i])
    
    print(np.mean(rgb_img[:,:,0]))
    print(np.mean(rgb_img[:,:,1]))
    print(np.mean(rgb_img[:,:,2]))
    return rgb_img



fitsFiles = []
fitsImages = []
fitsFalseColor = []
for i in range(8):
    fitsFiles.append(pyfits.open(f"data/object{i}.fits"))
    fitsImages.append(np.array(fitsFiles[i][0].data).T.swapaxes(0, 1))
    fitsFalseColor.append(makeFalseColorImage(fitsImages[i]))
    print(fitsFalseColor[-1].shape)
    def controlVisual(contrast=0.5, intensity=5):
        return fitsFalseColor[i]**(1/(contrast))*intensity*5
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    ax1 = plt.axes([0.25, 0.1, 0.65, 0.03])
    ax2 = plt.axes([0.25, 0.05, 0.65, 0.03])
    slider1 = Slider(ax1, label="Contrast", valmin=0.01, valmax=5, valinit=3)
    slider2 = Slider(ax2, label="Intensity", valmin=0.1, valmax=20)
    controls = iplt.imshow(controlVisual, contrast=slider1, intensity=slider2, ax=ax, origin="lower")
    plt.title(f"Object_{i}: RA:{(objects[i][0]):.4f} DEC:{(objects[i][1]):.4f}",fontsize=15)
    plt.show()