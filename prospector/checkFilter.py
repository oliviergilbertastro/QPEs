import numpy as np
from sedpy.observate import load_filters
import matplotlib.pyplot as plt


bands = "ugriz"

filters = load_filters([f"sdss_{b}0" for b in bands])
for filter in filters:
    ax1 = plt.subplot(111)
    ax1.plot(filter.wavelength, filter.transmission, label=filter.name)
ax1.set_xlabel(r"Wavelength [$\AA$]", fontsize=16)
ax1.set_ylabel(r"Transmission", fontsize=16)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
plt.legend(fontsize=14)
plt.show()

bands = ["H","J","Ks"]
filters = load_filters([f"twomass_{b}" for b in bands])
for filter in filters:
    ax1 = plt.subplot(111)
    ax1.plot(filter.wavelength, filter.transmission, label=filter.name)
ax1.set_xlabel(r"Wavelength [$\AA$]", fontsize=16)
ax1.set_ylabel(r"Transmission", fontsize=16)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
plt.legend(fontsize=14)
plt.show()