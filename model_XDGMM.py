#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 00:57:45 2018

@author: johnruan
"""

from numpy import *
from scipy import *
from pylab import *
from xdgmm import XDGMM
import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

alpha_r = []
alpha_r_err = []

inputfile = open('alpha_r.txt', 'r')
line = inputfile.read()
data1 = line.split()
for i in range(0, len(data1), 2):
    alpha_r.append(float(data1[i]))
    alpha_r_err.append(float(data1[i+1]))
alpha_r = asarray(alpha_r)
alpha_r_err = asarray(alpha_r_err)

alpha_r = alpha_r.reshape(alpha_r.shape[0], 1)    
alpha_r_err = alpha_r_err.reshape(alpha_r_err.shape[0], 1,1)


xdgmm1 = XDGMM()
param_range = np.array([1,2,3,4,5, 6, 7, 8, 9, 10])
bic, optimal_n_comp, lowest_bic = xdgmm1.bic_test(alpha_r, alpha_r_err, param_range)
print(bic, optimal_n_comp)

xdgmm1.n_components = optimal_n_comp
xdgmm1 = xdgmm1.fit(alpha_r, alpha_r_err)
#print(xdgmm1)
resampled_alpha_r = xdgmm1.sample(10000000)


ax = plt.subplot(111)

hist(alpha_r, histtype='step', log=True, bins=15, range=(-5,5), color='k', label='dC sample',  ls=':', lw=2, density=True)# weights=[1./float(len(RV_dC))]*len(RV_dC))  
hist(resampled_alpha_r, histtype='step', log=True, bins=1000, range=(-5,5), color='k', label='dC sample deconvolved', lw=2, density=True) #weights=[1./float(len(resampled_RV_dC))]*len(resampled_RV_dC))

print(mean(resampled_alpha_r), std(resampled_alpha_r), 0.126*0.126)

ax.set_xlabel(r'RV [km/s]', size=14)
ax.set_ylabel(r'PDF', size=14)
ax.tick_params(labelsize=12)
legend = plt.legend(loc=2, ncol=1, fontsize=9)

show()


