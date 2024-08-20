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
from utils import recombine_arrays
from paper_data import *

if not sys.warnoptions:
    warnings.simplefilter("ignore")

# load QPE and TDE data
qpe_load_0 = np.loadtxt("QPE_allRelevantData_0.txt")
qpe_load_1 = np.loadtxt("QPE_allRelevantData_1.txt")
qpe_load_2 = np.loadtxt("QPE_allRelevantData_2.txt")
QPE_fullData = recombine_arrays(qpe_load_0,qpe_load_1,qpe_load_2)

tde_load_0 = np.loadtxt("TDE_allRelevantData_0.txt")
tde_load_1 = np.loadtxt("TDE_allRelevantData_1.txt")
tde_load_2 = np.loadtxt("TDE_allRelevantData_2.txt")
TDE_fullData = recombine_arrays(tde_load_0,tde_load_1,tde_load_2)
TDE_fullData = TDE_fullData[:10,:,:]

french_tde_load_0 = np.loadtxt("french_TDE_allRelevantData_0.txt")
french_tde_load_1 = np.loadtxt("french_TDE_allRelevantData_1.txt")
french_tde_load_2 = np.loadtxt("french_TDE_allRelevantData_2.txt")
french_TDE_fullData = recombine_arrays(french_tde_load_0,french_tde_load_1,french_tde_load_2)
french_TDE_fullData = french_TDE_fullData[:,:,:]


# load reference catalog
refCat = np.loadtxt(f"referenceCatalog_with_uncertainties_final_9bins.txt")

QPE_stellar_masses, TDE_stellar_masses, french_TDE_stellar_masses = QPE_fullData[:,3,:], TDE_fullData[:,3,:], french_TDE_fullData[:,3,:]
QPE_mBH, TDE_mBH, french_TDE_mBH = QPE_fullData[:,4,:], TDE_fullData[:,4,:], french_TDE_fullData[:,4,:]
TDE_stellar_masses = np.concatenate((TDE_stellar_masses, french_TDE_stellar_masses))
TDE_redshifts = np.concatenate((TDE_redshifts, french_TDE_redshifts))
TDE_mBH = np.concatenate((TDE_mBH, french_TDE_mBH))
ref_mBH = refCat[:,100:102]


qpe = asarray(QPE_mBH[:,0])
qpe_err = asarray((QPE_mBH[:,1]+QPE_mBH[:,2])/2)
qpe = qpe.reshape(qpe.shape[0], 1)    
qpe_err = qpe_err.reshape(qpe_err.shape[0], 1,1)
xdgmm1 = XDGMM()
param_range = np.array([1,2,3,4,5, 6, 7, 8, 9])
bic, optimal_n_comp, lowest_bic = xdgmm1.bic_test(qpe, qpe_err, param_range)
print(bic, optimal_n_comp)
xdgmm1.n_components = optimal_n_comp
xdgmm1 = xdgmm1.fit(qpe, qpe_err)
resampled_qpe = xdgmm1.sample(10000)

tde = asarray(TDE_mBH[:,0])
tde_err = asarray((TDE_mBH[:,1]+TDE_mBH[:,2])/2)
tde = tde.reshape(tde.shape[0], 1)    
tde_err = tde_err.reshape(tde_err.shape[0], 1,1)
xdgmm1 = XDGMM()
param_range = np.array([1,2,3,4,5, 6, 7, 8, 9])
bic, optimal_n_comp, lowest_bic = xdgmm1.bic_test(tde, tde_err, param_range)
print(bic, optimal_n_comp)
xdgmm1.n_components = optimal_n_comp
xdgmm1 = xdgmm1.fit(tde, tde_err)
resampled_tde = xdgmm1.sample(10000)

ref = asarray(ref_mBH[:,0])
ref_err = asarray(ref_mBH[:,1])
ref = ref.reshape(ref.shape[0], 1)    
ref_err = ref_err.reshape(ref_err.shape[0], 1,1)
xdgmm1 = XDGMM()
param_range = np.array([1,2,3,4,5, 6, 7, 8, 9]) # [113544.58122997 112185.61955142 112228.41833261 112202.3906818
                                                #  112195.1773325  112154.5768182  111832.85299957 111867.68623161
                                                #  111898.1784371 ] 7

param_range = [7]
bic, optimal_n_comp, lowest_bic = xdgmm1.bic_test(ref, ref_err, param_range)
print(bic, optimal_n_comp)
xdgmm1.n_components = optimal_n_comp
xdgmm1 = xdgmm1.fit(ref, ref_err)
resampled_ref = xdgmm1.sample(10000)


ax = plt.subplot(111)

hist(qpe, histtype='step', log=True, bins=15, range=(4,10), color='blue', label='QPE sample',  ls=':', lw=2, density=True)# weights=[1./float(len(RV_dC))]*len(RV_dC))  
hist(resampled_qpe, histtype='step', log=True, bins=1000, range=(4,10), color='blue', label='QPE sample deconvolved', lw=2, density=True) #weights=[1./float(len(resampled_RV_dC))]*len(resampled_RV_dC))
hist(tde, histtype='step', log=True, bins=15, range=(4,10), color='red', label='TDE sample',  ls=':', lw=2, density=True)# weights=[1./float(len(RV_dC))]*len(RV_dC))  
hist(resampled_tde, histtype='step', log=True, bins=1000, range=(4,10), color='red', label='TDE sample deconvolved', lw=2, density=True) #weights=[1./float(len(resampled_RV_dC))]*len(resampled_RV_dC))
hist(ref, histtype='step', log=True, bins=15, range=(4,10), color='k', label='ref sample',  ls=':', lw=2, density=True)# weights=[1./float(len(RV_dC))]*len(RV_dC))  
hist(resampled_ref, histtype='step', log=True, bins=1000, range=(4,10), color='k', label='ref sample deconvolved', lw=2, density=True) #weights=[1./float(len(resampled_RV_dC))]*len(resampled_RV_dC))

print(mean(resampled_qpe), std(resampled_qpe))

ax.set_xlabel(r'$\log(M_\mathrm{BH}) [M_\odot]$', size=14)
ax.set_ylabel(r'PDF', size=14)
ax.tick_params(labelsize=12)
legend = plt.legend(loc=2, ncol=1, fontsize=9)

show()


