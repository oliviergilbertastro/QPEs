'''
File to organize the data retrieved from the litterature and simplify the code in other files (i.e. compareModels.py)
'''

import numpy as np
import pandas as pd

# Load TDE host galaxies black hole masses and Sérsic indices:
# Only load them if we are in the correct directory, otherwise don't bother, it's not like we're going to use them
import os
home = os.getcwd()
if home != "/Users/oliviergilbert":
    data = np.array(pd.read_csv("data/TDEsersic_mBH.csv"))
    #These are from the Law-Smith paper (https://arxiv.org/pdf/1707.01559) black hole mass vs Sérsic index plot

    # You can also obtain the sérsic indices from the Simard11 decomposition database directly from https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=+J%2FApJS%2F196%2F11%2Ftable3&-from=nav&-nav=cat%3AJ%2FApJS%2F196%2F11%26tab%3A%7BJ%2FApJS%2F196%2F11%2Ftable1%7D%26key%3Asource%3DJ%2FApJS%2F196%2F11%2Ftable1%26HTTPPRM%3A%26
    TDE_mBH = 10**data[:,0]
    TDE_sersicIndices = data[:,1]

TDE_mBH_litterature = [5.566206896551724, 6.173875862068965, 6.262041379310345, 5.6193103448275865, 6.524220689655173, 6.582924137931034, 6.514455172413793, 7.466013793103448, 7.0441379310344825, 8.436993103448277]
TDE_bt_ratio_litterature = [0.5505584522139586, 0.9102572585554349, 0.4201889738771405, 0.7601633009554561, 0.560761453563772, 0.1402416430669878, 0.37060318132493447, 0.2705780377418416, 0.5209485747558426, 0.81998517851944]


TDE_mBH_litterature[9] = 7.18 # French2020 review article
TDE_mBH = 10**np.array(TDE_mBH_litterature)

# French2020 review article
TDE_mBH_litterature = [(5.42,0.46,0.46), # ASASSN-14ae
                       (6.23,0.40,0.39), # ASASSN-14li
                       (6.31,0.39,0.39), # PTF-09ge
                       (5.25,0.62,0.67), # RBS 1032
                       (6.15,0.45,0.46), # SDSS J1323
                       (7.25,0.46,0.45), # SDSS J0748
                       (6.06,0.52,0.51), # SDSS J1342
                       (7.47,0.5,0.5), # SDSS J1350 v_disp = [13.50546932 40.5877037  40.5877037 ]
                       (7.04,0.5,0.5), # SDSS J0952 v_disp = [99.50531769 12.2647686  12.2647686 ]
                       (7.18,0.41,0.41), # SDSS J1201
                       ]
from utils import toLog
TDE_mBH_litterature = toLog(TDE_mBH_litterature, inverse=True)


if False:
    #TDE sersic indices from Review Paper (https://arxiv.org/pdf/2003.02863)
    TDE_sersicIndices = [
                    2.61,         #ASASSN-14ae
                    4.91,        #ASASSN-14li
                    4.03,          #PTF-09ge
                    2.24,        #RBS 1032
                    5.03,         #SDSS J1323
                    1.53,         #SDSS J0748
                    2.56,         #SDSS J1342
                    4.57,         #SDSS J1350
                    7.98,         #SDSS J0952
                    5.61,          #SDSS J1201
    ]

# log(Surface stellar densities) for the 10 TDE hosts used in Law-Smith paper, as found in the Graur paper. My own calculations of the stellar mass surface densities approximately
# give me the same thing as the Graur paper result, which is a good sanity check. https://iopscience.iop.org/article/10.3847/1538-4357/aaa3fd/pdf
TDE_stellarDensities = [
                (9.5, 0.2, 0.2), #ASASSN-14ae
                (10.1, 0.2, 0.2), #ASASSN-14li
                (9.2, 0.3, 0.2), #PTF-09ge
                (9.5, 0.2, 0.2), #RBS 1032
                (9.5, 0.2, 0.2), #SDSS J1323
                (9.5, 0.2, 0.2), #SDSS J0748
                (9.7, 0.2, 0.3), #SDSS J1342
                (9.3, 0.3, 0.3), #SDSS J1350
                (10.4, "upper limit", "upper limit"), #SDSS J0952
                (9.5, 0.3, 0.2), #SDSS J1201
                ]

#redshifts
QPE_redshifts = [
                (0.018),                                # GSN 069
                (0.02358),                              # RX J1301.9+2747
                (0.0505),                               # eRO-QPE1
                (0.0175),                               # eRO-QPE2
                (0.088),                                # AT 2019vcb
                (0.0186),                               # 2MASX J0249
                (0.024),                                # eRO-QPE3
                (0.044),                                # eRO-QPE4
                (0.01513),                              # AT 2019qiz
                ]

old_TDE_redshifts = [
                (0.0206),           #ASASSN-14li
                (0.064),            #PTF-09ge
                (0.0436),           #ASASSN-14ae
                ]

TDE_redshifts = [                   
                0.0436,         #ASASSN-14ae
                0.02058,        #ASASSN-14li
                0.064,          #PTF-09ge
                0.02604,        #RBS 1032
                0.0875,         #SDSS J1323
                0.0615,         #SDSS J0748
                0.0366,         #SDSS J1342
                0.0777,         #SDSS J1350
                0.0789,         #SDSS J0952
                0.146,          #SDSS J1201
                ]

hammerstein_TDE_redshifts = [
                0.075,
                0.051,
                0.088,
                0.212,
                0.046,
                0.138,
                0.091,
                0.340,
                0.193,
                0.121,
                0.022,
                0.051,
                0.074,
                0.148,
                0.152,
                0.117,
                0.015,
                0.087,
                0.068,
                0.088,
                0.160,
                0.070,
                0.159,
                0.070,
                0.093,
                0.345,
                0.435,
                0.027,
                0.057,
                0.277,
]

remaining_hammerstein_TDE_redshifts = [
                0.075,
                0.051,
                0.088,
                0.212,
                0.046,
                0.340,
                0.193,
                0.121,
                0.022,
                0.074,
                0.148,
                0.152,
                0.117,
                0.015,
                0.087,
                0.068,
                0.088,
                0.160,
                0.070,
                0.159,
                0.093,
                0.345,
                0.435,
                0.027,
                0.057,
                0.277,
]
# black hole masses from litterature in log M_sun from MOSFit https://arxiv.org/pdf/2203.01461
hammerstein_TDE_mBH = [
                6.53,
                6.65,
                6.50,
                6.34,
                6.58,
                6.60,
                6.83,
                7.55,
                6.41,
                6.80,
                7.43,
                6.86,
                6.78,
                6.64,
                6.68,
                6.37,
                6.18,
                6.05,
                7.98,
                np.log10(6.5E6),
                7.93,
                7.06,
                6.85,
                6.67,
                6.82,
                7.22,
                7.37,
                7.36,
                6.25,
                7.02,
]

remaining_hammerstein_TDE_mBH = [
                6.53,
                6.65,
                6.50,
                6.34,
                6.58,
                7.55,
                6.41,
                6.80,
                7.43,
                6.78,
                6.64,
                6.68,
                6.37,
                6.18,
                6.05,
                7.98,
                np.log10(6.5E6),
                7.93,
                7.06,
                6.85,
                6.82,
                7.22,
                7.37,
                7.36,
                6.25,
                7.02,
]

# All french data is from https://arxiv.org/pdf/2003.02863
french_TDE_redshifts = [
    0.079,
    0.037,
    0.046,
]

french_TDE_mBH = [
    (6.88,0.38,0.38), # iPTF15af
    (7.00,0.42,0.41), # AT2018dyk
    (5.68,0.52,0.51), # ASASSN18zj
]

TDE_r50s = [          #in kpc       (https://iopscience.iop.org/article/10.3847/1538-4357/aaa3fd/pdf)    
                1.3,            #ASASSN-14ae
                0.4,            #ASASSN-14li
                2.3,            #PTF-09ge
                0.7,            #RBS 1032
                2.7,            #SDSS J1323
                2.3,            #SDSS J0748
                0.9,            #SDSS J1342
                2.0,            #SDSS J1350
                1.0, #low-lim   #SDSS J0952
                3.6,            #SDSS J1201
                ]


#M_BH in solar masses
QPE_mBH = [
                (4E5, 0, 0),                                                        # GSN 069
                (1.8E6, 0.1E6, 0.1E6),                                              # RX J1301.9+2747
                (4E6, 0, 0),                                                        # eRO-QPE1
                (3E6, 0, 0),                                                        # eRO-QPE2
                (6.5E6, 1.5E6, 1.5E6),                                              # AT 2019vcb
                (8.5E4, 0, 0),#or 5E5 depending on paper                            # 2MASX J0249
                (5.3E6, 3.5E6, 0.7E6),                                              # eRO-QPE3
                (6.8E7, 3.2E7, 4.8E7),                                              # eRO-QPE4
                (10**6.18, 10**(6.18)-10**(6.18-0.44), 10**(6.18+0.44)-10**(6.18)), # AT 2019qiz
                
            ]


#M_star in solar masses
QPE_stellar_masses_litterature_sdssProspector = [
                (10**10.74, 10**10.74-10**(10.74-0.11), 10**(10.74+0.05)-10**10.74),# GSN 069               #Prospector (manual magnitudes)
                (10**10.44, 10**10.44-10**(10.44-0.21), 10**(10.44+0.15)-10**10.44),# RX J1301.9+2747       #Prospector
                (3.8E9, 1.9E9, 0.4E9),                                              # eRO-QPE1
                (1.01E9, 0.5E9, 0.01E9),                                            # eRO-QPE2
                (10**10.06, 10**10.06-10**(10.06-0.10), 10**(10.06+0.08)-10**10.06),# AT 2019vcb            #Prospector
                (10**9.56, 10**9.56-10**(9.56-0.04), 10**(9.56+0.05)-10**9.56),     # 2MASX J0249           #Prospector (manual magnitudes)
                (2.56E9, 1.40E9, 0.24E9),                                           # eRO-QPE3
                (1.6E10, 0.6E10, 0.7E10),                                           # eRO-QPE4
                (10**10.04, 10**10.04-10**(10.04-0.14), 10**(10.04+0.10)-10**10.04),# AT 2019qiz
                ]

QPE_stellar_masses_litterature = [
                (None),                                                             # GSN 069
                (None),# RX J1301.9+2747       #Prospector
                (3.8E9, 1.9E9, 0.4E9),                                              # eRO-QPE1
                (1.01E9, 0.5E9, 0.01E9),                                            # eRO-QPE2
                (None),# AT 2019vcb            #Prospector
                (None),     # 2MASX J0249           #Prospector
                (2.56E9, 1.40E9, 0.24E9),                                           # eRO-QPE3
                (1.6E10, 0.6E10, 0.7E10),                                           # eRO-QPE4
                (10**10.04, 10**10.04-10**(10.04-0.14), 10**(10.04+0.10)-10**10.04),# AT 2019qiz
                    ]


old_TDE_stellar_masses_litterature = [
                (10**9.3, 10**9.3-10**9.2, 10**9.4-10**9.3),                                #ASASSN-14li
                (10**9.87, 10**9.87-10**(9.87-0.17), 10**(9.87+0.13)-10**9.87),             #PTF-09ge
                (10**9.73, 10**9.73-10**(9.73-0.13), 10**(9.73+0.13)-10**9.73),           #ASASSN-14ae
                ]


TDE_stellar_masses_litterature = [          #in solar masses        
                (10**9.73, 10**9.73-10**(9.73-0.13), 10**(9.73+0.13)-10**9.73),            #ASASSN-14ae
                (10**9.3, 10**9.3-10**(9.3-0.1), 10**(9.3+0.1)-10**9.3),            #ASASSN-14li
                (10**9.87 , 10**9.87 -10**(9.87 -0.17), 10**(9.87 +0.13)-10**9.87),            #PTF-09ge
                (10**9.19, 10**9.19-10**(9.19-0.16), 10**(9.19+0.15)-10**9.19),            #RBS 1032
                (10**10.38, 10**10.38-10**(10.38-0.07), 10**(10.38+0.06)-10**10.38),            #SDSS J1323
                (10**10.18, 10**10.18-10**(10.18-0.09), 10**(10.18+0.06)-10**10.18),            #SDSS J0748
                (10**9.64, 10**9.64-10**(9.64-0.07), 10**(9.64+0.23)-10**9.64),            #SDSS J1342
                (10**9.94, 10**9.94-10**(9.94-0.2), 10**(9.94+0.17)-10**9.94),            #SDSS J1350
                (10**10.37, 10**10.37-10**(10.37-0.07), 10**(10.37+0.06)-10**10.37),   #SDSS J0952
                (10**10.61, 10**10.61-10**(10.61-0.16), 10**(10.61+0.08)-10**10.61),            #SDSS J1201
                ]


#Schlafy et al. extinction (https://irsa.ipac.caltech.edu/workspace/TMP_3ah3nH_9869/DUST/19.7861250_-34.1916944.v0001/extinction.html)
QPE_extinction = {
    # QPE hosts
    "GSN 069":          {"g":0.089,"r":0.062,"i":0.046,"z":0.034},
    "RX J1301.9+2747":  {"g":0.030,"r":0.021,"i":0.015,"z":0.011},
    "eRO-QPE1":         {"g":0.080,"r":0.055,"i":0.041,"z":0.030},
    "eRO-QPE2":         {"g":0.060,"r":0.042,"i":0.031,"z":0.023},
    "AT 2019vcb":       {"g":0.058,"r":0.040,"i":0.030,"z":0.022},
    "2MASX J0249":      {"g":0.101,"r":0.070,"i":0.052,"z":0.039},
    "eRO-QPE3":         {"g":0.192,"r":0.133,"i":0.099,"z":0.073},
    "eRO-QPE4":         {"g":0.280,"r":0.193,"i":0.144,"z":0.107},
    "AT 2019qiz":       {"g":0.372,"r":0.257,"i":0.191,"z":0.142},
    }

#Schlafy et al. extinction (https://irsa.ipac.caltech.edu/applications/DUST/)
TDE_extinction = {
    # TDE hosts
    "ASASSN-14ae":          {"g":0.058,"r":0.040,"i":0.030,"z":0.022},
    "ASASSN-14li":          {"g":0.085,"r":0.058,"i":0.043,"z":0.032},
    "PTF-09ge":             {"g":0.055,"r":0.038,"i":0.028,"z":0.021},
    "RBS 1032":             {"g":0.049,"r":0.034,"i":0.025,"z":0.019},
    "SDSS J1323":           {"g":0.027,"r":0.019,"i":0.014,"z":0.010},
    "SDSS J0748":           {"g":0.222,"r":0.154,"i":0.114,"z":0.085},
    "SDSS J1342":           {"g":0.082,"r":0.057,"i":0.042,"z":0.031},
    "SDSS J1350":           {"g":0.050,"r":0.035,"i":0.026,"z":0.019},
    "SDSS J0952":           {"g":0.093,"r":0.064,"i":0.048,"z":0.035},
    "SDSS J1201":           {"g":0.061,"r":0.042,"i":0.031,"z":0.023},
    }



#Schlafy et al. extinction (https://irsa.ipac.caltech.edu/applications/DUST/)
hammerstein_TDE_extinction = {
    # TDE hosts
    "AT2018zr":           {"g":0.198,"r":0.107,"i":0.079,"z":0.059},
    "AT2018bsi":          {"g":0.203,"r":0.140,"i":0.104,"z":0.078},
    "AT2018hco":          {"g":0.132,"r":0.091,"i":0.068,"z":0.050},
    "AT2018iih":          {"g":0.162,"r":0.112,"i":0.083,"z":0.062},
    "AT2018hyz":          {"g":0.114,"r":0.079,"i":0.058,"z":0.043},
    "AT2018jbv":          {"g":0.107,"r":0.074,"i":0.055,"z":0.041},
    "AT2019cho":          {"g":0.045,"r":0.031,"i":0.023,"z":0.017},
    "AT2019bhf":          {"g":0.082,"r":0.057,"i":0.042,"z":0.031},
    "AT2019azh":          {"g":0.146,"r":0.101,"i":0.075,"z":0.056},
    "AT2019ehz":          {"g":0.059,"r":0.041,"i":0.031,"z":0.023},
    "AT2019mha":          {"g":0.026,"r":0.018,"i":0.013,"z":0.010},
    "AT2019meg":          {"g":0.185,"r":0.128,"i":0.095,"z":0.071},
    "AT2019lwu":          {"g":0.122,"r":0.084,"i":0.063,"z":0.047},
    "AT2019qiz":          {"g":0.372,"r":0.257,"i":0.191,"z":0.142},
    "AT2019teq":          {"g":0.182,"r":0.126,"i":0.094,"z":0.070},
    "AT2020pj":           {"g":0.107,"r":0.074,"i":0.055,"z":0.041},
    "AT2019vcb":          {"g":0.058,"r":0.040,"i":0.030,"z":0.022},
    "AT2020ddv":          {"g":0.043,"r":0.029,"i":0.022,"z":0.016},
    "AT2020ocn":          {"g":0.062,"r":0.043,"i":0.032,"z":0.024},
    "AT2020opy":          {"g":0.214,"r":0.148,"i":0.110,"z":0.082},
    "AT2020mbq":          {"g":0.157,"r":0.109,"i":0.081,"z":0.060},
    "AT2020qhs":          {"g":0.071,"r":0.049,"i":0.037,"z":0.027},
    "AT2020riz":          {"g":0.285,"r":0.197,"i":0.147,"z":0.109},
    "AT2020wey":          {"g":0.157,"r":0.109,"i":0.081,"z":0.060},
    "AT2020zso":          {"g":0.224,"r":0.155,"i":0.115,"z":0.086},
    "AT2020ysg":          {"g":0.057,"r":0.039,"i":0.029,"z":0.022},
    }

#Schlafy et al. extinction (https://irsa.ipac.caltech.edu/applications/DUST/)
french_TDE_extinction = {
    # TDE hosts
    "PTF-15af":           {"g":0.111,"r":0.077,"i":0.057,"z":0.043},
    "AT2018dyk":          {"g":0.068,"r":0.047	,"i":0.035,"z":0.026},
    "ASASSN18zj":          {"g":0.114,"r":0.079,"i":0.058,"z":0.043},
    }