import numpy as np

# Distributions should be [mBH, mStar, bt_ratio, n_sersic, surface_mass_density]
# The program distributions.py creates the txt files of the QPE and TDE distributions

QPE_distribution = np.loadtxt("QPE_distribution.txt")
TDE_distribution = np.loadtxt("TDE_distribution.txt")
reference_distribution = np.loadtxt("referenceCatalog_modif.txt")
reference_distribution = reference_distribution[:,[67,67,12,60,68]]

def getStats(dist):
    n_data = dist.shape[0] # number of data points
    means = np.mean(dist, axis=0)
    stdevs = np.std(dist, axis=0)
    return n_data, means, stdevs

def z_statistic(mu1, sig1, ndata1, mu2, sig2, ndata2, verbose=False):
    """
    mu1, mu2: means of 1st and 2nd distributions
    sig1, sig2: standard deviations of 1st and 2nd distributions
    ndata1, ndata2: number of data points in 1st and 2nd distributions

    If the Z-statistic is less than 2, the two samples are the same.
    If the Z-statistic is between 2.0 and 2.5, the two samples are marginally different
    If the Z-statistic is between 2.5 and 3.0, the two samples are significantly different
    If the Z-statistic is more then 3.0, the two samples are highly signficantly different

    In many fields of social science, the cricital value of Z (Zcrit) is 2.33 which represents p =0.01 (shaded area below has
    1% of the total area of the probability distribution function). This means the two distributions have only a 1% probability
    of being actual the same.

    http://homework.uoregon.edu/pub/class/es202/ztest.html
    """
    if verbose:
        print("X1:", mu1)
        print("X2:", mu2)
        print("X1-X2:", mu1-mu2)
        print("sigma_1:", (sig1/np.sqrt(ndata1)))
        print("sigma_2:", (sig2/np.sqrt(ndata2)))
        print("sqrt of sig1^2+sig2^2:", np.sqrt((sig1/np.sqrt(ndata1))**2+(sig2/np.sqrt(ndata2))**2))
    z = (mu1-mu2)/np.sqrt((sig1/np.sqrt(ndata1))**2+(sig2/np.sqrt(ndata2))**2)
    #print(z)
    return np.abs(z)


N_qpe, mu_qpe, std_qpe = getStats(QPE_distribution)
N_tde, mu_tde, std_tde = getStats(TDE_distribution)
N_ref, mu_ref, std_ref = getStats(reference_distribution)


if __name__ == "__main__":
    different_params = ["B/T ratio", "n_sersic", "SMSD", "M_BH", "M_star"]
    for i in range(len(different_params)):
        print(f"\x1b[33m{different_params[i]}:\x1b[0m")
        print("QPE = TDE :", z_statistic(N_qpe, mu_qpe[i], std_qpe[i], N_tde, mu_tde[i], std_tde[i]) < 2.33)
        print("QPE = ref :", z_statistic(N_qpe, mu_qpe[i], std_qpe[i], N_ref, mu_ref[i], std_ref[i]) < 2.33)
        print("TDE = ref :", z_statistic(N_tde, mu_tde[i], std_tde[i], N_ref, mu_ref[i], std_ref[i]) < 2.33)
    #print(z_statistic(51.5,8,25,39.5,7,25))