import numpy as np
import matplotlib.pyplot as plt

# Distributions should be [bt_ratio, n_sersic, surface_mass_density, mStar, mBH]
# The program distributions.py creates the txt files of the QPE and TDE distributions

QPE_distribution = np.loadtxt("QPE_distribution.txt")
TDE_distribution = np.loadtxt("TDE_distribution.txt")
reference_distribution = np.loadtxt("referenceCatalog_modif.txt")
reference_distribution = reference_distribution[:,[12,60,68,63,67]]

def getStats(dist):
    n_data = dist.shape[0] # number of data points
    means = np.mean(dist, axis=0)
    stdevs = np.std(dist, axis=0)
    print(n_data, means, stdevs)
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
    different_params = ["B/T ratio", "n_sersic", "SMSD", "M_star", "M_BH"]
    for i in range(len(different_params)):
        print(f"\x1b[33m{different_params[i]}:\x1b[0m")
        z = z_statistic(N_qpe, mu_qpe[i], std_qpe[i], N_tde, mu_tde[i], std_tde[i])
        string = '\x1b[32mSame\x1b[0m' if z < 2.33 else '\x1b[31mDifferent\x1b[0m'
        print(f"QPE = TDE : {string} ({z})")
        z = z_statistic(N_qpe, mu_qpe[i], std_qpe[i], N_ref, mu_ref[i], std_ref[i])
        string = '\x1b[32mSame\x1b[0m' if z < 2.33 else '\x1b[31mDifferent\x1b[0m'
        print(f"QPE = ref : {string} ({z})")
        z = z_statistic(N_tde, mu_tde[i], std_tde[i], N_ref, mu_ref[i], std_ref[i])
        string = '\x1b[32mSame\x1b[0m' if z < 2.33 else '\x1b[31mDifferent\x1b[0m'
        print(f"TDE = ref : {string} ({z})")



def c(alpha):
    return np.sqrt(-np.log(alpha/2)/2)

def F(sample, x):
    n_sample = len(sample)
    subsample = np.array(sample)[np.array(sample) <= x]
    n_subsample = len(subsample)
    return n_subsample/n_sample

def KolmogorovSmirnov(sample1, sample2, if_plot=True):
    # Combine and sort the samples
    combined = list(np.concatenate((sample1, sample2)))
    combined.sort()
    absolute_differences = []
    for x in combined:
        absolute_differences.append(np.abs(F(sample1, x) - F(sample2, x)))
    D_nm = max(absolute_differences)
    if if_plot:
        index = absolute_differences.index(D_nm)
        plt.plot([combined[index],combined[index]], [F(sample1, combined[index]),F(sample2, combined[index])], "--", linewidth=3, color="black")
        plt.step(combined, [F(sample1, x) for x in combined], label="$F_1(x)$")
        plt.step(combined, [F(sample2, x) for x in combined], label="$F_2(x)$")
        plt.legend(fontsize=15)
        plt.xlabel("$x$", fontsize=17)
        plt.ylabel("Cumulative probability", fontsize=17)
        plt.show()
    return D_nm

def rejectNullHypothesis(D, n, m, alpha=0.2, verbose=True):
    # The null hypothesis is rejected if
    rejectBool = D > c(alpha)*np.sqrt((n+m)/(n*m))
    if verbose:
        string = '\x1b[32mSame\x1b[0m' if not rejectBool else '\x1b[31mDifferent\x1b[0m'
        print(string)
    return rejectBool

# Common alpha values:
# alpha is basically the probability the test is wrong
alphas = [0.20,0.15,0.10,0.05,0.025,0.01,0.005,0.001]



if __name__ == "__main__":
    print("*************************************")
    for i in range(len(different_params)):
        print(f"\x1b[33m{different_params[i]}:\x1b[0m")
        print(f"QPE = TDE :", end=" ")
        dnm = KolmogorovSmirnov(QPE_distribution[:,i], TDE_distribution[:,i])
        rejectNullHypothesis(dnm, len(QPE_distribution), len(TDE_distribution))
        print(f"QPE = ref :", end=" ")
        dnm = KolmogorovSmirnov(QPE_distribution[:,i], reference_distribution[:,i])
        rejectNullHypothesis(dnm, len(QPE_distribution), len(reference_distribution))
        print(f"TDE = ref :", end=" ")
        dnm = KolmogorovSmirnov(TDE_distribution[:,i], reference_distribution[:,i])
        rejectNullHypothesis(dnm, len(TDE_distribution), len(reference_distribution))


    # test from website
    #print(z_statistic(51.5,8,25,39.5,7,25))