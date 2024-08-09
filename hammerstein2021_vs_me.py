import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import seaborn as sns
from utils import myCombinedFinalPlot
import pandas as pd
from download_data import hammerstein2021_TDE_names, remaining_hammerstein_TDE_names

hammerstein_TDE_data = np.loadtxt("hammerstein2021_TDE_distribution.txt")
my_n_sersic = hammerstein_TDE_data[:,0]


# Hammerstein2021 vs Mine
if input("Hammerstein2021 vs Me? [y/n]") == "y":
    #hammerstein2011_m_star = [10.034548481266134, 10.615111341648774, 9.904298053029024, 10.761646843108398, 9.749318439528095, 9.943050495790008, 9.466651788945452, 10.100352240102294, 10.206235373591642, 9.853392120437164, 10.359435692055296, 9.663851962652899, 10.052793794783952, 9.692079422905257, 9.866209076214142, 10.021731525489155, 10.168297184491784, 9.91464209027962, 10.082167241670486]
    hammerstein2021_n_sersics = [1.6704634721131186, 3.3895522388059702, 3.4369206598586017, 3.4837391987431263, 1.8858601728201099, 2.287824037706206, 1.7633150039277297, 2.4139827179890023, 1.7309505106048704, 7.083346425765907, 4.3969363707776905, 4.421288295365279, 1.9984289080911233, 7.485624509033778, 3.287352710133543, 0.48774548311076205, 3.0940298507462685, 0.5113904163393559, 3.578790259230165]

    # Keep only the hammerstein sersic indices for the ones I have calculated and put them in the correct order:
    hammerstein2021_n_sersic = []
    for i in range(len(remaining_hammerstein_TDE_names)):
        try:
            index = hammerstein2021_TDE_names.index(remaining_hammerstein_TDE_names[i])
            hammerstein2021_n_sersic.append(hammerstein2021_n_sersics[index])
        except:
            pass

    # Correlation plot
    plt.plot(my_n_sersic, hammerstein2021_n_sersic, "o")
    x_plot = np.linspace(np.min(np.concatenate((my_n_sersic, hammerstein2021_n_sersic))),np.max(np.concatenate((my_n_sersic, hammerstein2021_n_sersic))),100)
    plt.plot(x_plot, x_plot, "--", color="red")
    plt.xlabel("My $n$")
    plt.ylabel("Hammerstein 2021 $n$")
    plt.show()




elif input("Simard2011 vs Me? [y/n]") == "y":
    my_n_sersic = np.loadtxt("hammerstein_TDE_allRelevantData_0.txt")[:,1]
    # Simard2011 vs Mine
    simard_names = ['AT2018zr', 'AT2018bsi', 'AT2018hyz', 'AT2019azh', 'AT2020pj', 'AT2020ocn', 'AT2020wey']
    simard_n_sersics = [0.99, 2.85, 2.54, 6.63, 1.11, 0.50, 4.29]
    my_n_sersics = []
    for i in range(len(simard_names)):
        try:
            index = remaining_hammerstein_TDE_names.index(simard_names[i])
            print(index, len(my_n_sersic))
            my_n_sersics.append(my_n_sersic[index])
        except:
            print(simard_names[i])
            pass

    print(my_n_sersics)
    plt.plot(my_n_sersics, simard_n_sersics, "o")
    x_plot = np.linspace(np.min(np.concatenate((my_n_sersics, simard_n_sersics))),np.max(np.concatenate((my_n_sersics, simard_n_sersics))),100)
    plt.plot(x_plot, x_plot, "--", color="red")
    plt.xlabel("My $n$")
    plt.ylabel("Hammerstein 2021 $n$")
    plt.show()