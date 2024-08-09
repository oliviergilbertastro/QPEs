import numpy as np
from download_data import *

ham_data = np.loadtxt("hammerstein_TDE_allRelevantData_0.txt")

for i in range(len(remaining_hammerstein_TDE_names)):
    print(f"{i} {remaining_hammerstein_TDE_names[i]} {np.around(ham_data[i,0], decimals=2)} {np.around(ham_data[i,1], decimals=2)}")