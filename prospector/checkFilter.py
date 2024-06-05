import os
import sys

#Little code to ensure we are in the user directory running this code so the env variable works correctly  
#home = os.getcwd()
#if home != "/Users/oliviergilbert":
#    os.system("cd ../../../ \nsource .bash_profile \npython3 /Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/checkFilter.py")
#    sys.exit()
#os.system("cd ../../../ \nsource .bash_profile")
#os.system("source ../../../.bash_profile")
#os.system("export SPS_HOME=/Users/oliviergilbert/Desktop/QPEs/fsps-master")
os.environ["SPS_HOME"] = "/Users/oliviergilbert/Desktop/QPEs/fsps-master"
import astropy.table
import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import astroquery
from sedpy.observate import load_filters
from prospect.utils.obsutils import fix_obs



bands = "ugriz"

filters = load_filters([f"sdss_{b}0" for b in bands])
print(filters[0])