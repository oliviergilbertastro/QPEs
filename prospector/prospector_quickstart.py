import os
import sys

#Little code to ensure we are in the user directory running this code so the env variable works correctly  
home = os.getcwd()
if home != "/Users/oliviergilbert":
    os.system("cd ../../../ \nsource .bash_profile \npython3 /Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/prospector_quickstart.py")
    sys.exit()

#Need to download fsps from https://github.com/cconroy20/fsps and put it locally, change the path in .bash_profile to point to it

import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import astroquery