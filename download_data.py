"""
DOWNLOADS THE NECESSARY FITS FILES (IMAGES) IN THE GRIZ BANDS
(You need an internet connection to run this code successfully)
"""

import requests
from tqdm import tqdm

def get_file(pos, filenumber, size=512, pixscale=0.262):
    ra, dec = pos
    url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&size={size}&layer=ls-dr10&pixscale={pixscale}'
    r = requests.get(url)
    open(r'data/object'+str(filenumber)+r'.fits' , 'wb').write(r.content)

#List of objects' RA and DEC to download
objects = [
    (19.7861250,-34.1916944),
    (195.500588,27.7827129),
    (37.9469167,-10.3361972),
    (38.7040417,-44.3254583),
    (189.734917,+33.165911),
    (42.322161985,-4.214494783),
    (210.2222,-28.7665),
    (71.3921,-10.2005),
    ]

if __name__ == "__main__":
    for i in tqdm(range(len(objects))):
        get_file(objects[i], i, size=1024, pixscale=0.131)