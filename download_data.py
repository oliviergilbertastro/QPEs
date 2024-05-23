"""
DOWNLOADS THE NECESSARY FITS FILES (IMAGES) IN THE GRIZ BANDS
(You need an internet connection to run this code successfully)
"""

import requests
from tqdm import tqdm

def get_file(pos, filenumber, size=512, pixscale=0.262):
    ra, dec = pos
    url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&size={size}&layer=ls-dr10&pixscale={pixscale}'
    #url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&layer=ls-dr10'
    r = requests.get(url)
    open(r'data/images/object'+str(filenumber)+r'.fits' , 'wb').write(r.content)

#List of objects' RA and DEC to download
#The ones that are commented out are replaced by a more accurate WCS position to their left
objects = [
    (19.7861250,-34.1916944),                                                   # GSN 069
    (195.500588,+27.7827129),                                                   # RX J1301.9+2747
    (37.9465, -10.3365), #(37.9469167,-10.3361972),                             # eRO-QPE1
    (38.7029, -44.3258), #(38.7040417,-44.3254583),                             # eRO-QPE2
    (189.734917,+33.165911),                                                    # AT 2019vcb
    (42.322161985,-4.214494783),                                                # 2MASX J0249
    (210.2222, -28.7670), #(210.2222,-28.7665),                                 # eRO-QPE3
    (71.3909, -10.2013), #(71.3921,-10.2005),                                   # eRO-QPE4
    ]


if __name__ == "__main__":
    for i in tqdm(range(len(objects))):
        get_file(objects[i], i, size=1024, pixscale=0.262)