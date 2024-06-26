"""
DOWNLOADS THE NECESSARY FITS FILES (IMAGES) IN THE GRIZ BANDS
(You need an internet connection to run this code successfully)
"""

import requests
from tqdm import tqdm

def get_file(pos, filenumber, size=512, pixscale=0.262, band="i", name="object"):
    ra, dec = pos
    url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&size={size}&layer=ls-dr10&pixscale={pixscale}&bands={band}&invvar'
    #url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&layer=ls-dr10'
    r = requests.get(url)
    open(r'data/images/'+f'{name}'+str(filenumber)+f"_{band}"+r'.fits' , 'wb').write(r.content)
    url = f'https://www.legacysurvey.org/viewer/coadd-psf/?ra={str(ra)}&dec={str(dec)}&layer=ls-dr10&bands={band}'
    #url = f'https://www.legacysurvey.org/viewer/fits-cutout?ra={str(ra)}&dec={str(dec)}&layer=ls-dr10'
    r = requests.get(url)
    open(r'data/images/'+f'{name}'+str(filenumber)+f"_{band}"+r'_PSF.fits' , 'wb').write(r.content)

#List of objects' RA and DEC to download
#The ones that are commented out are replaced by a more accurate WCS position to their left

#QPE hosts
objects = [
    (19.7861250,-34.1916944),                                                   # GSN 069
    (195.500588,+27.7827129),                                                   # RX J1301.9+2747
    (37.9465, -10.3365), #(37.9469167,-10.3361972),                             # eRO-QPE1
    (38.7029, -44.3258), #(38.7040417,-44.3254583),                             # eRO-QPE2
    (189.734917,+33.165911),                                                    # AT 2019vcb
    (42.322161985,-4.214494783),                                                # 2MASX J0249
    (210.2222, -28.7670), #(210.2222,-28.7665),                                 # eRO-QPE3
    (71.3909, -10.2013), #(71.3921,-10.2005),                                   # eRO-QPE4
    (71.657833, -10.226361)                                                     # AT 2019qiz
    ]

objects_names = [
    "GSN 069",
    "RX J1301.9+2747",
    "eRO-QPE1",
    "eRO-QPE2",
    "AT 2019vcb",
    "2MASX J0249",
    "eRO-QPE3",
    "eRO-QPE4",
    "AT 2019qiz"
]

objects_types = ["AGN", "AGN", "None", "None", "None", "AGN", "None", "AGN", "AGN"]
QPE_bands_list = ["i", "i", "i", "i", "z", "i", "i", "r", "r"]

TDE_names = [
    "ASASSN-14ae",
    "ASASSN-14li",
    "PTF-09ge",
    "RBS 1032",
    "SDSS J1323",
    "SDSS J0748",
    "SDSS J1342",
    "SDSS J1350",
    "SDSS J0952",
    "SDSS J1201",
]

TDE_coords = [
                (167.1671500, 34.0978417),                              #ASASSN-14ae
                (192.0634583, 17.7740139),                              #ASASSN-14li
                (224.2632500, 49.6113806),                              #PTF-09ge
                (176.8612, 49.7160),#(176.8616667, 49.7163889),         #RBS 1032
                (200.9248875, 48.4503500),                              #SDSS J1323
                (117.0861125, 47.2039528),                              #SDSS J0748
                (205.6850667, 5.5155944),                               #SDSS J1342
                (207.5062792, 29.2693639),                              #SDSS J1350
                (148.0398125, 21.7203444),                              #SDSS J0952
                (180.4001167, 30.0515333),                              #SDSS J1201
]

TDE_types = ["None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]

#TDE hosts
comparisons = [
    (192.0634, 17.7739),                                                        #ASASSN-14li
    (224.2632, 49.6113),                                                        #PTF-09ge
    (167.1672, 34.0978),                                                        #ASASSN-14ae
]

comparisons_names = [
    "ASASSN-14li",
    "PTF-09ge",
    "ASASSN-14ae",
]

if __name__ == "__main__":
    if input("Download 1024x1024 of all QPEs? [y/n]") == "y":
        band = input("Which band do you want to download? ['griz']")
        for i in tqdm(range(len(objects_names))):
            for b in band:
                get_file(objects[i], i, size=512*2, pixscale=0.262, band=b, name="qpe")
        print("\x1b[33mDownload finished\x1b[0m")
    elif input("Download 1024x1024 of all TDEs? [y/n]") == "y":
        band = input("Which band do you want to download? ['griz']")
        for i in tqdm(range(len(TDE_names))):
            for b in band:
                get_file(TDE_coords[i], i, size=512*2, pixscale=0.262, band=b, name="tde")
        print("\x1b[33mDownload finished\x1b[0m")
    elif input("Download 3000x3000 of all TDEs? [y/n]") == "y":
        band = input("Which band do you want to download? ['griz']")
        for i in tqdm(range(len(TDE_names))):
            for b in band:
                get_file(TDE_coords[i], i, size=3000, pixscale=0.262, band=b, name="tde")
        print("\x1b[33mDownload finished\x1b[0m")
    elif input("Download custom image of one TDE? [y/n]") == "y":
        objID = int(input("Which object?"))
        size = int(input("Which pixel size do you want? (i.e. 2048, 4096, 8192, etc.)"))
        custom_pos = input("Custom wcs position?\n")
        custom_pos = TDE_coords[objID] if custom_pos == "" else (float(custom_pos.split()[0]),float(custom_pos.split()[1]))
        b = input("Which band?")
        get_file(custom_pos, objID, size=size, pixscale=0.262, band=b, name="tde")
        print("\x1b[33mDownload finished\x1b[0m")
    elif input("Download 1024x1024 of all QPEs? [y/n]") == "y":
        band = input("Which band do you want to download? ['griz']")
        for i in tqdm(range(len(objects))):
            get_file(objects[i], i, size=512*2, pixscale=0.262, band=band)
        print("\x1b[33mDownload finished\x1b[0m")
    elif input("Download 2048x2048 of one QPE? [y/n]") == "y":
        objID = int(input("Which object?"))
        band = input("Which band do you want to download? ['griz']")
        get_file(objects[objID], objID, size=512*4, pixscale=0.262, band=band)
        print("\x1b[33mDownload finished\x1b[0m")