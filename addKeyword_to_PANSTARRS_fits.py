import astropy.io.fits as pyfits
from download_data import TDE_names
from tqdm import tqdm



panstarrIDs = ["2044.052","1724.046","2326.041","2318.040","2322.010","2246.076","1463.036","1976.039","1799.049","1970.059"]




key_add = input("Add keyword to header?") == "y"
for band in "griz":
    for i in tqdm(range(len(TDE_names))):
        data_repo = f"data/science/tde{i}/"
        img = data_repo+f"rings.v3.skycell.{panstarrIDs[i]}.stk.{band}.unconv.fits"
        exp_map = data_repo+f"rings.v3.skycell.{panstarrIDs[i]}.stk.{band}.unconv.exp.fits"
        img_hdul = pyfits.open(img, mode="update")
        exp_hdul = pyfits.open(exp_map, mode="update")
        if key_add:
            img_hdul[1].header.append(('RADESYS', 'FK5', 'wcs keyword'))
            exp_hdul[1].header.append(('RADESYS', 'FK5', 'wcs keyword'))
            img_hdul.writeto(img, overwrite=True)
            exp_hdul.writeto(exp_map, overwrite=True)
        print(img_hdul[1].header)