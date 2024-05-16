from galightFitting import galight_fit
from download_data import objects

# If the necessary files are not in the "data/science/" folder, download them from the folowing link (replacing RA and DEC by the ra and dec of the object):
# https://www.legacysurvey.org/viewer/data-for-radec/?ra=RA&dec=DEC&layer=ls-dr10
# You should go to the section "CCDs overlapping RA,Dec:" and download the filter band you want. You need the image and the weight map (oow)


objID = int(input("Enter the object ID you want to fit [0-7]:\n"))
band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
type = input("What fitting model do you want to use?\n")

data_repo = "data/science/"

#-------------------------------------------------------------------object0---------------------------------------------------------
if objID == 0:
    if band == "i" or band == 'I' or band == "2":
        img_path = f"image-decam-498255-N24-i.fits.gz"
        oow_path = f"iv-decam-498255-N24-i.fits.gz"
    elif band == "z" or band == 'Z' or band == "3":
        img_path = f"image-decam-498254-N24-z.fits.gz"
        oow_path = f"iv-decam-498254-N24-z.fits.gz"
    galight_fit(ra_dec=objects[objID],
                img_path = data_repo+img_path,
                oow_path = data_repo+oow_path,
                type = type,
                median_noise = 0,#863.3368,
                PSF_pos_list=[[19.7915, -34.2054], [19.7957, -34.1628], [19.7597, -34.1971]], #We find stars in the image online, click on them and copy their WCS coordinates here
                band=band,
                )

#-------------------------------------------------------------------object1---------------------------------------------------------
elif objID == 1:
    #object1 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "data/science/image-decam-989756-N29-i.fits.gz",
                oow_path = "data/science/iv-decam-989756-N29-i.fits.gz",
                type = type,
                median_noise = 3944.6768+0.0128,
                PSF_pos_list = [[ 981., 1927.], [513., 1919.]],
                band=band,
                )

#-------------------------------------------------------------------object2---------------------------------------------------------
elif objID == 2:
    #object1 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "data/science/image-decam-989756-N29-i.fits.gz",
                oow_path = "data/science/iv-decam-989756-N29-i.fits.gz",
                type = type,
                median_noise = 3944.6768+0.0128,
                PSF_pos_list = [[ 981., 1927.], [513., 1919.]],
                band=band,
                )