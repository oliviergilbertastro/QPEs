from galightFitting import galight_fit
from download_data import objects

# If the necessary files are not in the "data/science/" folder, download them from the folowing link (replacing RA and DEC by the ra and dec of the object):
# https://www.legacysurvey.org/viewer/data-for-radec/?ra=RA&dec=DEC&layer=ls-dr10
# You should go to the section "CCDs overlapping RA,Dec:" and download the filter band you want. You need the image and the weight map (oow)


objID = int(input("Enter the object ID you want to fit [0-7]:\n"))
type = input("What fitting model do you want to use?\n")
if objID == 0:
    #object0 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "data/science/image-decam-498255-N24-i.fits.gz",
                oow_path = "data/science/iv-decam-498255-N24-i.fits.gz",
                type = type,
                median_noise = 863.3368,
                PSF_pos_list = [[ 1617., 1027.], [1502., 667.], [1150., 1437.]],
                )

elif objID == 1:
    #object1 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "data/science/image-decam-989756-N29-i.fits.gz",
                oow_path = "data/science/iv-decam-989756-N29-i.fits.gz",
                type = type,
                median_noise = 3944.6768+0.0128,
                PSF_pos_list = [[ 981., 1927.], [513., 1919.]],
                )

elif objID == 2:
    #object1 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "data/science/image-decam-989756-N29-i.fits.gz",
                oow_path = "data/science/iv-decam-989756-N29-i.fits.gz",
                type = type,
                median_noise = 3944.6768+0.0128,
                PSF_pos_list = [[ 981., 1927.], [513., 1919.]],
                )