from galightFitting import galight_fit
from download_data import objects

objID = int(input("Enter the object ID you want to fit [0-7]:\n"))
type = input("What fitting model do you want to use?\n")
if objID == 0:
    #object0 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "/Users/oliviergilbert/Downloads/image-decam-498255-N24-i.fits.gz",
                oow_path = "/Users/oliviergilbert/Downloads/iv-decam-498255-N24-i.fits.gz",
                type = type,
                median_noise = 863.3368,
                PSF_pos_list = [[ 1617., 1027.], [1502., 667.], [1150., 1437.]],
                )
    

if objID == 1:
    #object1 in i-band:
    galight_fit(ra_dec=objects[objID],
                img_path = "/Users/oliviergilbert/Downloads/image-decam-989756-N29-i.fits.gz",
                oow_path = "/Users/oliviergilbert/Downloads/iv-decam-989756-N29-i.fits.gz",
                type = "AGN",
                median_noise = 3944.6768+0.0128,
                PSF_pos_list = [[ 981., 1927.], [513., 1919.]],
                )