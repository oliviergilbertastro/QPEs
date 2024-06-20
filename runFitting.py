from galightFitting import galight_fit
from download_data import objects, comparisons, objects_names, objects_types, TDE_names, TDE_coords

# If the necessary files are not in the "data/science/" folder, download them from the folowing link (replacing RA and DEC by the ra and dec of the object):
# https://www.legacysurvey.org/viewer?ra=RA&dec=DEC&layer=ls-dr10&zoom=15
# You should go to the section "CCDs overlapping RA,Dec:" and download the filter band you want. You need the image (ooi) and the weight map (oow)

#The QPE host galaxies are named "obj" while the TDE host galaxies are named "comp"

if input("Fit a PANSTARRS TDE host galaxy? [y/n]\n") == "y":
    objID = int(input(f"Enter the object ID you want to fit [0-{len(objects)-1}]:\n"))
    band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")
    data_repo = f"data/science/tde{objID}/"
    survey = "PS"
    savename = f"{TDE_names[objID]}_{band}-band_{type}_{survey}"
    pixel_scale = 0.250
    panstarrID = ["2044.052",None,None,None,None,None,None,None,None,None][objID]
    stars = [[[167.1107, 34.1297], [167.1150, 34.1387]],
                    [],
                    None,
                    [[176.8422, 49.6955],[176.8693, 49.7345]],
                    [[200.9640, 48.4047], [200.9628, 48.3902], [201.0256, 48.4088]],
                    [[117.1183, 47.1922],[117.0543, 47.1920]],
                    None,
                    None,
                    None,
                    None,
                    ][objID]


    galight_fit(ra_dec=objects[objID],
            img_path = data_repo+f"rings.v3.skycell.{panstarrID}.stk.{band}.unconv.fits",
            exp_path = data_repo+f"rings.v3.skycell.{panstarrID}.stk.{band}.unconv.exp.fits",
            type = type,
            PSF_pos_list=stars, #We find stars in the image online, click on them and copy their WCS coordinates here
            band=band,
            survey=survey,
            savename=savename,
            pixel_scale=pixel_scale,
            )


elif input("Fit a  CO-ADDED image TDE host galaxy? [y/n]\n") == "y":
    for i, name in enumerate(TDE_names):
        print(f"{i}: {name}")
    objID = int(input(f"Enter the object ID you want to fit [0-{len(objects)-1}]:\n"))
    band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")
    img_path = f"data/images/tde{objID}_{band}.fits"
    #Find stars manually, or leave None if you want galight to search for them itself
    list_of_PSFs = [None,#167.20061582  34.16038963 with 20px radius
                    None,#192.16658903  17.76412554
                    None,#224.33033068  49.59058316
                    None,#176.99564392  49.75402138
                    None,#201.02569809  48.40885903
                    None,#117.01561576  47.25207367
                    None,#205.73562533   5.50056366
                    None,#207.42361071  29.20576714
                    None,
                    None,
                    ]
    galight_fit(ra_dec=TDE_coords[objID],
                    img_path = img_path,
                    oow_path = None,
                    type = type,
                    PSF_pos_list=list_of_PSFs[objID], #We find stars in the image online, click on them and copy their WCS coordinates here
                    band=band,
                    survey="COADDED_DESI",
                    savename=f"{TDE_names[objID]}_{band}-band_{type}_DESI",
                    threshold=5,
                    nsigma=15,
                    exp_sz_multiplier=1,
                    fitting_level="deep",
                    )



elif input("Fit a  CO-ADDED image QPE host galaxy? [y/n]\n") == "y":
    objID = int(input(f"Enter the object ID you want to fit [0-{len(objects)-1}]:\n"))
    band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")
    img_path = f"data/images/object{objID}_{band}.fits"
    #Find stars manually, or leave None if you want galight to search for them itself
    list_of_PSFs = [None,
                    None,
                    None,
                    [[38.7643, -44.3211], [38.6746, -44.3405], [38.6773, -44.3321]],
                    None,   #[[189.6887, 33.2096], [189.7192, 33.1615], [189.7368, 33.1792]],
                    [[42.3057, -4.1936], [42.3714, -4.1891], [42.3621, -4.1921]],
                    None,
                    None,
                    None,
                    ]
    #obj2 ra and dec should be (37.9466, -10.3362)
    #obj5 ra and dec should be (42.3222,-4.2146)
    #obj7 ra and dec should be (71.3910, -10.2014)
    #obj8 ra and dec should be (71.6579, -10.2264)
    galight_fit(ra_dec=objects[objID],
                    img_path = img_path,
                    oow_path = None,
                    type = type,
                    PSF_pos_list=list_of_PSFs[objID], #We find stars in the image online, click on them and copy their WCS coordinates here
                    band=band,
                    survey="COADDED_DESI",
                    savename=f"{objects_names[objID]}_{band}-band_{type}_DESI",
                    threshold=5,
                    nsigma=10,
                    exp_sz_multiplier=1,
                    fitting_level="deep",
                    )




#CHANGE ALL THAT TO FIT PANSTARRS IMAGES

elif input("Fit a PANSTARRS QPE host galaxy? [y/n]\n") == "y":
    objID = int(input(f"Enter the object ID you want to fit [0-{len(objects)-1}]:\n"))
    band = "r"
    type = objects_types[objID]
    data_repo = f"data/science/obj{objID}/"
    survey = "PANSTARRS"
    savename = f"{objects_names[objID]}_{band}-band_{type}_{survey}"
    pixel_scale = 0.250
    panstarrID = [None,"1893.099","1063.041",None,"2049.024","1153.090","0681.083","1072.048","1072.047"][objID]
    stars = [None,
             None,
             None,
             None,
             None,
             [[42.3057, -4.1935]],
             None,
             None,
             None][objID]

    if objID == 0:
        raise ValueError(f"No PANSTARRS image for {objects_names[objID]}")

    galight_fit(ra_dec=objects[objID],
            img_path = data_repo+f"cutout_rings.v3.skycell.{panstarrID}.stk.r.unconv.fits",
            exp_path = data_repo+f"cutout_rings.v3.skycell.{panstarrID}.stk.r.unconv.exp.fits",
            type = type,
            PSF_pos_list=stars, #We find stars in the image online, click on them and copy their WCS coordinates here
            band=band,
            survey=survey,
            savename=savename,
            pixel_scale=pixel_scale,
            )







































elif input("Fit a QPE host galaxy? [y/n]\n") == "y":
    objID = int(input(f"Enter the object ID you want to fit [0-{len(objects)-1}]:\n"))
    band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")

    data_repo = f"data/science/obj{objID}/"

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
                PSF_pos_list=[[19.7915, -34.2054], [19.7957, -34.1628], [19.7597, -34.1971]], #We find stars in the image online, click on them and copy their WCS coordinates here
                band=band,
                nsigma=3,
                )

    #-------------------------------------------------------------------object1---------------------------------------------------------
    elif objID == 1:
        if band == "i" or band == 'I' or band == "2":
            img_path = f"image-decam-989756-N29-i.fits.gz"
            oow_path = f"iv-decam-989756-N29-i.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = [[195.5072, 27.7707], [195.4652, 27.7651], [195.5622, 27.7844], [195.5709, 27.7705], [195.4288, 27.7695]],
                    band=band,
                    )

    #-------------------------------------------------------------------object2---------------------------------------------------------
    elif objID == 2:
        if  band == "r" or band == 'R' or band == "1":
            img_path = f"image-decam-499464-N20-r.fits.gz"
            oow_path = f"iv-decam-499464-N20-r.fits.gz"
        elif band == "i" or band == 'I' or band == "2":
            img_path = f"image-decam-499463-N20-i.fits.gz"
            oow_path = f"iv-decam-499463-N20-i.fits.gz"
        elif band == "z" or band == 'Z' or band == "3":
            img_path = f"image-decam-495390-N20-z.fits.gz"
            oow_path = f"iv-decam-495390-N20-z.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = [[37.9556, -10.3576], [38.0116, -10.3474]],
                    band=band,
                    )
        

    #-------------------------------------------------------------------object3---------------------------------------------------------
    elif objID == 3:
        if band == "i" or band == 'I' or band == "2":
            img_path = f"image-decam-667213-N27-i.fits.gz"
            oow_path = f"iv-decam-667213-N27-i.fits.gz"
        elif band == "z" or band == 'Z' or band == "3":
            img_path = f"image-decam-513359-N22-z.fits.gz"
            oow_path = f"iv-decam-513359-N22-z.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = [[38.7643, -44.3211], [38.8033, -44.3424], [38.8122, -44.3468], [38.8552, -44.3693]],
                    band=band,
                    nsigma=3,
                    )

    #-------------------------------------------------------------------object4---------------------------------------------------------
    elif objID == 4:
        if band == "z" or band == 'Z' or band == "3":
            img_path = f"image-decam-830049-N24-z.fits.gz"
            oow_path = f"iv-decam-830049-N24-z.fits.gz"
            survey = None
            exppath = None
            stars = [[189.6886, 33.2096], [189.6298, 33.2235], [189.7368, 33.1792]]
        if band == "r" or band == 'R' or band == "2":
            img_path = f"image-90prime-74030109-CCD4-r.fits.gz"
            oow_path = f"iv-90prime-74030109-CCD4-r.fits.gz"
            survey = None
            exppath = None
            stars = [[189.6886, 33.2096], [189.6298, 33.2235], [189.7368, 33.1792]]
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = stars,
                    band=band,
                    exp_path=exppath,
                    survey="DESI"
                    )

    #-------------------------------------------------------------------object5---------------------------------------------------------
    elif objID == 5:
        if band == "i" or band == 'I' or band == "2":
            img_path = f"image-decam-791708-N5-i.fits.gz"
            oow_path = f"iv-decam-791708-N5-i.fits.gz"
        elif band == "z" or band == 'Z' or band == "3":
            img_path = f"image-decam-611386-N13-z.fits.gz"
            oow_path = f"iv-decam-611386-N13-z.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = None, #[[42.3057, -4.1935], [42.3336, -4.1881], [42.2590, -4.1946]],
                    band=band,
                    nsigma=3,
                    )

    #-------------------------------------------------------------------object6---------------------------------------------------------
    elif objID == 6:
        if band == "i" or band == 'I' or band == "2":
            img_path = f"image-decam-977812-N10-i.fits.gz"
            oow_path = f"iv-decam-977812-N10-i.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = [[210.2263, -28.7413], [210.2434, -28.7409]],
                    band=band,
                    nsigma=10,
                    radius=35
                    )

        
    #-------------------------------------------------------------------object7---------------------------------------------------------
    elif objID == 7:
        if band == "r" or band == 'R' or band == "1":
            img_path = f"image-decam-384902-N2-r.fits.gz"
            oow_path = f"iv-decam-384902-N2-r.fits.gz"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list = [[71.4243, -10.2004], [71.4131, -10.2141], [71.4304, -10.2135]],
                    band=band,
                    radius=45,
                    )
        
    #-------------------------------------------------------------------object8---------------------------------------------------------
    elif objID == 8:
        survey = "DESI"
        exp_path = None
        if band == "r" or band == "R" or band == "1":
            img_path = f"image-decam-781092-S28-r.fits.gz"
            oow_path = f"iv-decam-781092-S28-r.fits.gz"
        elif band == "z" or band == 'Z' or band == "3":
            img_path = f"image-decam-780799-S28-z.fits.gz"
            oow_path = f"iv-decam-780799-S28-z.fits.gz"
        elif band == "g" or band == "G" or band == "0":
            img_path = f"cutout_rings.v3.skycell.1072.047.stk.g.unconv.fits"
            oow_path = f"cutout_rings.v3.skycell.1072.047.stk.g.unconv.wt.fits"
            survey = "PANSTARRS"
            exp_path = data_repo+f"cutout_rings.v3.skycell.1072.047.stk.g.unconv.exp.fits"
        galight_fit(ra_dec=objects[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    exp_path = exp_path,
                    type = type,
                    PSF_pos_list = [],
                    band=band,
                    radius=60,
                    nsigma=15,
                    exp_sz_multiplier=2,
                    npixels=10,
                    survey=survey,
                    pixel_scale=0.250,
                    )








elif input("Fit a TDE host r-band raw? [y/n]\n") == "y":
    for i in range(len(TDE_names)):
        print(f"{i}: {TDE_names[i]}")
    objID = int(input(f"Enter the object ID you want to fit [0-{len(TDE_names)-1}]:\n"))
    band = "r"
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")
    data_repo = f"data/science/tde{objID}/"
    survey = "DESI"
    savename = f"{TDE_names[objID]}_{band}-band_{type}_{survey}_raw"

    desi_id = [
        f"90prime-74210166-CCD3",
        f"decam-722768-S14",
        f"90prime-75500054-CCD3",
        f"90prime-78650135-CCD4",
        f"90prime-75210068-CCD3",
        f"90prime-78630035-CCD3",
        f"decam-635311-S30",
        f"decam-720906-N8",
        f"decam-433134-N12",
        f"decam-535150-N14",
    ]
    pixel_scale = [0.453,0.262,0.453,0.453,0.453,0.453,0.262,0.262,0.262,0.262][objID]

    img_path = f"image-{desi_id[objID]}-{band}.fits.gz"
    oow_path = f"iv-{desi_id[objID]}-{band}.fits.gz"

    galight_fit(ra_dec=TDE_coords[objID],
                img_path = data_repo+img_path,
                oow_path = data_repo+oow_path,
                type = type,
                PSF_pos_list=None,#[[192.0682, 17.7950]], #We find stars in the image online, click on them and copy their WCS coordinates here
                band=band,
                survey=survey,
                savename=savename,
                pixel_scale=pixel_scale,
                threshold=2,
                )




elif input("Fit a TDE host raw (outdated)? [y/n]\n") == "y":
    objID = int(input(f"Enter the object ID you want to fit [0-{len(comparisons)-1}]:\n"))
    band = input("Enter the filter band you want to fit [g,r,i,z]:\n")
    type = input("What extra-component fitting model do you want to use [None, AGN, Bulge, Bulge+AGN]?\n")

    data_repo = f"data/science/comp{objID}/"

    #-------------------------------------------------------------------object0---------------------------------------------------------
    if objID == 0:
        if band == "r" or band == 'R' or band == "1":
            observation = input("Which observation?\n    1. 2018-02-18 08:48:35\n    2. 2017-03-27 05:12:26\n    3. 2015-04-12 04:14:11\n    4. 2018-03-16 06:48:37\n")
            if observation == "1":
                #2018-02-18 08:48:35	 Current images:
                img_path = f"image-decam-722768-S14-r.fits.gz"
                oow_path = f"iv-decam-722768-S14-r.fits.gz"
                star = [[191.83713361,17.78680452]]

            elif observation == "2":
                #2017-03-27 05:12:26
                img_path = f"image-decam-634095-S15-r.fits.gz"
                oow_path = f"iv-decam-634095-S15-r.fits.gz"
                star = [[191.83713361,17.78680452]]

            elif observation == "3":
                #2015-04-12 04:14:11
                img_path = f"image-decam-432059-N31-r.fits.gz"
                oow_path = f"iv-decam-432059-N31-r.fits.gz"
                star = [[192.21190281,17.77426765]]
            
            elif observation == "4":
                #2018-03-16 06:48:37
                img_path = f"image-decam-730868-S14-r.fits.gz"
                oow_path = f"iv-decam-730868-S14-r.fits.gz"
                star = None

        galight_fit(ra_dec=comparisons[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    PSF_pos_list=star,#[[192.0682, 17.7950]], #We find stars in the image online, click on them and copy their WCS coordinates here
                    band=band,
                    )
        
    elif objID == 1:
        if band == "r" or band == "R" or band == "1":
            img_path = f"image-90prime-75500054-CCD3-r.fits.gz"
            oow_path = f"iv-90prime-75500054-CCD3-r.fits.gz"
        galight_fit(ra_dec=comparisons[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    pixel_scale=0.453,
                    PSF_pos_list=[[224.2723, 49.6162]], #We find stars in the image online, click on them and copy their WCS coordinates here
                    band=band,
                    )
        
    elif objID == 2:
        if band == "g" or band == "G" or band == "0":
            img_path = f"image-90prime-78920101-CCD3-g.fits.gz"
            oow_path = f"iv-90prime-78920101-CCD3-g.fits.gz"
        galight_fit(ra_dec=comparisons[objID],
                    img_path = data_repo+img_path,
                    oow_path = data_repo+oow_path,
                    type = type,
                    pixel_scale=0.453,
                    PSF_pos_list=[[167.1150, 34.1386], [167.1107, 34.1297], [167.1312, 34.1600]], #We find stars in the image online, click on them and copy their WCS coordinates here
                    band=band,
                    )