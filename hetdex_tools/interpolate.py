import os.path as op
import numpy as np
import tables as tb

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

from hetdex_api.config import HDRconfig
from hetdex_api.extract import Extract

import datetime

try:  # using HDRconfig
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    config = HDRconfig(LATEST_HDR_NAME)
    CONFIG_HDR2 = HDRconfig("hdr2.1")
    CONFIG_HDR3 = HDRconfig("hdr3")
    CONFIG_HDR4 = HDRconfig("hdr4")
    CONFIG_HDR5 = HDRconfig("hdr5")

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr4"

OPEN_DET_FILE = None
DET_HANDLE = None
DET_FILE = None

current_hdr = LATEST_HDR_NAME
surveyh5 = None


def make_narrowband_image(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=None,
    imsize=None,
    wave_range=None,
    convolve_image=False,
    ffsky=False,
    ffsky_rescor=False,
    subcont=False,
    dcont=50.0,
    include_error=False,
    survey=LATEST_HDR_NAME,
    extract_class=None,
    fiber_flux_offset=None,
    interp_kind="linear",
    apply_mask=False,
    fill_value=np.nan,
):
    """
    Function to make narrowband image from either a detectid or from a
    coordinate/shotid combination.

    Paramaters
    ----------
    detectid: int
        detectid from the continuum or lines catalog. Default is
        None. Provide a coords/shotid combo if this isn't given
    coords: SkyCoords object
        coordinates to define the centre of the data cube
    pixscale: astropy angle quantity
         plate scale
    imsize: astropy angle quantity
        image size
    wave_range: list or None
        start and stop value for the wavelength range in Angstrom.
        If not given, the detectid linewidth is used
    convolve_image: bool
        option to convolve image with shotid seeing
    ffsky: bool
        option to use full frame calibrated fibers. Default is
        True.
    ffsky_rescor: bool
        option to use full frame calibrated fibers with residual
        correction. Default False
    subcont: bool
        option to subtract continuum. Default is False. This
        will measure the continuum 50AA below and above the
        input wave_range
    dcont
        width in angstrom to measure the continuum. Default is to
        measure 50 AA wide regions on either side of the line
    include_error bool
        option to include error array
    extract   Extract class object
        option to include a preloaded Extract class object.
        Default is to load extract class according to detection info
    fiber_flux_offset: 1036 array
        array of values in units of 10**-17 ergs/s/cm2/AA to add
        to each fiber spectrum used in the extraction. Defaults
        to None
    interp_kind: str
        Kind of interpolation to pixelated grid from fiber intensity.
        Options are 'linear', 'cubic', 'nearest'. Default is linear.
    apply_mask: bool
        Apply HETDEX fiber mask model. This will mask all fibers contributing
        to the spectral extraction before summation. Masked in place as NaNs
    fill_value: float, optional
        Value used to fill in for requested points outside of coverage or in a mask
        region. If not provided, then the default is nan.

    Returns
    -------
    hdu: PrimaryHDU object
        the 2D summed data array and associated 2d header
        Units are '10^-17 erg cm-2 s-1'
        If include_error=True will include addiional hdu

    Examples
    --------

    For a specific detectid:
    >>> hdu = make_narrowband_image(detectid=2101046271)

    For a SkyCoords object. You must provide shotid and
    wavelength range

    >>> coords = SkyCoord(188.79312, 50.855747, unit='deg')
    >>> wave_obj = 4235.84 #in Angstrom
    >>> hdu = make_narrowband_image(coords=coords,
                                    shotid=20190524021,
                                    wave_range=[wave_obj-10, wave_obj+10])
    """
    global config, current_hdr, surveyh5
    global CONFIG_HDR2, CONFIG_HDR3, CONFIG_HDR4, CONFIG_HDR5, OPEN_DET_FILE, DET_HANDLE
    global DET_FILE

    if survey != current_hdr:
        config = HDRconfig(survey)
        current_hdr = survey

        try:
            if surveyh5 is not None:
                try:
                    surveyh5.close()
                except:
                    pass
            surveyh5 = tb.open_file(config.surveyh5, "r")
        except:
            pass
    else:
        # check to see if survey file is closed and open if so
        try:
            surveyh5.root
        except AttributeError:
            surveyh5 = tb.open_file(config.surveyh5, "r")

    if pixscale is None:
        pixscale = 0.25 * u.arcsec

    if imsize is None:
        imsize = 30.0 * u.arcsec

    if detectid is not None:
        if (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.detecth5
        elif (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.contsourceh5
        elif (detectid >= 3000000000) * (detectid < 3090000000):
            DET_FILE = CONFIG_HDR3.detecth5
        elif (detectid >= 3090000000) * (detectid < 3100000000):
            DET_FILE = CONFIG_HDR3.contsourceh5
        elif (detectid >= 4000000000) * (detectid < 4090000000):
            DET_FILE = CONFIG_HDR4.detecth5
        elif (detectid >= 4090000000) * (detectid < 4100000000):
            DET_FILE = CONFIG_HDR4.contsourceh5
        elif (detectid >= 5000000000) * (detectid < 5090000000):
            DET_FILE = CONFIG_HDR5.detecth5
        elif (detectid >= 5090000000) * (detectid < 5100000000):
            DET_FILE = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:
            DET_HANDLE = tb.open_file(DET_FILE, "r")
            OPEN_DET_FILE = DET_FILE
        else:
            if DET_FILE == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = DET_FILE
                try:
                    DET_HANDLE = tb.open_file(DET_FILE, "r")
                except:
                    print("Could not open {}".format(self.det_file))

        detectid_obj = detectid
        det_info = DET_HANDLE.root.Detections.read_where("detectid == detectid_obj")[0]

        shotid_obj = det_info["shotid"]
        wave_obj = det_info["wave"]
        linewidth = det_info["linewidth"]
        wave_range = [wave_obj - 2.0 * linewidth, wave_obj + 2.0 * linewidth]
        coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

    elif coords is not None:
        if shotid is not None:
            shotid_obj = shotid
        else:
            print("Provide a shotid")
        if wave_range is None:
            print(
                "Provide a wavelength range to collapse. \
            Example wave_range=[4500,4540]"
            )
    else:
        print("Provide a detectid or both a coords and shotid")

    fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
    pa = surveyh5.root.Survey.read_where("shotid == shotid_obj")["pa"][0]

    if extract_class is None:
        E = Extract()
        E.load_shot(shotid_obj, fibers=False, survey=survey)
    else:
        E = extract_class
    # get spatial dims:
    ndim = int(imsize / pixscale)
    center = int(ndim / 2)

    rad = imsize.to(u.arcsec).value  # convert to arcsec value, not quantity

    info_result = E.get_fiberinfo_for_coord(
        coords,
        radius=rad,
        ffsky=ffsky,
        ffsky_rescor=ffsky_rescor,
        fiber_flux_offset=fiber_flux_offset,
        add_mask=apply_mask,
    )

    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

    if include_error:
        zarray = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            error=error,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=wave_range,
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
        )
        imslice = zarray[0]
        imerror = zarray[1]
    else:
        zarray = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=wave_range,
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
        )

        imslice = zarray[0]

    if subcont:
        zarray_blue = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=[wave_range[0] - dcont - 10, wave_range[0] - 10],
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
        )

        zarray_red = E.make_narrowband_image(
            ifux_cen,
            ifuy_cen,
            ifux,
            ifuy,
            data,
            mask,
            seeing_fac=fwhm,
            scale=pixscale.to(u.arcsec).value,
            boxsize=imsize.to(u.arcsec).value,
            wrange=[wave_range[1] + 10, wave_range[1] + dcont + 10],
            convolve_image=convolve_image,
            interp_kind=interp_kind,
            fill_value=fill_value,
        )

        dwave = wave_range[1] - wave_range[0]
        im_cont = (zarray_blue[0] + zarray_red[0]) / (2 * dcont)

        imslice = zarray[0] - dwave * im_cont

    w = wcs.WCS(naxis=2)
    imsize = imsize.to(u.arcsec).value
    w.wcs.crval = [coords.ra.deg, coords.dec.deg]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value]

    # get rotation:
    sys_rot = 1.55
    rot = 360.0 - (90.0 + pa + sys_rot)
    rrot = np.deg2rad(rot)

    #    w.wcs.crota = [ 0, rot]

    w.wcs.pc = [[np.cos(rrot), np.sin(rrot)], [-1.0 * np.sin(rrot), np.cos(rrot)]]

    header = w.to_header()

    # add chosen variable info
    header["FILETIME"] = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    header["SURVEY"] = str(survey)
    header["MASK"] = str(apply_mask)
    header["CONVOLVE"] = str(convolve_image)
    header["INTERP"] = interp_kind
    header["SUBCONT"] = str(subcont)
    header["FFSKY"] = str(ffsky)
    header["FFSKYRC"] = str(ffsky_rescor)
    
    # Copy Shot table info
    shot_info_table = Table(surveyh5.root.Survey.read_where("shotid == shotid_obj"))
    for col in shot_info_table.colnames:
        if col.startswith("x") or col.startswith("y"):
            continue
        if "nstars" in col:
            continue
        if "darktime" in col:
            continue
        if "response" in col:
            header["RESPONSE"] = shot_info_table[col][0]
        elif "virus" in col:
            headername = col.replace("_virus", "")
            try:
                header[headername.upper()] = shot_info_table[col][0]
            except:
                header[headername.upper() + "1"] = shot_info_table[col][0][0]
                header[headername.upper() + "2"] = shot_info_table[col][0][1]
                header[headername.upper() + "3"] = shot_info_table[col][0][2]
        else:
            headername = str
            try:
                header[col.upper()] = shot_info_table[col][0]
            except:
                header[col.upper() + "1"] = shot_info_table[col][0][0]
                header[col.upper() + "2"] = shot_info_table[col][0][1]
                header[col.upper() + "3"] = shot_info_table[col][0][2]

    hdu = fits.PrimaryHDU(imslice, header=header)

    if extract_class is None:
        E.close()

    if include_error:
        hdu_error = fits.ImageHDU(imerror, header=header)
        hdu_x = fits.ImageHDU(zarray[2], header=header)
        hdu_y = fits.ImageHDU(zarray[3], header=header)
        return fits.HDUList([hdu, hdu_error, hdu_x, hdu_y])
    else:
        return hdu


def make_data_cube(
    detectid=None,
    coords=None,
    shotid=None,
    pixscale=None,
    imsize=None,
    wave_range=[3470, 5540],
    dwave=2.0,
    dcont=50.0,
    convolve_image=False,
    ffsky=False,
    ffsky_rescor=False,
    subcont=False,
    survey=LATEST_HDR_NAME,
    fiber_flux_offset=None,
    interp_kind="linear",
    apply_mask=False,
    fill_value=np.nan,
    included_error=False,
):
    """
    Function to make a datacube from either a detectid or from a
    coordinate/shotid combination.

    Paramaters
    ----------
    detectid: int
        detectid from the continuum or lines catalog. Default is
        None. Provide a coords/shotid combo if this isn't given
    coords: SkyCoords object
        coordinates to define the centre of the data cube
    pixscale: astropy angle quantity
        plate scale
    imsize: astropy angle quantity
        spatial length of cube (equal dims is only option)
    wave_range: list
        start and stop value for the wavelength range in Angstrom
    dwave
        step in wavelength range in Angstrom
    convolve_image: bool
         option to convolve image with shotid seeing
    ffsky: bool
        option to use full frame calibrated fibers. Default False
    ffsky_rescor: bool
        option to use full frame calibrated fibers. Default is False
    subcont: bool
        option to subtract continuum. Default is False. This
        will measure the continuum 50AA below and above the
        input wave_range
    dcont: float
        width in angstrom to measure the continuum. Default is to
        measure 50 AA wide regions on either side of the line
    fiber_flux_offset: 1036 array
        array of values in units of 10**-17 ergs/s/cm2/AA to add
        to each fiber spectrum used in the extraction. Defaults
        to None
    interp_kind: str
        Kind of interpolation to pixelated grid from fiber intensity.
        Options are 'linear', 'cubic', 'nearest'. Default is linear.
    apply_mask: bool                                                                                            Apply HETDEX fiber mask model. This will mask all fibers contributing
         to the spectral extraction before summation. Masked in place as NaNs                               fill_value: float, optional                                                                                Value used to fill in for requested points outside of coverage or in a mask                            region. If not provided, then the default is nan.

    Returns
    -------
    hdu: PrimaryHDU object
        the data cube 3D array and associated 3d header
        Units are '10^-17 erg cm-2 s-1 per spaxel'

    Examples
    --------

    Can either pass in a detectid:

    >>> detectid_obj=2101602788
    >>> hdu = make_data_cube( detectid=detectid_obj)
    >>> hdu.writeto( str(detectid_obj) + '.fits', overwrite=True)

    or can put in an SkyCoord object:

    >>> star_coords = SkyCoord(9.625181, -0.043587, unit='deg')
    >>> hdu = make_data_cube( coords=star_coords[0], shotid=20171016108, dwave=2.0)
    >>> hdu.writeto( 'star.fits', overwrite=True)

    """
    global config, current_hdr, surveyh5
    global CONFIG_HDR2, CONFIG_HDR3, OPEN_DET_FILE, DET_HANDLE
    global DET_FILE

    if pixscale is None:
        pixscale = 0.25 * u.arcsec

    if imsize is None:
        imsize = 30.0 * u.arcsec
    
    if survey != current_hdr:
        config = HDRconfig(survey)
        current_hdr = survey

        try:
            if surveyh5 is not None:
                try:
                    surveyh5.close()
                except:
                    pass
            surveyh5 = tb.open_file(config.surveyh5, "r")
        except:
            pass
    else:
        # check to see if survey file is closed and open if so                                                                   
        try:
            surveyh5.root
        except AttributeError:
            surveyh5 = tb.open_file(config.surveyh5, "r")

    if detectid is not None:
        if (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.detecth5
        elif (detectid >= 2100000000) * (detectid < 2190000000):
            DET_FILE = CONFIG_HDR2.contsourceh5
        elif (detectid >= 3000000000) * (detectid < 3090000000):
            DET_FILE = CONFIG_HDR3.detecth5
        elif (detectid >= 3090000000) * (detectid < 3100000000):
            DET_FILE = CONFIG_HDR3.contsourceh5
        elif (detectid >= 4000000000) * (detectid < 4090000000):
            DET_FILE = CONFIG_HDR4.detecth5
        elif (detectid >= 4090000000) * (detectid < 4100000000):
            DET_FILE = CONFIG_HDR4.contsourceh5
        elif (detectid >= 5000000000) * (detectid < 5090000000):
            DET_FILE = CONFIG_HDR5.detecth5
        elif (detectid >= 5090000000) * (detectid < 5100000000):
            DET_FILE = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:
            OPEN_DET_FILE = DET_FILE
            DET_HANDLE = tb.open_file(DET_FILE, "r")

        else:
            if DET_FILE == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = DET_FILE
                try:
                    DET_HANDLE = tb.open_file(DET_FILE, "r")
                except Exception:
                    print("Could not open {}".format(DET_FILE))

        detectid_obj = detectid
        det_info = DET_HANDLE.root.Detections.read_where("detectid == detectid_obj")[0]

        shotid = det_info["shotid"]
        coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

    if coords is None or shotid is None:
        print("Provide a detectid or both a coords and shotid")

    E = Extract()
    E.load_shot(shotid, fibers=False, survey=survey)

    # get spatial dims:
    ndim = int(imsize / pixscale)
    center = int(ndim / 2)

    # get wave dims:
    nwave = int((wave_range[1] - wave_range[0]) / dwave + 1)

    w = wcs.WCS(naxis=3)
    w.wcs.crval = [coords.ra.deg, coords.dec.deg, wave_range[0]]
    w.wcs.crpix = [center, center, 1]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    w.wcs.cdelt = [-pixscale.to(u.deg).value, pixscale.to(u.deg).value, dwave]

    rad = imsize.to(u.arcsec).value

    info_result = E.get_fiberinfo_for_coord(
        coords,
        radius=rad,
        ffsky=ffsky,
        ffsky_rescor=ffsky_rescor,
        fiber_flux_offset=fiber_flux_offset,
        add_mask=apply_mask,
    )

    ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

    # get ifu center:
    ifux_cen, ifuy_cen = E.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coords.ra.deg, coords.dec.deg
    )

    # get FWHM and PA
    shotid_obj = shotid

    pa = surveyh5.root.Survey.read_where("shotid == shotid_obj")["pa"][0]

    if convolve_image:
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
    else:
        fwhm = 1.8  # just a dummy variable as convolve_image=False

    # add in rotation
    sys_rot = 1.55
    rot = 360.0 - (90.0 + pa + sys_rot)
    w.wcs.crota = [0, rot, 0]
    #    rrot = np.deg2rad(rot)
    #    w.wcs.pc = [[np.cos(rrot),
    #                 np.sin(rrot),0],
    #                [-1.0*np.sin(rrot),
    #                 np.cos(rrot),0], [0,0,0]]

    im_cube = np.zeros((nwave, ndim, ndim))

    wave_i = wave_range[0]
    i = 0

    while wave_i <= wave_range[1]:
        try:
            im_src = E.make_narrowband_image(
                ifux_cen,
                ifuy_cen,
                ifux,
                ifuy,
                data,
                mask,
                scale=pixscale.to(u.arcsec).value,
                wrange=[wave_i, wave_i + dwave],
                nchunks=1,
                seeing_fac=fwhm,
                convolve_image=convolve_image,
                boxsize=imsize.to(u.arcsec).value,
                interp_kind=interp_kind,
                fill_value=fill_value,
            )

            im_slice = im_src[0]

            if subcont:
                zarray_blue = E.make_narrowband_image(
                    ifux_cen,
                    ifuy_cen,
                    ifux,
                    ifuy,
                    data,
                    mask,
                    seeing_fac=fwhm,
                    scale=pixscale.to(u.arcsec).value,
                    boxsize=imsize.to(u.arcsec).value,
                    nchunks=2,
                    wrange=[wave_i - dcont, wave_i],
                    convolve_image=convolve_image,
                    interp_kind=interp_kind,
                    fill_value=fill_value,
                )
                zarray_red = E.make_narrowband_image(
                    ifux_cen,
                    ifuy_cen,
                    ifux,
                    ifuy,
                    data,
                    mask,
                    seeing_fac=fwhm,
                    nchunks=2,
                    scale=pixscale.to(u.arcsec).value,
                    boxsize=imsize.to(u.arcsec).value,
                    wrange=[wave_i + dwave, wave_i + dwave + dcont],
                    convolve_image=convolve_image,
                    interp_kind=interp_kind,
                    fill_value=fill_value,
                )

                im_cont = (zarray_blue[0] + zarray_red[0]) / (2 * dcont)
                im_slice = im_src[0] - dwave * im_cont

            im_cube[i, :, :] = im_slice

        except Exception:
            im_cube[i, :, :] = np.zeros((ndim, ndim))
        wave_i += dwave
        i += 1

    # CREATE HEADER
    header = w.to_header()

    # add chosen variable info
    header["FILETIME"] = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    header["SURVEY"] = str(survey)
    header["MASK"] = str(apply_mask)
    header["CONVOLVE"] = str(convolve_image)
    header["INTERP"] = interp_kind
    header["SUBCONT"] = str(subcont)
    header["FFSKY"] = str(ffsky)
    header["FFSKYRC"] = str(ffsky_rescor)
    
    # Copy Shot table info
    shot_info_table = Table(surveyh5.root.Survey.read_where("shotid == shotid_obj"))
    for col in shot_info_table.colnames:
        if col.startswith("x") or col.startswith("y"):
            continue
        if "nstars" in col:
            continue
        if "darktime" in col:
            continue
        if "response" in col:
            header["RESPONSE"] = shot_info_table[col][0]
        elif "virus" in col:
            headername = col.replace("_virus", "")
            try:
                header[headername.upper()] = shot_info_table[col][0]
            except:
                header[headername.upper() + "1"] = shot_info_table[col][0][0]
                header[headername.upper() + "2"] = shot_info_table[col][0][1]
                header[headername.upper() + "3"] = shot_info_table[col][0][2]
        else:
            headername = str
            try:
                header[col.upper()] = shot_info_table[col][0]
            except:
                header[col.upper() + "1"] = shot_info_table[col][0][0]
                header[col.upper() + "2"] = shot_info_table[col][0][1]
                header[col.upper() + "3"] = shot_info_table[col][0][2]

    hdu = fits.PrimaryHDU(im_cube, header=header)

    E.close()

    return hdu
