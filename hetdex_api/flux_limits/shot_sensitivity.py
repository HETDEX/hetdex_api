"""

Unlike the sensitivity cubes, this module
generates flux limit estimates on the fly
from the Fiber class.

AUTHOR: Daniel Farrow

"""

import logging 
from os.path import isfile, join
from tempfile import NamedTemporaryFile

from numpy import (arange, meshgrid, ones, array, sum, sqrt, square, newaxis,
                   convolve, repeat, loadtxt, where, argmin, invert, logical_not,
                   isnan, newaxis, ceil, nansum, savetxt, transpose)
from numpy.ma import MaskedArray

import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from pyhetdex.coordinates.tangent_projection import TangentPlane
from pyhetdex.het.fplane import FPlane

from hetdex_api.survey import Survey
from hetdex_api.config import HDRconfig
from hetdex_api.extract import Extract
from hetdex_api.flux_limits.generate_simulation_inputs import create_sensitivity_cube_from_astrom
from hetdex_api.mask import create_gal_ellipse, amp_flag_from_fiberid, meteor_flag_from_coords, create_dummy_wcs
from hetdex_api.flux_limits.sensitivity_cube import WavelengthException 
from hetdex_api.flux_limits.flim_models import return_flux_limit_model

def bad_central_mask(weights, badmask, index, wmin = 0.1):
    """
    Karl Gebhardt's false positive cut. If any fiber
    with a weight above wmin is flagged within the
    central three wave bins of a source remove it.

    Parameters
    ----------
    weights : array
        the weights of the fibers
    badmask : 2d numpy array (nfibers x nwaves)
        True where data is bad
    index : int
        the wavelength index of the 
    wmin : float (optional)
        the minimum weight of fiber to
        consider when looking for flagged
        values

    Returns
    -------
    maskval : boolean
        False if bad, True if data is
        good

    """

 
    # Weights of fibers at the source index
    sel = weights[:, index] > wmin


    #print(weights[:, index])
    #print(badmask[sel, index - 1 : index + 2])

    if any(badmask[sel, index - 1 : index + 2].ravel()):
        return False
    else:
        return True



class ShotSensitivity(object):
    """
    Generate completeness estimates for a shot
    on the fly, using the Extract class. This
    is written to be as backward compatible as
    possible with the user interface of 
    hdf5_sensitivity_cubes:SensitivityCubeHDF5Container 
    
    A lot of this is adapted from 
   `hetdex_tools/get_spec.py` `hetdex_api/extract.py`
    and scripts written by Erin Mentuch Cooper.

    The aperture correction stuff was taken
    from 'Remedy' by Greg Zeimann:
    grzeimann Remedy 
   
    Parameters
    ----------
    datevshot : str
        the name of the shot in the form
        YYYYMMDDvSSS
    release : str (Optional)
        The name of the release e.g. hdr2.1
        defaults to latest
    flim_model : str
        The flux limit model to convert
        the noise to completeness, defaults
        to the latest (which might not be
        compatible with the release, see README)
    rad : float
        A radius in arcseconds to grab fibers
        from when computing the flux limit, 
        default 3.5
    ffsky : boolean
        Full frame sky subtraction (default: False)
    wavenpix : int
        Number of wave pixels either side of the
        pixel the source is on to add in
        quadrature, or sets the size of tophat 
        convolution when producing data cubes as 
        2*wavenpix + 1 (default 3)
    d25scale : float
        Sets the multiplier for the galaxy masks
        applied (default 3.0)
    sclean_bad : bool
        Replace bad data using the sclean
        tool (see hetdex_api.extract:Extract)
    verbose : bool
        Print information about the flux limit model
        to the screen
    """
    def __init__(self, datevshot, release=None, flim_model=None, rad=3.5, 
                 ffsky=False, wavenpix=3, d25scale=3.0, verbose=False,
                 sclean_bad = True, log_level="WARNING"): 

        self.conf = HDRconfig()
        self.extractor = Extract()
        self.shotid = int(datevshot.replace("v", ""))
        self.date = datevshot[:8]
        self.rad = rad
        self.ffsky = ffsky
        self.wavenpix = wavenpix
        self.sclean_bad = sclean_bad

        logger = logging.getLogger(name="ShotSensitivity")
        logger.setLevel(log_level)

        if verbose:
            raise DeprecationWarning("Using verbose is deprecated, set log_level instead")    
            logger.setLevel("DEBUG")        
   
        logger.info("shotid: {:d}".format(self.shotid))
 
        if not release:
            self.release = self.conf.LATEST_HDR_NAME
        else:
            self.release = release
           
        logger.info("Data release: {:s}".format(self.release)) 
        self.survey = Survey(survey=self.release)
   
        # Set up flux limit model
        self.f50_from_noise, self.sinterp, interp_sigmas \
                                       = return_flux_limit_model(flim_model, 
                                                                 cache_sim_interp=False,
                                                                 verbose=verbose)
 
        # Generate astrometry for this shot
        survey_sel = (self.survey.shotid == self.shotid)
        self.shot_pa = self.survey.pa[survey_sel][0]
        self.shot_ra = self.survey.ra[survey_sel][0]
        self.shot_dec = self.survey.dec[survey_sel][0]
        rot = 360.0 - (self.shot_pa + 90.)
        self.tp = TangentPlane(self.shot_ra, self.shot_dec, rot)
    
        #Set up masking
        logger.info("Using d25scale {:f}".format(d25scale))
        self.setup_mask(d25scale)
 
 
        # Set up spectral extraction
        if release == "hdr1":
            fwhm = self.survey.fwhm_moffat[survey_sel][0]
        else:
            fwhm = self.survey.fwhm_virus[survey_sel][0]

        logger.info("Using Moffat PSF with FWHM {:f}".format(fwhm))
        self.moffat = self.extractor.moffat_psf(fwhm, 3.*rad, 0.25)
        self.extractor.load_shot(self.shotid, fibers=True, survey=self.release)
        
        # Set up the focal plane astrometry
        fplane_table = self.extractor.shoth5.root.Astrometry.fplane

        # Bit of a hack to avoid changing pyhetdex
        with NamedTemporaryFile(mode='w') as tpf:
            for row in fplane_table.iterrows():
                tpf.write("{:03d} {:8.5f} {:8.5f} {:03d} {:03d} {:03d} {:8.5f} {:8.5f}\n".format(
                            row['ifuslot'], row['fpx'], row['fpy'], row['specid'],
                            row['specslot'], row['ifuid'], row['ifurot'], row['platesc']))
            tpf.seek(0)
            self.fplane = FPlane(tpf.name)

        
    def setup_mask(self, d25scale):
        """
        Setup the masking, to speed up checking
        if sources are in the mask later. This
        is run at initialisation, so you need
        only run again to change `d25scale`

        Parameters
        ----------
        d25scale : float
            Sets the multiplier for the galaxy masks
            applied (default 3.5)
        """      

        logger = logging.getLogger(name="ShotSensitivity")
 
        # see if this is a bad shot
        #print("Bad shot from ", self.conf.badshot)
        badshot = loadtxt(self.conf.badshot, dtype=int)
        badtpshots = loadtxt(self.conf.lowtpshots, dtype=int)
        if (self.shotid in badshot) or (self.shotid in badtpshots):
            logger.warn("Shot is in bad. Making mask zero everywhere")
            self.badshot = True
        else:
            self.badshot = False
        
        # set up bad amps
        logger.info("Bad amps from {:s}".format(self.conf.badamp))
        self.bad_amps = Table.read(self.conf.badamp)
        sel_shot = (self.bad_amps["shotid"] == self.shotid)
        self.bad_amps = self.bad_amps[sel_shot]
        
        # set up galaxy mask
        logger.info("Galaxy mask from {:s}".format(self.conf.rc3cat))
        galaxy_cat = Table.read(self.conf.rc3cat, format='ascii')
        gal_coords = SkyCoord(galaxy_cat['Coords'], frame='icrs')
        shot_coords = SkyCoord(ra=self.shot_ra, dec=self.shot_dec,
                               unit="deg")
        sel_reg = where(shot_coords.separation(gal_coords) < 1.*u.deg)[0]

        self.gal_regions = []
        if len(sel_reg) > 0:
            for idx in sel_reg:
                self.gal_regions.append(create_gal_ellipse(galaxy_cat, 
                                                           row_index=idx, 
                                                           d25scale=d25scale))
                
        # set up meteor mask
        # check if there are any meteors in the shot:
        logger.info("Meteors from {:s}".format(self.conf.meteor))
        self.met_tab = Table.read(self.conf.meteor, format="ascii")
        self.met_tab = self.met_tab[self.shotid == self.met_tab["shotid"]]
            
    def extract_ifu_sensitivity_cube(self, ifuslot, nx=31, ny=31, 
                                    ifusize=62, generate_sigma_array=True,
                                    cache_sim_interp=True):
        """       
        Extract the sensitivity cube
        from IFU `ifuslot`

        Parameters
        ----------
        ifuslot : string
            the IFU slot to extract
        nx, ny : int
            the dimensions in pixels 
            of the cube (default 31,31)
        ifusize : float
            the length of the side of
            the cube in arcseconds,
            default is 62 arcseconds
        generate_sigma_array: bool
            this fills the 3D array
            of noise, this makes it quite
            slow to run, so if you want
            to just use the cube for 
            the astrometry do not use
            this option (default: True)
        cache_sim_interp : bool
            cache the simulation interpolator
            to speed up execution (default: True)

        Returns
        -------
        scube : hetdex_api.flux_limits.sensitivity_cube:SensitivityCube
            the sensitivity cube
 
        """

        waves = self.extractor.get_wave()
        wrange = [waves[0], waves[-1]]
        nz = len(waves)
        
        pa = self.shot_pa
        ifu = self.fplane.difus_ifuslot[ifuslot.replace("ifuslot_", "")]
        ra_ifu, dec_ifu = self.tp.xy2raDec(ifu.y, ifu.x)
        
        scube = create_sensitivity_cube_from_astrom(float(ra_ifu), float(dec_ifu), 
                                                    pa, nx, ny, nz, ifusize, 
                                                    wrange=wrange, 
                                                    cache_sim_interp=cache_sim_interp)
        
        if generate_sigma_array:
            
            ix, iy = meshgrid(arange(0, nx, 1.0),
                              arange(0, ny, 1.0))
        

            all_ra, all_dec, junk = scube.wcs.all_pix2world(ix.ravel(), iy.ravel(), 
                                                            [500.], 0)
            noises, norm, mask = self.get_f50(all_ra, all_dec, None, 1.0, 
                                              direct_sigmas = True)
            sigmas = noises.ravel(order="F").reshape(nz, ny, nx)

            mask = logical_not(mask.reshape(ny, nx))
            mask3d = repeat(mask[newaxis, :, :], sigmas.shape[0], axis=0)
            scube.sigmas = MaskedArray(sigmas, mask=mask3d, fill_value=999.0)
        
        return scube
       

    def get_f50(self, ra, dec, wave, sncut, direct_sigmas = False, 
                nmax = 5000, return_amp = False, linewidth=None):
        """
        Return flux at 50% for the input positions
        most of this is cut and paste from 
        `hetdex_tools/get_spec.py` and 
        the old sensitivity_cube:SensitivityCube 
        class.

        This class splits the data up into sets of
        up to `nmax` to save memory

        Parameters
        ----------
        ra, dec : array 
            right ascension and dec in degrees
        wave : array 
            wavelength in Angstroms. If None,
            then return flux limits for
            all wave bins.
        sncut : float
            cut in detection significance 
            that defines this catalogue
        direct_sigmas : bool
            return the noise values directly
            without passing them through
            the noise to 50% completeness 
            flux (default = False)
        nmax : int
            maximum number of sources to 
            consider at once, otherwise split
            up and loop.
        return_amp : bool
            if True return amplifier information
            for the closest fiber to each source
            (default = False)
        linewidth : array
            optionally pass the linewidth of
            the source (in AA) to activate the linewidth
            dependent part of the completeness
            model (default = None).

        Returns
        -------
        f50s : array
            50% completeness. If outside
            of cube return 999. If None
            was passed for wave this is 
            a 2D array of ra, dec and
            all wavelengths
        mask : array
            Only returned if `wave=None`. This
            mask is True where the ra/dec
            positions passed are in good
            regions of data
        amp : array
            Only returned in `return_amp=True`,
            it's an array of amplifier information
            for the closest fiber to each source
        """
        if type(wave) != type(None):
            wave_passed = True
        else:
            wave_passed = False
 

        try:
            nsrc = len(ra)

            if wave_passed:   
                # Trim stuff very far away
                gal_coords = SkyCoord(ra=ra, dec=dec, unit="deg")
                shot_coords = SkyCoord(ra=self.shot_ra, dec=self.shot_dec,
                                       unit="deg")
        
                sel = array(shot_coords.separation(gal_coords) < 2.0*u.deg)
    
                ra_sel = array(ra)[sel]
                dec_sel = array(dec)[sel]
                wave_sel = array(wave)[sel]
                nsel = len(ra_sel)
            else:
                # If not passing wave always loop
                # over all ra/dec in range
                ra_sel = ra
                dec_sel = dec
                wave_sel = None
                nsel = len(ra)
 
            nsplit = int(ceil(float(nsel)/float(nmax)))

        except TypeError as e:

            # If the user does not pass arrays
            nsplit = 1
            nsrc = 1
            nsel = 1
            sel = True
            ra_sel = array([ra])
            dec_sel = array([dec])
            wave_sel = array([wave])

        # Array to store output actually in the shot
        f50s_sel = []
        mask_sel = []
        amp_sel = []
        norm_sel = []

        wave_rect = self.extractor.get_wave()
        pixsize_aa = wave_rect[1] - wave_rect[0]

        # This will give 999 once the noise is scaled suitably
        badval = 999*1e17/pixsize_aa

        # Arrays to store full output
        f50s = badval*ones(nsrc)
        mask = ones(nsrc)
        norm = ones(nsrc)
        amp = array(["notinshot"]*nsrc)

        if nsel > 0:
            for i in range(nsplit):  
     
                tra = ra_sel[i*nmax : (i+1)*nmax]
                tdec = dec_sel[i*nmax : (i+1)*nmax]
    
                if wave_passed:
                    twave = wave_sel[i*nmax : (i+1)*nmax]
                    if not self.badshot:
                        tf50s, tamp, tnorm = self._get_f50_worker(tra, tdec, twave, sncut, 
                                                                  direct_sigmas = direct_sigmas,
                                                                  linewidth = linewidth)
                    else:
                        tamp = ["bad"]*len(tra)
                        tf50s = [badval]*len(tra)
                        tnorm = [1.0]*len(tra)
                else:
                    # if bad shot then the mask is all set to zero
                    tf50s, tmask, tamp, tnorm = \
                                      self._get_f50_worker(tra, tdec, None, sncut, 
                                                          direct_sigmas = direct_sigmas,
                                                          linewidth = linewidth) 
    
                    mask_sel.extend(tmask)
          
                f50s_sel.extend(tf50s)
                amp_sel.extend(tamp)
                norm_sel.extend(tnorm)
                 

        if return_amp:
            if wave_passed:

                # copy to output
                f50s[sel] = f50s_sel
                amp[sel] = amp_sel
                norm[sel] = norm_sel

                return f50s, norm, amp       
            else:
                return array(f50s_sel), array(norm_sel), array(mask_sel), array(amp_sel)
        else:
            if wave_passed:
                f50s[sel] = f50s_sel
                norm[sel] = norm_sel

                return f50s, norm
            else:
                return array(f50s_sel), array(norm_sel), array(mask_sel)
 

    def _get_f50_worker(self, ra, dec, wave, sncut, 
                        direct_sigmas = False, linewidth = None):
        """
        Return flux at 50% for the input positions
        most of this is cut and paste from 
        `hetdex_tools/get_spec.py` and 
        the old sensitivity_cube:SensitivityCube 
        class.

        Parameters
        ----------
        ra, dec : array 
            right ascension and dec in degrees
        wave : array 
            wavelength in Angstroms. If None,
            then return flux limits for
            all wave bins.
        sncut : float
            cut in detection significance 
            that defines this catalogue
        direct_sigmas : bool
            return the noise values directly
            without passing them through
            the noise to 50% completeness 
            flux
        linewidth : array
            optionally pass the linewidth of
            the source (in AA) to activate the linewidth
            dependent part of the completeness
            model (default = None).


        Returns
        -------
        f50s : array
            50% completeness. If outside
            of cube return 999. If None
            was passed for wave this is 
            a 2D array of ra, dec and
            all wavelengths
        norm_all : array
            the aperture corrections
        mask : array
            Only returned if `wave=None`. This
            mask is True where the ra/dec
            positions passed are in good
            regions of data
        amp : array
            Only returned in `return_amp=True`,
            it's an array of amplifier information
            for the closest fiber to each source 
        """

        logger = logging.getLogger(name="ShotSensitivity")

        try:
            [x for x in ra]
        except TypeError:
            ra = array([ra])        
            dec = array([dec])       
            wave = array([wave])       

        coords = SkyCoord(ra=ra, dec=dec, unit="deg")
        wave_rect = self.extractor.get_wave()
        pixsize_aa = wave_rect[1] - wave_rect[0]

        # This will give 999 once the noise is scaled suitably
        badval = 999*1e17/pixsize_aa

        if type(wave) != type(None):
            wave_passed = True
        else:
            wave_passed = False
            convolution_filter = ones(2*self.wavenpix + 1)        
            mask = True*ones(len(coords), dtype=int)
            
        noise = []
        
        info_results = self.extractor.get_fiberinfo_for_coords(
                                coords,
                                radius=self.rad,
                                ffsky=self.ffsky,
                                return_fiber_info=True,
                                fiber_lower_limit=2, 
                                verbose=False
                                )

        id_, aseps, aifux, aifuy, axc, ayc, ara, adec, adata, aerror, afmask, afiberid, \
                    amultiframe = info_results
               
 
        I = None
        fac = None
        norm_all = []
        amp = []       
        nan_fib_mask = []

        for i, c in enumerate(coords):
                
            sel = (id_ == i)

            if type(wave) != type(None):
                logger.debug("Running on source {:f} {:f} {:f}".format(ra[i], dec[i], wave[i]))
            else:
                logger.debug("Running on position {:f} {:f}".format(ra[i], dec[i]))

            logger.debug("Found {:d} fibers".format(sum(sel)))

            if sum(sel) > 0:
                
                # fiber properties    
                xc = axc[sel][0]
                yc = ayc[sel][0]
                ifux = aifux[sel]
                ifuy = aifuy[sel]
                data = adata[sel]
                error = aerror[sel]
                fmask = afmask[sel]
                fiberid = afiberid[sel]
                multiframe = amultiframe[sel]
                seps = aseps[sel]

                # test if the data are all zero
                zero_flag = True
                if max((fmask*data).ravel()) < 1e-30:
                    zero_flag = False

                iclosest = argmin(seps)

                amp.append(fiberid[iclosest])

                if len(self.bad_amps) > 0:
                    amp_flag = amp_flag_from_fiberid(fiberid[iclosest], 
                                                         self.bad_amps)
                else:
                    amp_flag = True
                    
                # XXX Could be faster - reloads the file every run
                meteor_flag = meteor_flag_from_coords(c, self.shotid)

                if not (amp_flag and meteor_flag and zero_flag):
                    logger.debug("The data here are bad, position is masked")
                    if wave_passed:
                        noise.append(badval)
                        norm_all.append(1.0)
                        # value doesn't matter as in amp flag
                        nan_fib_mask.append(True)
                        continue
                    else:
                        mask[i] = False
                        
                weights, I, fac = self.extractor.build_weights(xc, yc, ifux, ifuy, self.moffat, 
                                                               I=I, fac=fac, return_I_fac = True)
                

                # (See Greg Zeimann's Remedy code)
                # normalized in the fiber direction
                norm = sum(weights, axis=0) 
                weights = weights/norm

                result = self.extractor.get_spectrum(data, error, fmask, weights,
                                                     remove_low_weights = False,
                                                     sclean_bad = self.sclean_bad,
                                                     return_scleaned_mask = True)
                
                spectrum_aper, spectrum_aper_error, scleaned = [res for res in result] 


                #for w, s in zip(wave_rect, spectrum_aper/norm):
                #    print(w, s)
 
                if wave_passed:
                    
                    index = where(wave_rect >= wave[i])[0][0]
                    ilo = index - self.wavenpix
                    ihi = index + self.wavenpix + 1

                    # Output lots of information for very detailed debugging
                    if logger.getEffectiveLevel() == logging.DEBUG: 
                        logger.debug("Table of fibers:")
                        logger.debug("# fiberid    wave_index ifux ifuy  weight     noise")
                        for fibidx, fid in enumerate(fiberid):
                            for wi, (tw, tnoise) in enumerate(zip((weights*norm)[fibidx, ilo:ihi], 
                                                              error[fibidx, ilo:ihi]), 
                                                              ilo): 
                                logger.debug("{:s} {:d} {:f} {:f} {:f} {:f}".format(fid, wi, ifux[fibidx], 
                                                                                    ifuy[fibidx], tw, tnoise))


                    # Mask source if bad values within the central 3 wavebins
                    nan_fib = bad_central_mask(weights*norm, scleaned, index)   
                    nan_fib_mask.append(nan_fib)

                    # Account for NaN and masked spectral bins
                    bad = isnan(spectrum_aper_error[ilo:ihi])
                    goodfrac = 1.0 - sum(bad)/(ihi - ilo)


                    if all(isnan(spectrum_aper_error[ilo:ihi])):
                        sum_sq = badval
                    else:
                        sum_sq = \
                            sqrt(nansum(square(spectrum_aper_error[ilo:ihi]/goodfrac)))


                    norm_all.append(nansum(norm[ilo:ihi])/len(norm[ilo:ihi]))

                    noise.append(sum_sq)
                else:
                    convolved_variance = convolve(convolution_filter,
                                                  square(spectrum_aper_error),
                                                  mode='same')
                    std = sqrt(convolved_variance)

                    # XXX This will be slow might be able to speed up
                    # Mask wavelengths with too many bad pixels
                    wunorm = weights*norm
                    for index in range(len(convolved_variance)):
                        if not bad_central_mask(wunorm, scleaned, index):
                            std[index] = badval

                    noise.append(std)
                    norm_all.append(norm)
                    
            else:
                if wave_passed:
                    noise.append(badval)
                    norm_all.append(1.0)
                    amp.append("000")
                    nan_fib_mask.append(True)
                else:
                    noise.append(badval*ones(len(wave_rect)))
                    norm_all.append(ones(len(wave_rect)))
                    amp.append("000")
                    mask[i] = False


                   
        # Apply the galaxy mask     
        gal_mask = ones(len(coords), dtype=int)
        for gal_region in self.gal_regions:
            dummy_wcs = create_dummy_wcs(gal_region.center,
                                         imsize=2*gal_region.height)
            # zero if near galaxy
            gal_mask = gal_mask & invert(gal_region.contains(coords, dummy_wcs))

        noise = array(noise)
        snoise = pixsize_aa*1e-17*noise

        if wave_passed:

            bad = (gal_mask < 0.5) | (snoise > 998) | isnan(snoise) | invert(nan_fib_mask)

            normnoise = snoise/norm_all

            if not direct_sigmas:
                normnoise = self.f50_from_noise(normnoise, wave, sncut,
                                                linewidth = linewidth)

            
            normnoise[bad] = 999.

            return normnoise, amp, norm_all

        else:
            mask[gal_mask < 0.5] = False
 
            if self.badshot:
                mask[:] = False

            bad = snoise > 998
            normnoise = snoise/norm_all

            if not direct_sigmas:
                normnoise = self.f50_from_noise(normnoise, wave, sncut,
                                                linewidth = linewidth)

            normnoise[bad] = 999

            return normnoise, mask, amp, norm_all

    def return_completeness(self, flux, ra, dec, lambda_, sncut, f50s=None,
                            linewidth = None):
        """
        Return completeness at a 3D position as an array. 
        If for whatever reason the completeness is NaN, it's
        replaced by 0.0. This is cut and paste from
        sensitivity_cube:SensitivityCube

        Parameters
        ----------
        flux : array
            fluxes of objects
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstrom
        sncut : float
            the detection significance (S/N) cut
            applied to the data
        f50s : array (optional)
            optional array of precomputed
            50% completeness fluxes. Otherwise
            the method will compute them itself
            from the ra/dec/linewidth (default:None)
        linewidth : array (optional)
            optionally pass the linewidth of
            the source (in AA) to activate the linewidth
            dependent part of the completeness
            model (default = None). Only does
            anything when you don't pass the f50s
            (default: None)

        Return
        ------
        fracdet : array
            fraction detected 

        Raises
        ------
        WavelengthException :
            Annoys user if they pass
            wavelength outside of
            VIRUS range
        """

        try: 
            if lambda_[0] < 3000.0 or lambda_[0] > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        except TypeError as e:
             if lambda_ < 3000.0 or lambda_ > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        

        if type(f50s) == type(None):
            f50s, norm = self.get_f50(ra, dec, lambda_, sncut, linewidth = linewidth)

        try:
            # to stop bad values breaking interpolation
            bad = (f50s > 998)
            f50s[bad] = 1e-16
            fracdet = self.sinterp(flux, f50s, lambda_, sncut)
            #print(min(flux), max(flux), min(f50s), max(f50s))

            # check to see if we're passing multiple fluxes
            # for one f50 value
            if any(bad):
                logger.debug("There are bad values here to mask")
                if len(f50s) == 1:
                    logger.debug("Just one ra/dec/wave passed.")
                    fracdet[:] = 0.0
                    f50s[:] = 999.0
                else:
                    fracdet[bad] = 0.0
                    f50s[bad] = 999.0

        except IndexError as e:
            print("Interpolation failed!")
            print(min(flux), max(flux), min(f50s), max(f50s))
            print(min(lambda_), max(lambda_))
            raise e

        try:
            fracdet[isnan(fracdet)] = 0.0
        except TypeError:
            if isnan(fracdet):
                fracdet = 0.0

        return fracdet
    
    def close(self):
        """ 
        Close the Extractor object 
        (especially if it has a Shot HDF
        file open)
        
        """
        self.extractor.close()

    def __enter__(self):
        """ Added to support using the `with` statement """
        return self

    def __exit__(self, type_, value, traceback):
        """ Support tidying up after using the `with` statement """
        self.close()


