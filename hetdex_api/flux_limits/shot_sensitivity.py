"""

Unlike the sensitivity cubes, this module
generates flux limit estimates on the fly
from the Fiber class.

AUTHOR: Daniel Farrow

"""

from os.path import isfile, join
from tempfile import NamedTemporaryFile

from numpy import (arange, meshgrid, ones, array, sum, sqrt, square, newaxis,
                   convolve, repeat, loadtxt, where, argmin, invert, logical_not,
                   isnan, newaxis)
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
        applied (default 3.5)
    verbose : bool
        Print information about the flux limit model
        to the screen
    """
    def __init__(self, datevshot, release=None, flim_model=None, rad=3.5, 
                 ffsky=False, wavenpix=3, d25scale=3.5, verbose=True):

        self.conf = HDRconfig()
        self.extractor = Extract()
        self.shotid = int(datevshot.replace("v", ""))
        self.date = datevshot[:8]
        self.rad = rad
        self.ffsky = ffsky
        self.wavenpix = wavenpix
        self.verbose = verbose
        
        if not release:
            self.release = self.conf.LATEST_HDR_NAME
        else:
            self.release = release
            
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
        self.setup_mask(d25scale)
        
        # Set up spectral extraction
        if release == "hdr1":
            fwhm = self.survey.fwhm_moffat[survey_sel][0]
        else:
            fwhm = self.survey.fwhm_virus[survey_sel][0]

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
        # see if this is a bad shot
        badshot = loadtxt(self.conf.badshot, dtype=int)
        badtpshots = loadtxt(self.conf.lowtpshots, dtype=int)
        if (self.shotid in badshot) or (self.shotid in badtpshots):
            print("Shot is in bad. Making mask zero everywhere")
            self.badshot = True
        else:
            self.badshot = False
        
        # set up bad amps
        self.bad_amps = Table.read(self.conf.badamp)
        sel_shot = (self.bad_amps["shotid"] == self.shotid)
        self.bad_amps = self.bad_amps[sel_shot]
        
        # set up galaxy mask
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
            
            ix, iy = meshgrid(arange(0.5, nx + 0.5, 1.0),
                              arange(0.5, ny + 0.5, 1.0))
        

            all_ra, all_dec, junk = scube.wcs.all_pix2world(ix.ravel(), iy.ravel(), 
                                                            [500.], 0)
            noises, mask = self.get_f50(all_ra, all_dec, None, 
                                        direct_sigmas = True)
            sigmas = noises.ravel(order="F").reshape(nz, ny, nx)

            mask = logical_not(mask.reshape(ny, nx))
            mask3d = repeat(mask[newaxis, :, :], sigmas.shape[0], axis=0)
            scube.sigmas = MaskedArray(sigmas, mask=mask3d, fill_value=999.0)
        
        return scube
       
    
    def get_f50(self, ra, dec, wave, direct_sigmas = False):
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

        Returns
        -------
        f50s : array
            50% completeness. If outside
            of cube return 999. If None
            was passed for wave this is 
            a 2D array of ra, dec and
            all wavelengths
        """
 
        try:
            [x for x in ra]
        except TypeError:
            ra = array([ra])        
            dec = array([dec])       
            wave = array([wave])       

        coords = SkyCoord(ra=ra, dec=dec, unit="deg")
        wave_rect = self.extractor.get_wave()
        pixsize_aa = wave_rect[1] - wave_rect[0]

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
                                verbose=self.verbose
                                )

        id_, aseps, aifux, aifuy, axc, ayc, ara, adec, adata, aerror, afmask, afiberid, \
                    amultiframe = info_results
        
        
        I = None
        fac = None
        norm_all = []
        
        for i, c in enumerate(coords):
                
            sel = (id_ == i)
            
            if sum(sel) > 0:
                    
                xc = axc[sel][0]
                yc = ayc[sel][0]
                ifux = aifux[sel]
                ifuy = aifuy[sel]
                ra = ara[sel]
                dec = adec[sel]
                data = adata[sel]
                error = aerror[sel]
                fmask = afmask[sel]
                fiberid = afiberid[sel]
                multiframe = amultiframe[sel]
                seps = aseps[sel]
                    
                if len(self.bad_amps) > 0:
                    iclosest = argmin(seps)
                    amp_flag = amp_flag_from_fiberid(fiberid[iclosest], 
                                                         self.bad_amps)
                else:
                    amp_flag = True
                    
                # XXX Could be faster - reloads the file every run
                meteor_flag = meteor_flag_from_coords(c, self.shotid)
                
                if not (amp_flag and meteor_flag):
                    if wave_passed:
                        noise.append(999.)
                        norm_all.append(1.0)
                        continue
                    else:
                        mask[i] = False
                        
                weights, I, fac = self.extractor.build_weights(xc, yc, ifux, ifuy, self.moffat, 
                                                               I=I, fac=fac, return_I_fac = True)
                

                # See Greg Zeimann's Remedy code
                norm = sum(weights, axis=0)
                weights = weights/norm
                

                result = self.extractor.get_spectrum(data, error, fmask, weights,
                                                     remove_low_weights = False)
                
                spectrum_aper, spectrum_aper_error = [res for res in result] 
                if wave_passed:
                    index = where(wave_rect >= wave[i])[0][0]
                    ilo = index - self.wavenpix 
                    ihi = index + self.wavenpix + 1
                    sum_sq = \
                        sqrt(sum(square(spectrum_aper_error[ilo:ihi])))
                    norm_all.append(sum(norm[ilo:ihi])/len(norm[ilo:ihi]))
                    noise.append(sum_sq)
                else:
                    convolved_variance = convolve(convolution_filter,
                                                  square(spectrum_aper_error),
                                                  mode='same')
                    noise.append(sqrt(convolved_variance))
                    norm_all.append(norm)
                    
            else:
                if wave_passed:
                    noise.append(999.)
                    norm_all.append(1.0)
                else:
                    noise.append(999.*ones(len(wave_rect)))
                    norm_all.append(ones(len(wave_rect)))

                    
        # Apply the galaxy mask     
        gal_mask = ones(len(coords), dtype=int)
        for gal_region in self.gal_regions:
            dummy_wcs = create_dummy_wcs(gal_region.center,
                                         imsize=2*gal_region.height)
            gal_mask = gal_mask & invert(gal_region.contains(coords, dummy_wcs))

        noise = array(noise)
        snoise = pixsize_aa*1e-17*noise

        if wave_passed:
            bad = (gal_mask < 0.5) | (noise > 998)

            if not direct_sigmas:
                snoise = self.f50_from_noise(snoise)

            snoise[bad] = 999.
            return snoise/norm_all

        else:
            mask[(gal_mask < 0.5) & mask] = True
 
            if not direct_sigmas:
                snoise = self.f50_from_noise(snoise)

            return snoise/norm_all, mask

    def return_completeness(self, flux, ra, dec, lambda_, sncut):
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
        

        f50s = self.get_f50(ra, dec, lambda_, sncut)
        fracdet = self.sinterp(flux, f50s, lambda_, sncut)

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


