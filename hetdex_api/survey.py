# -*- coding: utf-8 -*-
"""

Initiates survey class and provides an API to query the survey class based
on astropy coordinates

Added 2020/03/29 = FiberIndex class to query all fibers in the survey

Created on Tue Jan 22 11:02:53 2019
@author: gregz/Erin Mentuch Cooper
"""
from __future__ import print_function

import numpy as np
import tables as tb
import random
import numpy
import copy
import time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, join, hstack
from astropy import wcs

from regions import EllipseSkyRegion, EllipsePixelRegion

import matplotlib

matplotlib.use("agg")

import healpy as hp
from hetdex_api.config import HDRconfig


try:
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr5"


class Survey:
    def __init__(self, survey=LATEST_HDR_NAME):
        """
        Initialize the Survey class for a given data release

        Parameters
        ----------
        survey : string
            Data release you would like to load, i.e., 'hdr1','HDR2'
            This is case insensitive.

        Returns
        -------
        Survey class object with the following attributes
        """

        global config
        config = HDRconfig(survey=survey.lower())

        self.filename = config.surveyh5
        self.hdfile = tb.open_file(self.filename, mode="r")
        colnames = self.hdfile.root.Survey.colnames

        for name in colnames:
            if name == "ra_flag":
                setattr(
                    self,
                    name,
                    getattr(self.hdfile.root.Survey.cols, name)[:].astype(str),
                )
            elif isinstance(getattr(self.hdfile.root.Survey.cols, name)[0], np.bytes_):
                setattr(
                    self,
                    name,
                    getattr(self.hdfile.root.Survey.cols, name)[:].astype(str),
                )
            else:
                setattr(self, name, getattr(self.hdfile.root.Survey.cols, name)[:])

        # set the SkyCoords
        self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame="icrs")

        # append flux limits
        if survey == "hdr2.1":
            flim = Table.read(
                config.flim_avg,
                format="ascii",
                names=["datevobs", "col2", "fluxlimit_4540"],
            )
            fluxlimit = []
            for datevobs in self.datevobs:
                sel = flim["datevobs"] == datevobs
                if np.sum(sel) == 1:
                    fluxlimit.extend(flim["fluxlimit_4540"][sel])
                elif np.sum(sel) > 1:
                    print("Found two fluxlimits for ", datevobs)
                    fluxlimit.extend(flim["fluxlimit_4540"][sel][0])
                else:
                    fluxlimit.append(np.nan)

            self.fluxlimit_4540 = np.array(fluxlimit)
        else:
            flim = Table.read( config.flim_avg, format='ascii')

            fluxlimit = []

            for shotid in self.shotid:
                sel = flim["shotid"] == shotid
                if np.sum(sel) == 1:
                    fluxlimit.extend(flim["f1sigma_biweight"][sel])
                elif np.sum(sel) > 1:
                    print("Found two fluxlimits for ", shotid)
                    fluxlimit.extend(flim["f1sigma_biweight"][sel][0])
                else:
                    fluxlimit.append(np.nan)

            self.fluxlimit_4600 = np.array(fluxlimit).astype(np.float32)
            
        # added 2024-09-19 by EMC
        # manually change field entry for 20190802012 and 20190803012 incorrectly given a fall objid

        for s in [20190802012, 20190803012]:
            self.field[self.shotid == s] = "dex-spring"

        # added 2025-02-14 by EMC                                                                                        
        # manually change field entry for 20190802012 and 20190803012 incorrectly given a fall objid                      
        for s in [20231014013, 20231017015, 20231104011, 20231104015, 20231108013, 20231204011, 20231206015]:
            self.field[self.shotid == s] = "dex-fall"

        self.field[ self.field == 'egs'] = 'dex-spring'

    def __getitem__(self, indx):
        """
        This allows for slicing of the survey class
        object so that a mask can be applied to
        every attribute automatically by:

        Example
        -------
        survey_sliced = survey[indx]
        """

        p = copy.copy(self)
        attrnames = self.__dict__.keys()
        for attrname in attrnames:
            try:
                setattr(p, attrname, getattr(self, attrname)[indx])
            except TypeError:
                setattr(p, attrname, getattr(self, attrname))
        return p

    def slice(self):
        """
        This will slice the survey class upon
        initialization to remove
        any bad shots from the survey so that
        only one class variable is used.

        Example
        -------
        survey = Survey('hdr1').slice()

        """

        maskshots = self.remove_shots()
        return self[maskshots]

    def remove_shots(self, badshot=True, badtp=True, other=True):
        """
        Function to remove shots that should be exluded from most analysis.

        Parameters
        ----------
        self
           Survey class object

        Returns
        -------
        ind_good
           Boolean mask of good shots

        Examples
        --------
        S = Survey('hdr2')
        ind_good_shots = S.remove_shots()
        S_good = S[ind_good_shots]

        """
        global config

        mask = np.zeros(np.size(self.shotid), dtype=bool)
        badshots = np.loadtxt(config.badshot, dtype=int)

        if badshot:
            for shot in badshots:
                maskshot = self.shotid == shot
                mask = mask | maskshot

        mask = np.invert(mask)
        
        if badtp:
            mask = mask * (self.response_4540 >= 0.08)

        if other:
            mask = mask * (self.field != 'other')
            
        notvalid = self.shotid > 20170000
        mask = mask * notvalid

        return mask

    def get_shotlist(self, coords, radius=None, width=None, height=None):
        """
        Returns a list of shots for defined region of sky.

        Parameters
        ----------
        self
            Survey Class object
        coords
            astropy coordinate object
        radius
            an astropy Quantity object, or a string that can be parsed into
            one.  e.g., '1 degree' or 1*u.degree.  If radius is specified,
            the shape is assumed to be a circle
        width
            in degress.  Specifies the edge length of a square box. If a radius
            is not specified a width and height of a region must be given.
        height
            in degrees.  Specifies the height of a rectangular box.  Must be
            passed with width.

        Examples
        --------
        S = Survey('hdr1')
        coords = SkyCoord(150.025513 * u.deg, 2.087767 * u.deg, frame='icrs')
        shotlist = S.get_shotlist(coords, radius=0.5*u.deg)

        or for a rectangular region around the point defined in coords

        shotlist = S.get_shotlist(coords, width=0.5, height=0.1)
        """

        try:
            self.coords
        except:
            self.coords = SkyCoord(
                self.ra * u.degree, self.dec * u.degree, frame="icrs"
            )

        if radius is not None:
            try:
                idx = self.coords.separation(coords) < radius
            except:
                print("Assuming radius in degrees")
                idx = self.coords.separation(coords) < radius * u.degree
        else:
            try:
                idx1 = abs(self.ra - coords.ra.value) < width / 2.0
                idx2 = abs(self.dec - coords.dec.value) < height / 2.0
                idx = idx1 * idx2
            except:
                print("Provide both width and height of sky region in degrees.")

        return self.shotid[idx]

    def return_astropy_table(self, return_good=True):
        """
        Function to return an astropy table that is machine readable

        Parameters
        ----------
        self
            Survey Class object
        return_good
            Boolean flag to only exclude shots in config.badshot

        Returns
        -------
        survey_table
            astropy table of the Survey class object
        """

        survey_table = Table(self.hdfile.root.Survey.read())

        # added 2025-05-17 by EMC                                                                                        
        for s in [20190802012, 20190803012]:
            survey_table['field'][survey_table['shotid'] == s] = "dex-spring"
        # added 2025-05-17 by EMC                                                                      
        for s in [20231014013, 20231017015, 20231104011, 20231104015, 20231108013, 20231204011, 20231206015]:
            survey_table['field'][survey_table['shotid'] == s] = "dex-fall"

        survey_table['field'][ survey_table['field']=='egs'] = 'dex-spring'

        good_shots = self.remove_shots()

        survey_table["shot_flag"] = good_shots

        try:
            survey_table["mjd"] = self.mjd[:, 0]
        except:
            pass

        try:
            survey_table["exptime"] = np.mean(self.exptime, axis=1)
        except:
            pass

        try:
            survey_table["fluxlimit_4540"] = self.fluxlimit_4540
        except:
            pass
        try:
            survey_table["fluxlimit_4600"] = self.fluxlimit_4600
        except:
            pass


        for col in survey_table.colnames:
            try:
                if np.shape(survey_table[col])[1] == 3:
                    survey_table.remove_column(col)
            except:
                pass

        if return_good:
            return survey_table[good_shots]
        else:
            return survey_table

    def close(self):
        """
        Close the HDF5 file when you are done using
        it to release anything that might be in memory

        Example
        -------

        S.close()
        """

        self.hdfile.close()


class FiberIndex:
    def __init__(
        self,
        survey=LATEST_HDR_NAME,
        load_fiber_table=False,
        loadall=False,
        keep_mask=True,
        shot_h5=None
    ):
        """
        Initialize the Fiber class for a given data release

        Parameters
        ----------
        survey : string
            Data release you would like to load, i.e., 'hdr1','HDR2'
            This is case insensitive.
        load_fiber_table : bool
            Option to read in all fibers. This takes about a minute
            and will use a lot of memory.

        Returns
        -------
        FiberIndex class object
        """

        #will need to know if this is None or set
        self.shot_h5 = shot_h5

        if self.shot_h5 is None: #the usual HETDEX, survey level
            self.survey = survey

            if self.survey == "hdr1":
                print("Sorry there is no FiberIndex for hdr1")
                return None

            global config
            config = HDRconfig(survey=survey.lower())

            self.filename = config.fiberindexh5
            self.hdfile = tb.open_file(self.filename, mode="r")
            self.fiber_table = None
            self.fibermaskh5 = None

            if keep_mask:
                try:
                    self.fibermaskh5 = tb.open_file(config.fibermaskh5, "r")
                except:
                    print("Could not find fiber mask file in {}".format(config.fibermaskh5))
                    self.fibermaskh5 = None

            if load_fiber_table:
                self.fiber_table = Table(self.hdfile.root.FiberIndex.read())
                self.coords = SkyCoord(
                    self.fiber_table["ra"] * u.degree,
                    self.fiber_table["dec"] * u.degree,
                    frame="icrs",
                )

                # add masking info if found
                if self.fibermaskh5 is not None:
                    if keep_mask:
                        self.mask_table = Table(self.fibermaskh5.root.Flags.read())
                        self.fiber_table = hstack([self.fiber_table, self.mask_table])

                        for row in self.fiber_table:
                            if row["fiber_id_1"] == row["fiber_id_2"]:
                                continue
                            else:
                                print(
                                    "Something is wrong. Mismatcheded fiber:{} and {}".format(
                                        row["fiber_id_1"], row["fiber_id_2"]
                                    )
                                )
                        self.fiber_table.rename_column("fiber_id_1", "fiber_id")
                        self.fiber_table.remove_column("fiber_id_2")

        else: #shot_h5 was specified
            if keep_mask:
                try:
                    self.fibermaskh5 = tb.open_file(self.shot_h5, "r")
                    self.mask_table = Table(self.fibermaskh5.root.FiberIndex.read()) #the "Flags" are here
                    self.fiber_table = self.mask_table #for compatibility, and these are now the same
                except:
                    print("Could not find fiber mask file in {}".format(shot_h5))
                    self.fibermaskh5 = None

            self.hdfile = tb.open_file(self.shot_h5, mode="r")


    def return_shot(self, shotid_obs, return_index=False, return_flags=True):
        """
        Function to return the fiber index and mask for a single shot

        Parameters
        ----------
        self
            the FiberIndex class for a specific survey
        shotid
            Specific shotid (dtype=int) you want
        return_index: bool
            Option to return row index values for slicing. Default is False.
        return_flags: bool
            Option to return mask info. Default is True

        Returns
        -------
        fiber_table: astropy table
            An astropy table of Fiber infomation in queried aperture
        table_index: optional
            an optional array of row coordinates corresponding to the
            retrieved fiber table
        """

        tab_idx = self.hdfile.root.FiberIndex.get_where_list("(shotid == shotid_obs)")

        if len(tab_idx) == 0:
            print("Something is wrong, no rows returned from Fiber Index")

        tab = Table(self.hdfile.root.FiberIndex.read_coordinates(tab_idx))

        if return_flags:
            if self.fibermaskh5 is None:
                print("No fiber mask file found")
            else:
                if self.shot_h5 is None:
                    mask_table = Table(
                        self.fibermaskh5.root.Flags.read_coordinates(tab_idx)
                    )
                else:
                    mask_table = Table(
                        self.fibermaskh5.root.FiberIndex.read_coordinates(tab_idx)
                    )

            fiber_table = hstack([tab, mask_table])

            # quick check to make sure fibers match
            for row in fiber_table[
                np.random.random_integers(0, high=len(fiber_table), size=50)
            ]:
                if row["fiber_id_1"] == row["fiber_id_2"]:
                    continue
                else:
                    print(
                        "Something is wrong. Mismatcheded fiber:{} and {}".format(
                            row["fiber_id_1"], row["fiber_id_2"]
                        )
                    )

            fiber_table.rename_column("fiber_id_1", "fiber_id")
            fiber_table.remove_column("fiber_id_2")
        else:
            fiber_table = tab

        if return_index:
            try:
                return fiber_table, tab_idx
            except TypeError:
                return None, None
        else:
            return fiber_table

    def query_region(
        self,
        coords,
        radius=None,
        width=None,
        shotid=None,
        return_index=False,
        return_flags=True,
    ):
        """
        Function to retrieve the indexes of the FiberIndex table
        for a specific region

        Parameters
        ----------
        self
            the FiberIndex class for a specific survey
        coords
            center coordinate you want to search. This should
            be an astropy SkyCoord object
        width
            width of rectangle search if provided. An astropy quantity.
            This will create a cutout of fibers width x width centered at coords.
        radius
            radius you want to search. An astropy quantity object
        shotid
            Specific shotid (dtype=int) you want
        return_index: bool
            Option to return row index values for slicing. Default
            is False
        return_flags: bool
            Option to return mask info. Default is True

        Returns
        -------
        table: astropy table
            An astropy table of Fiber infomation in queried aperture
        table_index: optional
            an optional array of row coordinates corresponding to the
            retrieved fiber table
        """

        Nside = 2**15

        ra_obj = coords.ra.deg
        dec_obj = coords.dec.deg

        if width is not None:
            if radius is not None:
                print('Extracting square cutout of width={:3.2f}. Ignoring radius parameter')
            # retrieve square cutout region
            width_deg = width.to(u.deg).value
            ra_low = coords.ra.deg - 0.5*width_deg
            ra_high = coords.ra.deg + 0.5*width_deg

            table_index = self.hdfile.root.FiberIndex.get_where_list( '(ra >= ra_low) & (ra <= ra_high)')    

            seltab = Table( self.hdfile.root.FiberIndex.read_coordinates( table_index))

            # now select in declination
            idx = np.abs( seltab['dec'] - dec_obj) <= 0.5*width_deg
                          
            selected_index = np.array(table_index)[idx]
                
        else:
            if radius is None:
                radius = 3.5 * u.arcsec

            # if search radius is greater than 1 arcmin, use RA indexing to down-select
            if radius >= 1*u.arcmin:

                rad_deg = radius.to(u.deg).value

                ra_low = coords.ra.deg - rad_deg/np.cos( coords.dec)
                ra_high = coords.ra.deg + rad_deg/np.cos( coords.dec)
                table_index = self.hdfile.root.FiberIndex.get_where_list( '(ra >= ra_low) & (ra <= ra_high)')
                seltab =  Table( self.hdfile.root.FiberIndex.read_coordinates( table_index))
                fib_coords = SkyCoord( ra=seltab['ra']*u.deg, dec=seltab['dec']*u.deg)
                idx = fib_coords.separation( coords) <= radius
                selected_index = np.array(table_index)[idx]
                
            else:
                ra_sep = radius.to(u.degree).value + 4.0 / 3600.0

                vec = hp.ang2vec(ra_obj, dec_obj, lonlat=True)

                pix_region = hp.query_disc(Nside, vec, (ra_sep * np.pi / 180))

                table_index = []

                for hpix in pix_region:
                    if shotid:
                        h_tab, h_tab_index = self.get_fib_from_hp(
                            hpix, shotid=shotid, return_index=True
                        )
                    else:
                        h_tab, h_tab_index = self.get_fib_from_hp(hpix, return_index=True)
                    
                    table_index.extend(h_tab_index)

                seltab = Table( self.hdfile.root.FiberIndex.read_coordinates( table_index))

                fibcoords = SkyCoord(
                    seltab["ra"] * u.degree, seltab["dec"] * u.degree, frame="icrs"
                )

                idx = coords.separation(fibcoords) <= radius

                selected_index = np.array(table_index)[idx]
            
        if return_flags:
            if self.fibermaskh5 is None:
                print("No fiber mask file found")
            else:
                if self.shot_h5 is None:
                    mask_table = Table(
                        self.fibermaskh5.root.Flags.read_coordinates(selected_index)
                    )
                else:
                    mask_table = Table(
                        self.fibermaskh5.root.FiberIndex.read_coordinates(selected_index)
                    )

            fiber_table = hstack([seltab[idx], mask_table])
            # check fibers match
            for row in fiber_table:
                if row["fiber_id_1"] == row["fiber_id_2"]:
                    continue
                else:
                    print(
                        "Something is wrong. Mismatcheded fiber:{} and {}".format(
                            row["fiber_id_1"], row["fiber_id_2"]
                        )
                    )
            fiber_table.rename_column("fiber_id_1", "fiber_id")
            fiber_table.remove_column("fiber_id_2")
        else:
            fiber_table = seltab[idx]

        if return_index:
            try:
                return fiber_table, selected_index
            except TypeError:
                return None, None
        else:
            return fiber_table

    def get_fib_from_hp(self, hp, shotid=None, astropy=True, return_index=False, return_flags=True):
        """
        Return rows with corresponding healpix value

        Parameters
        ----------
        hp:int
           healpix integer you want to search
        shotid:int
           optional shotid to search. Default is none.
        astropy: bool
           returned table is an astropy table object. Defaults to True.
        return_index: bool
           option to return fiber row index values for additional slicing
        return_flags: bool
           will return latest mask flags. Defaults to True

        Returns
        -------
        table
            table with fiber information including coordinates, healpix, fiber_id name
            and more
        table_index
            row indices to access h5 table coordinates
        """

        if shotid is None:
            tab_idx = self.hdfile.root.FiberIndex.get_where_list("(healpix == hp)")
        else:
            tab_hp = self.hdfile.root.FiberIndex.get_where_list("(healpix == hp)")
            tab_shotid = np.where(
                self.hdfile.root.FiberIndex.read_coordinates(tab_hp)["shotid"] == shotid
            )
            tab_idx = tab_hp[tab_shotid]

        tab = self.hdfile.root.FiberIndex.read_coordinates(tab_idx)

        if return_flags:
            if self.fibermaskh5 is None:
                print("No fiber mask file found")
            else:
                if self.shot_h5 is None:
                    mask_table = Table(
                        self.fibermaskh5.root.Flags.read_coordinates(tab_idx)
                    )
                else:
                    mask_table = Table(
                        self.fibermaskh5.root.FiberIndex.read_coordinates(tab_idx)
                    )
                # force an astropy table format
                astropy=True
            
        if astropy:
            if return_index:
                if return_flags:
                    return Table(tab), tab_idx
                else:
                    return hstack( [Table(tab), mask_table]), tab_idx
            else:
                if return_flags:
                    return hstack( [Table(tab), mask_table])
                else:
                    return Table(tab)
        else:
            if return_index:
                return tab, tab_idx
            else:
                return tab

    def get_closest_fiberid(self, coords, shotid=None, maxdistance=None):
        """
        Function to retrieve the closest fiberid in a shot

        Parameters
        ----------
        self
            the FiberIndex class for a specific survey
        coords
            center coordinate you want to search. This should
            be an astropy SkyCoord object
        shotid
            Specific shotid (dtype=int) you want
        maxdistance
            The max distance (an angle quantity) you want to search for a nearby fiber.
            Default is 8.*u.arcsec

        Returns
        -------
        fiberid
            unique fiber identifier for the shot

        """

        # start searching at small radius to search more efficiently
        search = 2.0 * u.arcsec

        if maxdistance is None:
            maxdistance = 8.0 * u.arcsec

        while search <= maxdistance:
            fiber_table = self.query_region(coords, radius=search, shotid=shotid)
            search = search + 4.0 * u.arcsec
            if np.size(fiber_table) > 0:
                break

        if np.size(fiber_table) > 0:
            fibcoords = SkyCoord(
                fiber_table["ra"] * u.degree,
                fiber_table["dec"] * u.degree,
                frame="icrs",
            )

            idx = np.argmin(coords.separation(fibcoords))
            fiberid = fiber_table["fiber_id"][idx]

            return fiber_table["fiber_id"][idx]
        else:
            return None

    def get_amp_flag(self):
        """
        Generate amp flag for each fiber in the survey

        Parameters
        ----------
        FiberIndex Class
        
        Returns
        -------
        fiber_id: str
        unique fiber identifier string
        amp_flag: bool
        True if fiber is on a good quality amplifier
        """
        global config

        print("Adding amplifier flags")
        t0 = time.time()

        badamps = Table.read(config.badamp)
        self.fiber_table["row_index"] = np.arange(0, len(self.fiber_table))
        join_tab = join(
            self.fiber_table, badamps, keys=["shotid", "multiframe"], join_type="left"
        )
        join_tab.rename_column("flag", "amp_flag")
        t1 = time.time()

        join_tab.sort("row_index")

        # quick check to make sure columns match

        for idx in np.random.random_integers(0, high=len(self.fiber_table), size=500):
            if self.fiber_table["fiber_id"][idx] != join_tab["fiber_id"][idx]:
                print("Something went wrong. fiber_id columns don't match")

        print("Done adding amplifier flags in {:4.3} minutes".format((t1 - t0) / 60))
        
        return np.array(join_tab["amp_flag"], dtype=bool)

    def get_shot_flag(self):
        """
        Generate shot flag for each fiber in the survey

        Parameters
        ----------
        FiberIndex Class

        Returns
        -------
        shot_flag: bool
            True if fiber is on a good shot

        """
        global config

        print("Adding shot flag")
        t0 = time.time()

        badshots = np.loadtxt(config.badshot, dtype=int)

        shot_flag = np.ones_like(self.fiber_table["fiber_id"], dtype=bool)
        for badshot in badshots:
            sel = self.fiber_table["shotid"] == badshot
            shot_flag[sel] = False
        t1 = time.time()
        print("Shot flags added in {:4.2f}".format((t1 - t0) / 60))

        return shot_flag

    def get_badfiber_flag(self):
        """
        Flag bad fibers

        Parameters
        ----------
        FiberIndex Class

        Returns
        -------
        badfiber_mask
        """

        global config

        print("Adding bad fiber flag")

        t0 = time.time()
        badfib = Table.read(config.badfib, format="ascii")

        badfib_transient = Table.read(
            config.badfib_transient,
            format="ascii",
            names=["multiframe", "fibnum", "date_start", "date_end"],
        )

        mask = np.ones(len(self.fiber_table), dtype=bool)

        for row in badfib:
            sel_fib = (self.fiber_table["multiframe"] == row["multiframe"]) * (
                self.fiber_table["fibnum"] == row["fibnum"]
            )

            if np.sum(sel_fib) > 0:
                # print('Adding badfiber flag', row['multiframe'], row['fibnum'])
                mask = mask * np.invert(sel_fib)

        for row in badfib_transient:
            sel_fib = (self.fiber_table["multiframe"] == row["multiframe"]) * (
                self.fiber_table["fibnum"] == row["fibnum"]
            )
            sel_date = (self.fiber_table["date"] >= row["date_start"]) * (
                self.fiber_table["date"] <= row["date_end"]
            )
            if np.sum(sel_fib * sel_date) > 0:
                mask = mask * np.invert(sel_fib)

        t1 = time.time()
        print("Bad fiber flags added in {:4.2f}".format((t1 - t0) / 60))

        return mask

    def get_throughput_flag(self, throughput_threshold=0.08):
        """
        Generate throughput flag for each fiber in the survey

        Parameters
        ----------
        FiberIndex Class

        Returns
        -------
        throughput_flag: bool
            True if fiber observation is above throughput_threshold
        """
        global config

        print("Adding throughput flag")

        t0 = time.time()

        S = Survey(self.survey)

        sel_lowtp = S.response_4540 < throughput_threshold

        throughput_flag = np.ones_like(self.fiber_table["fiber_id"], dtype=bool)

        for badshot in S.shotid[sel_lowtp]:
            sel = self.fiber_table["shotid"] == badshot
            throughput_flag[sel] = False
        t1 = time.time()
        print("Throughput flags added in {:4.2f}".format((t1 - t0) / 60))

        return throughput_flag

    def get_gal_flag(self, d25scale=3.0):
        """
        Returns boolean mask with detections landing within
        galaxy defined by d25scale flagged as False.

        Based on check_all_large_gal from hetdex_tools/galmask.py
        written by John Feldmeier
        """

        global config
        t0 = time.time()

        S = Survey(self.survey)
        galaxy_cat = Table.read(config.rc3cat, format="ascii")

        mask = np.ones(len(self.fiber_table), dtype=bool)
        # Loop over each galaxy

        for idx in np.arange(len(galaxy_cat)):
            row = galaxy_cat[idx]
            gal_coord = SkyCoord(row["Coords"], frame="icrs")
            rlimit = 1.1 * d25scale * row["SemiMajorAxis"] * u.arcmin

            shots = S.get_shotlist(gal_coord, radius=rlimit)
            if len(shots) == 0:
                continue
            down_select = self.coords.separation(gal_coord) < rlimit
            if np.any(down_select):
                galregion = create_gal_ellipse(
                    galaxy_cat, row_index=idx, d25scale=d25scale
                )

                dummy_wcs = create_dummy_wcs(
                    galregion.center, imsize=2 * galregion.height
                )
                galflag = np.zeros(len(self.fiber_table), dtype=bool)
                galflag[down_select] = galregion.contains(
                    self.coords[down_select], dummy_wcs
                )

                mask = mask * np.invert(galflag)

        t1 = time.time()
        S.close()
        print("Galaxy mask array generated in {:3.2f} minutes".format((t1 - t0) / 60))

        return mask

    def get_meteor_flag(self, streaksize=None):
        """
        Returns boolean mask with detections landing on meteor
        streaks masked. Use np.invert(mask) to find meteors
        """

        if streaksize is None:
            streaksize = 12.0 * u.arcsec

        t0 = time.time()
        global config

        S = Survey(self.survey)
        # meteors are found with +/- X arcsec of the line DEC=a+RA*b in this file

        met_tab = Table.read(config.meteor, format="ascii")

        mask = np.ones(len(self.fiber_table), dtype=bool)

        for row in met_tab:
            sel_shot = np.where(
                (self.fiber_table["shotid"] == row["shotid"])
                & (np.isfinite(self.fiber_table["ra"]))
                & (np.isfinite(self.fiber_table["dec"]))
            )[0]

            if np.sum(sel_shot) > 0:
                a = row["a"]
                b = row["b"]

                sel_shot_survey = S.shotid == row["shotid"]
                shot_coords = S.coords[sel_shot_survey]

                # add in special handling for vertical streak
                if np.abs(b) < 500:
                    ra_met = shot_coords.ra + np.arange(-1500, 1500, 0.1) * u.arcsec
                else:  # for very large slopes we needed higher resolution
                    ra_met = shot_coords.ra + np.arange(-1500, 1500, 0.001) * u.arcsec

                dec_met = (a + ra_met.deg * b) * u.deg

                sel_met = (dec_met > -90 * u.deg) * (dec_met < 90.0 * u.deg)
                met_coords = SkyCoord(ra=ra_met[sel_met], dec=dec_met[sel_met])

                idxc, idxcatalog, d2d, d3d = met_coords.search_around_sky(
                    self.coords[sel_shot], streaksize
                )
                mask_index = [sel_shot][0][idxc]
                mask[mask_index] = False

        S.close()
        t1 = time.time()

        print("Meteor mask array generated in {:3.2f} minutes".format((t1 - t0) / 60))

        return mask

    def get_satellite_flag(self, streaksize=None):
        """
        Returns boolean mask with detections landing on meteor
        streaks masked. Use np.invert(mask) to find meteors
        """

        if streaksize is None:
            streaksize = 6.0 * u.arcsec

        t0 = time.time()
        global config

        S = Survey(self.survey)
        # satellites are found with +/- X arcsec of the line DEC=a+RA*b in this file

        sat_tab = Table.read(
            config.satellite,
            format="ascii",
            names=["shotid", "expnum", "slope", "intercept"],
        )

        mask = np.ones(len(self.fiber_table), dtype=bool)

        for row in sat_tab:
            sel_shot = np.where(
                (self.fiber_table["shotid"] == row["shotid"])
                & (np.isfinite(self.fiber_table["ra"]))
                & (np.isfinite(self.fiber_table["dec"]))
            )[0]

            if np.sum(sel_shot) > 0:
                a = row["slope"]
                b = row["intercept"]

                sel_shot_survey = S.shotid == row["shotid"]
                shot_coords = S.coords[sel_shot_survey]

                ra_sat = shot_coords.ra + np.arange(-1500, 1500, 0.1) * u.arcsec
                dec_sat = (b + ra_sat.deg * a) * u.deg

                sat_coords = SkyCoord(ra=ra_sat, dec=dec_sat)

                idxc, idxcatalog, d2d, d3d = sat_coords.search_around_sky(
                    self.coords[sel_shot], streaksize
                )
                mask_index = [sel_shot][0][idxc]
                mask[mask_index] = False

        S.close()
        t1 = time.time()

        print(
            "Satellite mask array generated in {:3.2f} minutes".format((t1 - t0) / 60)
        )

        return mask

    def get_fiber_flags(self, coord=None, shotid=None):
        """
        Returns boolean mask values for a coord/shotid combo

        Parameters
        ----------
        coord: SkyCoord object
            sky coordinates
        shotid: int
            shotid

        Returns
        -------
        flag_dict :  dictionary
        1 if good, 0 if flagged bad
        """

        if self.fibermaskh5 is None:
            print("No fiber mask file found")
            return None

        table, table_index = self.query_region(coord, return_index=True, shotid=shotid)
        if self.shot_h5 is None:
            mask_table = Table(self.fibermaskh5.root.Flags.read_coordinates(table_index))
        else:
            mask_table = Table(self.fibermaskh5.root.FiberIndex.read_coordinates(table_index))

        flag_dict = {
            "flag": 1,
            "flag_badamp": 1,
            "flag_badfib": 1,
            "flag_meteor": 1,
            "flag_satellite": 1,
            "flag_largegal": 1,
            "flag_shot": 1,
            "flag_throughput": 1,
        }

        for col in mask_table.columns:
            if col == "fiber_id":
                continue
            else:
                flag_dict[col] = int(np.all(mask_table[col]))

        return flag_dict

    def close(self):
        """
        Close the hdfile when done
        """
        self.hdfile.close()
        if self.fibermaskh5 is not None:
            self.fibermaskh5.close()


def create_gal_ellipse(galaxy_cat, row_index=None, pgcname=None, d25scale=3.0):
    """
    Similar to galmask.py/ellreg but can take a galaxy name as input.

    Create galaxy ellipse region using an input galaxy catalog_table (likely need
    to change this API as it heavily follows the RC3 catalog format)

    Parameters
    ----------
    galaxy_catalog_table: an astropy table
        table of nearby galaxy positions and sizes. Must have a central
        coordinate and the SemiMajorAxis, SemiMinorAxis, Position Angle
        info as table column names
    row_index: int
        a row index in the galaxy catalog to create the ellipse for
    pgcname: str
        the PGCNAME in the RC3 cat. This is a string. eg. "PGC 43255" for NGC 4707
    d25scale: float
        how many times D25 should to scale the region
    """
    if row_index is not None:
        index = row_index
    elif pgcname is not None:
        index = np.where(galaxy_cat["PGC"] == pgcname)[0][0]

    coords = SkyCoord(galaxy_cat["Coords"][index], frame="icrs")

    # The ellipse region uses the major and minor axes, so we have to multiply by
    # two first, before applying any user scaling.

    major = (galaxy_cat["SemiMajorAxis"][index]) * d25scale * u.arcmin
    minor = (galaxy_cat["SemiMinorAxis"][index]) * d25scale * u.arcmin
    pa = (galaxy_cat["PositionAngle"][index]) * u.deg
    ellipse_reg = EllipseSkyRegion(center=coords, height=major, width=minor, angle=pa)

    return ellipse_reg


def create_dummy_wcs(coords, pixscale=None, imsize=None):
    """
    Create a simple fake WCS in order to use the regions subroutine.
    Adapted from John Feldmeiers galmask.py

    Parameters
    ----------
    coords: a SkyCoord object
        center coordinates of WCS
    pixscale: astropy quantity
        pixel scale of WCS in astropy angle quantity units. Will default
        to 0.5 arcsec if not set
    imsize: astropy quantity
        size of WCS in astropy angle quanity units. Will default to
        60 arcsec if not set
    """

    if pixscale is None:
        pixscale = 0.5 * u.arcsec
    if imsize is None:
        imsize = 60.0 * u.arcsec

    gridsize = imsize.to_value("arcsec")
    gridstep = pixscale.to_value("arcsec")

    # Create coordinate center
    ra_cen = coords.ra.deg
    dec_cen = coords.dec.deg

    ndim = int(2 * gridsize / gridstep + 1)
    center = ndim / 2
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.cdelt = [-gridstep / gridsize, gridstep / gridsize]
    w.array_shape = [ndim, ndim]

    return w
