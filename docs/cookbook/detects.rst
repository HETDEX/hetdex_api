Detections Database
===================

   This notebook demonstrates how to access the HDF5 container for the
   HETDEX line detections database through the API. Querying of the
   database through an interactive GUI follows in Notebook 11 - Querying
   Detections GUI. This database is a catalog of line emission
   detections and their associated 1D, aperture summed, psf-weighted
   spectra. There are three tables contained within this HDF5 file:

   #. Detections - this is the main database of line detection sources.
      It provides the position and central wavelength of each detection
      and corresponding line fluxes. A source detection corresponds to
      an emission line so it is possible to have multiple line
      detections at different wavelengths for a single source. There can
      also be multiple observations of the same line if it has been
      observed in multiple shots or if it is associated with a large
      source.

   #. Fibers - for each source detection, this table lists information
      about each fiber used to extract the flux measurment and weighted
      spectrum. This allows a user to return to the processed data
      products (ie. the shot HDF5 files) to investigate the source
      further.

   #. Spectra - for each source, this table contains arrays of
      wavelength and 1D flux-weighted aperture summed spectral data and
      corresponding errors. Non-calibrated spectra is also provided in
      counts

.. container:: cell code

   .. code:: python

      %matplotlib inline
      import sys
      import os
      import os.path
      import subprocess
      import numpy as np
      import tables as tb
      import matplotlib.pyplot as plt

      from astropy.io import ascii
      from astropy.table import Table, Column
      from astropy.coordinates import SkyCoord
      import astropy.units as u

      from hetdex_api import config
      from hetdex_api.detections import Detections
      from hetdex_api.elixer_widgets import *

.. container:: cell markdown

   .. rubric:: Initiate the API
      :name: initiate-the-api

.. container:: cell markdown

   When you call ``Detections()`` you intiate the Detections Class
   object which takes columns from the Detections Table in the HDF5 file
   and adds them as array attributes to the Detections class object. It
   also converts ra/dec into astropy skycoords in the ``coords``
   attribute, calculates an approximate gband magnitude using the 1D
   spectra and adds elixer probabilities for each detection. If you
   append the call with ``refine()`` then a number of downselections are
   applied to the database to return a more robust list of line
   emitters. ``refine()`` removes spurious detections found in bad amps
   or at the edges of the CCD or in shots that are not deemed
   appropriate for HETDEX analysis. It can also remove all bright
   objects above a specific gband magnitude if desired (default to None
   if no option is given).

.. container:: cell code

   .. code:: python

      detects = Detections('hdr1').refine()
      # detects = Detections('hdr1').refine(gmagcut=22)
      ## or if you want to open the continuum source catalog:
      ## detects = Detections('cont_sources')

.. container:: cell markdown

   Here are a list of attributes built into the Detections class:

.. container:: cell code

   .. code:: python

      detects.__dict__.keys()

   .. container:: output execute_result

      ::

         dict_keys(['survey', 'filename', 'hdfile', 'detectid', 'detectname', 'inputid', 'date', 'obsid', 'ra', 'dec', 'wave', 'wave_err', 'flux', 'flux_err', 'linewidth', 'linewidth_err', 'continuum', 'continuum_err', 'sn', 'sn_err', 'chi2', 'chi2_err', 'multiframe', 'fibnum', 'x_raw', 'y_raw', 'amp', 'ifuid', 'ifuslot', 'shotid', 'specid', 'coords', 'hdfile_elix', 'detectid_elix', 'plae_poii_hetdex', 'aperture_mag', 'plae_poii_aperture', 'aperture_filter', 'mag_match', 'cat_filter', 'plae_poii_cat', 'dec_match', 'dist_match', 'ra_match', 'z_prelim', 'field', 'fwhm', 'fluxlimit_4550', 'throughput', 'n_ifu', 'vis_class', 'gmag', 'plae_poii_hetdex_gmag'])

.. container:: cell markdown

   Most of these columns come directly from the Detections table in the
   detect_hdr1.h5 file. You can find information for each column in the
   HETDEX DR1 Document
   (/work/03946/hetdex/hdr1/doc/Hetdex_Data_Release_1.pdf). Some values
   including field (survey.field), fwhm (guider
   seeing/survey.fwhm_moffat), throughput (survey.response_4540),
   flux_limit (survey.fluxlimit_4550), n_ifu (survey.n_ifu) are copied
   in from the Survey HDF5 file for additional analysis. We also attach
   probabilities for a source to be an LAE compared to an OII emitter
   (plae_poii) calculated by ELiXer using the VIRUS continuum
   ('plae_poii_hetdex'), a matched catalog magnitude ('plae_poii_cat'),
   or from a magnitude measured from associated optical imaging
   ('plae_poii_aperture').

.. container:: cell markdown

   .. rubric:: Querying by sky coordinates
      :name: querying-by-sky-coordinates

.. container:: cell markdown

   Upon initialization of the Detections Class, sky coordinates are
   converted to an Astropy sky coordinates array to allow for easy
   querying:

.. container:: cell code

   .. code:: python

      detects.coords

   .. container:: output execute_result

      ::

         <SkyCoord (ICRS): (ra, dec) in deg
             [(149.77184,  1.948503), (149.78958,  1.900873),
              (149.84552,  1.935277), ..., (189.64575, 50.83938 ),
              (170.32608, 51.549652), (170.32336, 51.58739 )]>

.. container:: cell markdown

   To query a region of the sky, you can use the Detections function
   ``query_by_coords`` which takes an astropy coords objects as an
   argument as well as a radius represented by an astropy quantity. It
   returns a boolean mask to index the Detections class object.

.. container:: cell code

   .. code:: python

      obj_coords = SkyCoord(199.35704 * u.deg, 51.06718 * u.deg, frame='icrs')

.. container:: cell code

   .. code:: python

      maskregion = detects.query_by_coords(obj_coords, 10. * u.arcsec)

.. container:: cell markdown

   The Detections class allows slicing so that a boolean mask applied to
   the class will slice each array attribute accordingly:

.. container:: cell code

   .. code:: python

      detects_in_region = detects[maskregion]
      print(np.size(detects_in_region.detectid))

   .. container:: output stream stdout

      ::

         12

.. container:: cell markdown

   For this example, we have found 12 detections in this region, we can
   examine these via the ELiXer reports using the ``ElixerWidget()``
   class from ``elixer_widgets.py``. To do so we need to save the
   detectid list to examine in the widget.

.. container:: cell code

   .. code:: python

      np.savetxt('detects_obj.txt', detects_in_region.detectid)

.. container:: cell markdown

   You can the run the elixer_widget to scan through the ELiXer reports
   for this object. Use the "Next DetectID" button to scan the list. The
   "DetectID" text widget will give access to all reports interactively
   and scans in increasing single digit increments, but the green Next
   DetectID button will go in order of the ingest list from
   'detects_obj.txt'.

.. container:: cell code

   .. code:: python

      elix_widget = ElixerWidget(detectfile='detects_obj.txt')

   .. container:: output display_data

      .. code:: json

         {"model_id":"c8eb54e7d388405598f37edff94bfc2d","version_major":2,"version_minor":0}

.. container:: cell markdown

   For more information on using the Elixer Widgets GUI go to Notebook
   12. We will discuss team classification efforts there. But for quick
   investigation its helpful to pull the GUI up to just scan through a
   detection list.

.. container:: cell markdown

   .. rubric:: Accessing 1D Spectra
      :name: accessing-1d-spectra

.. container:: cell markdown

   Spectra in counts and flux-calibrated units are stored in the Spectra
   Table of the Detection HDF5 file, it can be accessed directly through
   the Detections class object which stores the detect HDF5 as an
   attribute:

.. container:: cell code

   .. code:: python

      print(detects.hdfile)

   .. container:: output stream stdout

      ::

         /work/03946/hetdex/hdr1/detect/detect_hdr1.h5 (File) 'HDR1 Detections Database'
         Last modif.: 'Mon Apr  1 12:39:46 2019'
         Object Tree: 
         / (RootGroup) 'HDR1 Detections Database'
         /Detections (Table(690868,)) 'HETDEX Line Detection Catalog'
         /Fibers (Table(9397618,)) 'Fiber info for each detection'
         /Spectra (Table(690797,)) '1D Spectra for each Line Detection'

.. container:: cell code

   .. code:: python

      spectra = detects.hdfile.root.Spectra

.. container:: cell markdown

   This is a very large table so its not advised to read it in all at
   once. The columns are:

.. container:: cell code

   .. code:: python

      spectra.cols

   .. container:: output execute_result

      ::

         /Spectra.cols (Cols), 8 columns
           detectid (Column(690797,), int64)
           wave1d (Column(690797, 1036), ('<f4', (1036,)))
           spec1d (Column(690797, 1036), ('<f4', (1036,)))
           spec1d_err (Column(690797, 1036), ('<f4', (1036,)))
           counts1d (Column(690797, 1036), ('<f4', (1036,)))
           counts1d_err (Column(690797, 1036), ('<f4', (1036,)))
           apsum_counts (Column(690797, 1036), ('<f4', (1036,)))
           apsum_counts_err (Column(690797, 1036), ('<f4', (1036,)))

.. container:: cell markdown

   Flux calibrated, psf-weighted 1D spectra can be retrieved via the API
   for a single detectid through the function ``get_spectrum``:

.. container:: cell code

   .. code:: python

      detectid_nice_lae = 1000208773
      spec_table = detects.get_spectrum(detectid_nice_lae) 

.. container:: cell code

   .. code:: python

      detects.plot_spectrum(detectid_nice_lae)

   .. container:: output display_data

      |image0|

.. container:: cell markdown

   or if we want to zoom in on the emission line:

.. container:: cell code

   .. code:: python

      cw = detects.wave[detects.detectid == detectid_nice_lae]
      detects.plot_spectrum(detectid_nice_lae, xlim=(cw-50, cw+50))

   .. container:: output display_data

      |image1|

.. container:: cell markdown

   You can also save the spectrum to a text file. It is automatically
   saved as spec_##detectid##.dat, but you can also use the argument
   ``outfile``

.. container:: cell code

   .. code:: python

      detects.save_spectrum(detectid_nice_lae)
      # or
      # detects.save_spectrum(detectid_nice_lae, outfile='tmp.txt')

.. container:: cell markdown

   .. rubric:: Example: Finding average number of sources per IFU
      :name: example-finding-average-number-of-sources-per-ifu

.. container:: cell markdown

   To reach our survey goal we need to obtain a critical number of
   detections per IFU on average. Here we show how the number of
   detections based on the signal-to-noise requirement.

.. container:: cell code

   .. code:: python

      ndets_ifu = []

      sn_array = np.arange(start = 5, stop = 10, step = 1)

      # only choose detections that lied on good shots

      for sn_i in sn_array:
          sel = (detects.sn > sn_i) * (detects.n_ifu > 0) * (detects.chi2 < 3) * (detects.chi2 >0.1)
          detifu = 1./(detects.n_ifu[sel])
          ndets_ifu.append(np.sum(detifu)/np.size(np.unique(detects.shotid)))
          
      ndets_ifult2 = []

      for sn_i in sn_array:
          sel = (detects.sn > sn_i) * (detects.n_ifu > 0) * (detects.chi2 < 2) * (detects.chi2 >0.1)
          detifu = 1./(detects.n_ifu[sel])
          ndets_ifult2.append(np.sum(detifu)/np.size(np.unique(detects.shotid)))
          
      # only choose detections that lied on good shots
      sel = (detects.throughput > 0.095) * (detects.fwhm < 2.5)
      detects_good_shots = detects[sel]

      ndets_ifu_gs =[]
      for sn_i in sn_array:
          sel = (detects_good_shots.sn > sn_i) * (detects_good_shots.n_ifu > 0) * (detects_good_shots.chi2 < 3) * (detects_good_shots.chi2 >0.1)
          detifu = 1./(detects_good_shots.n_ifu[sel])
          ndets_ifu_gs.append(np.sum(detifu)/np.size(np.unique(detects_good_shots.shotid)))
          
      ndets_ifu_gs_lt2 = []
      for sn_i in sn_array:
          sel = (detects_good_shots.sn > sn_i) * (detects_good_shots.n_ifu > 0) * (detects_good_shots.chi2 < 2) * (detects_good_shots.chi2 >0.1)
          detifu = 1./(detects_good_shots.n_ifu[sel])
          ndets_ifu_gs_lt2.append(np.sum(detifu)/np.size(np.unique(detects_good_shots.shotid)))
          

.. container:: cell code

   .. code:: python

      plt.rcParams.update({'font.size': 18})
      plt.figure(figsize=(9,9))
      plt.scatter(sn_array, ndets_ifu, label='Chi2 < 3')
      plt.scatter(sn_array, ndets_ifu_gs, label='Chi2 < 3, tp > 0.095, fwhm < 2.5')
      plt.scatter(sn_array, ndets_ifult2, label='Chi2 < 2')
      plt.scatter(sn_array, ndets_ifu_gs_lt2, label='Chi2 < 2, tp > 0.095, fwhm < 2.5')
      plt.xlabel('SN')
      plt.ylabel('N detections per IFU')
      plt.legend(fontsize='small')
      plt.savefig('ndetsperifu_vs_sn.png')

   .. container:: output display_data

      |image2|

.. container:: cell markdown

   .. rubric:: Saving to a file
      :name: saving-to-a-file

.. container:: cell markdown

   If you want to just save a subset of columns for a subset of
   detections, use the ``return_astropy_table()`` function to return all
   column attributes of the Detections class into an astropy table which
   you may then save.

.. container:: cell code

   .. code:: python

      detects = Detections('hdr1').refine(gmagcut=21)
      sel = (detects.throughput > 0.09) * (detects.fwhm < 2.6) * (detects.chi2 < 1.6) * (detects.chi2 < 1.1+0.9*(detects.sn-5.2)/(8-5.2)) 
      detects_sel = detects[sel]
      table_sel = detects_sel.return_astropy_table()

.. container:: cell code

   .. code:: python

      ascii.write(table_sel, 'HDR1_source_catalog_20190628.dat', overwrite=True)

.. container:: cell markdown

   .. rubric:: Getting Fiber information for a detection
      :name: getting-fiber-information-for-a-detection

.. container:: cell markdown

   You can find a list of all fibers used in the measurement in the
   Fibers table. The Fibers table and its associated columns can be
   accessed similar to the Spectra table by searching for a match in the
   the detectid column.

.. container:: cell code

   .. code:: python

      fibers = detects.hdfile.root.Fibers
      fibers.cols

   .. container:: output execute_result

      ::

         /Fibers.cols (Cols), 20 columns
           detectid (Column(9397618,), int64)
           ra (Column(9397618,), float32)
           dec (Column(9397618,), float32)
           multiframe (Column(9397618,), |S20)
           fibnum (Column(9397618,), int32)
           x_ifu (Column(9397618,), float32)
           y_ifu (Column(9397618,), float32)
           date (Column(9397618,), int32)
           obsid (Column(9397618,), int32)
           expnum (Column(9397618,), int32)
           distance (Column(9397618,), float32)
           timestamp (Column(9397618,), |S17)
           wavein (Column(9397618,), float32)
           flag (Column(9397618,), int32)
           weight (Column(9397618,), float32)
           ADC (Column(9397618, 5), ('<f4', (5,)))
           amp (Column(9397618,), |S2)
           ifuid (Column(9397618,), |S3)
           ifuslot (Column(9397618,), |S3)
           specid (Column(9397618,), |S3)

.. container:: cell markdown

   Access the fiber table for the above source:

.. container:: cell code

   .. code:: python

      fiber_table = fibers.read_where("detectid == detectid_nice_lae") 

.. container:: cell code

   .. code:: python

      Table(fiber_table)

   .. container:: output execute_result

      ::

         <Table length=16>
          detectid      ra       dec         multiframe      ... ifuid  ifuslot specid
           int64     float32   float32        bytes20        ... bytes3  bytes3 bytes3
         ---------- --------- --------- -------------------- ... ------ ------- ------
         1000208773 199.35771 51.066715 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35664 51.066483 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35793 51.067413 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35686  51.06718 multi_025_076_032_RU ...    032     076    025
         1000208773  199.3558 51.066948 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35707 51.067875 multi_025_076_032_RU ...    032     076    025
         1000208773   199.356 51.067642 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35814 51.067024 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35707  51.06679 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35728 51.067486 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35622 51.067257 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35728  51.06641 multi_025_076_032_RU ...    032     076    025
         1000208773  199.3575 51.067104 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35643  51.06687 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35771   51.0678 multi_025_076_032_RU ...    032     076    025
         1000208773 199.35664 51.067566 multi_025_076_032_RU ...    032     076    025

.. container:: cell markdown

   When you are done with the HDF5 file, close it. The data that you
   extracted into tables and arrays will remain.

.. container:: cell code

   .. code:: python

      detects.hdfile.close()

.. container:: cell markdown

   .. rubric:: Accessing the ELiXer Classifications
      :name: accessing-the-elixer-classifications

.. container:: cell code

   .. code:: python

      file_elix = tb.open_file(config.elixerh5)

.. container:: cell code

   .. code:: python

      file_elix.root.Classifications

   .. container:: output execute_result

      ::

         /Classifications (Table(690868,)) ''
           description := {
           "detectid": Int64Col(shape=(), dflt=0, pos=0),
           "plae_poii_hetdex": Float32Col(shape=(), dflt=0.0, pos=1),
           "aperture_mag": Float32Col(shape=(), dflt=0.0, pos=2),
           "plae_poii_aperture": Float32Col(shape=(), dflt=0.0, pos=3),
           "aperture_filter": StringCol(itemsize=15, shape=(), dflt=b'', pos=4),
           "mag_match": Float32Col(shape=(), dflt=0.0, pos=5),
           "cat_filter": StringCol(itemsize=15, shape=(), dflt=b'', pos=6),
           "plae_poii_cat": Float32Col(shape=(), dflt=0.0, pos=7),
           "dec": Float32Col(shape=(), dflt=0.0, pos=8),
           "dec_match": Float32Col(shape=(), dflt=0.0, pos=9),
           "dist_match": Float32Col(shape=(), dflt=0.0, pos=10),
           "ra": Float32Col(shape=(), dflt=0.0, pos=11),
           "ra_match": Float32Col(shape=(), dflt=0.0, pos=12),
           "z_prelim": Float32Col(shape=(), dflt=0.0, pos=13)}
           byteorder := 'little'
           chunkshape := (799,)
           autoindex := True
           colindexes := {
             "detectid": Index(9, full, shuffle, zlib(1)).is_csi=True}

.. container:: cell markdown

   Note: these are also appended to the Detections() class object. Each
   column in the above table can be accessed as an attribute of the
   Detections() class object. For example, the probability of LAE to OII
   measured from the HETDEX continuum is:

.. container:: cell code

   .. code:: python

      detects.plae_poii_hetdex

   .. container:: output execute_result

      ::

         array([9.9900000e+02, 7.2290176e-01, 9.9900000e+02, ..., 1.4184024e+00,
                8.7493747e-01, 2.0698796e-01], dtype=float32)

.. container:: cell markdown

   or the nearest neighbour magnitude in an ancillary photometric
   catalog is:

.. container:: cell code

   .. code:: python

      detects.mag_match

   .. container:: output execute_result

      ::

         array([23.479797, 99.9     , 99.9     , ..., 99.9     , 99.9     ,
                99.9     ], dtype=float32)

.. container:: cell markdown

   and this comes from the filter:

.. container:: cell code

   .. code:: python

      detects.cat_filter

   .. container:: output execute_result

      ::

         array(['r', '?', '?', ..., '?', '?', '?'], dtype='<U15')

.. |image0| image:: images/68ce842e38ec438396778621fb6847f390b730c5.png
.. |image1| image:: images/4b0f69464bddc450c2dadb63ed71e67d5b1ba89a.png
.. |image2| image:: images/19c8b4ebadce2360bfdeea697d9d1cd601c65198.png
