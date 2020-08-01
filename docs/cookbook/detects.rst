Detections Database and API
===========================

.. container:: cell markdown

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
      import pickle
      import subprocess
      import numpy as np
      import tables as tb
      import matplotlib.pyplot as plt

      from astropy.io import ascii
      from astropy.table import Table, Column
      from astropy.coordinates import SkyCoord
      import astropy.units as u

      from hetdex_api.config import HDRconfig
      from hetdex_api.detections import Detections
      from hetdex_api.elixer_widget_cls import ElixerWidget

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

      #old way still works
      #detects = Detections('hdr1').refine()

      # but this is the new way
      detects = Detections(survey='hdr2.1', catalog_type='lines')

      # or if you want to open the continuum source catalog:
      # detects = Detections(survey='hdr2.1', catalog_type='continuum')

.. container:: cell markdown

   .. rubric:: Note if you do not want to load the whole table, but just
      access spectra for a specific detectid:
      :name: note-if-you-do-not-want-to-load-the-whole-table-but-just-access-spectra-for-a-specific-detectid

.. container:: cell code

   .. code:: python

      det_object = Detections('hdr2.1', loadtable=False)

.. container:: cell code

   .. code:: python

      spec = det_object.get_spectrum(2100191119)

.. container:: cell code

   .. code:: python

      spec

   .. container:: output execute_result

      ::

         <Table length=1036>
          wave1d             spec1d                     spec1d_err         
         Angstrom 1e-17 erg / (Angstrom cm2 s) 1e-17 erg / (Angstrom cm2 s)
         float32            float32                      float32           
         -------- ---------------------------- ----------------------------
           3470.0                 -0.016949153                         9.75
           3472.0                 -0.016949153                         9.75
           3474.0                 -0.016949153                         9.75
           3476.0                 -0.016949153                         9.75
           3478.0                 -0.016949153                         9.75
           3480.0                 -0.016949153                         9.75
           3482.0                 -0.016949153                         9.75
           3484.0                        0.875                     9.254767
           3486.0                    0.5010593                     9.188029
              ...                          ...                          ...
           5522.0                  -0.05632306                     1.176408
           5524.0                  -0.05632306                     1.176408
           5526.0                  -0.05632306                     1.176408
           5528.0                  -0.05632306                     1.176408
           5530.0                  -0.05632306                     1.176408
           5532.0                  -0.05632306                     1.176408
           5534.0                  -0.05632306                     1.176408
           5536.0                  -0.05632306                     1.176408
           5538.0                  -0.05632306                     1.176408
           5540.0                  -0.05632306                     1.176408

.. container:: cell markdown

   Here are a list of attributes built into the Detections class:

.. container:: cell code

   .. code:: python

      detects.__dict__.keys()

   .. container:: output execute_result

      ::

         dict_keys(['survey', 'filename', 'hdfile', 'loadtable', 'detectid', 'shotid', 'ra', 'dec', 'date', 'obsid', 'wave', 'wave_err', 'flux', 'flux_err', 'linewidth', 'linewidth_err', 'continuum', 'continuum_err', 'sn', 'sn_err', 'chi2', 'chi2_err', 'multiframe', 'fibnum', 'x_raw', 'y_raw', 'amp', 'chi2fib', 'detectname', 'expnum', 'fiber_id', 'ifuid', 'ifuslot', 'inputid', 'noise_ratio', 'specid', 'weight', 'x_ifu', 'y_ifu', 'combined_continuum', 'combined_continuum_err', 'combined_plae', 'combined_plae_err', 'mag_sdss_g', 'mag_sdss_g_err', 'plae_classification', 'plae_sdss_g', 'plae_sdss_g_max', 'plae_sdss_g_min', 'gmag', 'gmag_err', 'field', 'fwhm', 'fluxlimit_4540', 'throughput', 'n_ifu', 'vis_class', 'coords'])

.. container:: cell markdown

   If you prefer working in astropy tables, you can grab it this way:

.. container:: cell code

   .. code:: python

      detect_table = detects.return_astropy_table()

.. container:: cell code

   .. code:: python

      detect_table

   .. container:: output execute_result

      ::

         <Table length=2591424>
          detectid         fwhm        ... plae_sdss_g_max plae_sdss_g_min
           int64         float64       ...     float32         float32    
         ---------- ------------------ ... --------------- ---------------
         2100000000 2.3224666118621826 ...          1000.0          1000.0
         2100000001 2.3224666118621826 ...          1000.0          1000.0
         2100000002 2.3224666118621826 ...          1000.0          1000.0
         2100000003 2.3224666118621826 ...          1000.0          1000.0
         2100000004 2.3224666118621826 ...          1000.0          1000.0
         2100000005 2.3224666118621826 ...           0.001           0.001
         2100000006 2.3224666118621826 ...           0.001           0.001
         2100000007 2.3224666118621826 ...     0.013710628           0.001
         2100000008 2.3224666118621826 ...      0.03905671    0.0040173857
         2100000009 2.3224666118621826 ...          1000.0          1000.0
                ...                ... ...             ...             ...
         2102591414 1.2000000476837158 ...          1000.0          1000.0
         2102591415 1.2000000476837158 ...          1000.0          1000.0
         2102591416 1.2000000476837158 ...          1000.0          1000.0
         2102591417 1.2000000476837158 ...          1000.0          1000.0
         2102591418 1.2000000476837158 ...          1000.0          1000.0
         2102591419 1.2000000476837158 ...          1000.0          1000.0
         2102591420 1.2000000476837158 ...          1000.0          1000.0
         2102591421 1.2000000476837158 ...          1000.0          1000.0
         2102591422 1.2000000476837158 ...        715.8892        56.65825
         2102591423 1.2000000476837158 ...      0.05794077     0.043867636

.. container:: cell markdown

   .. rubric:: How we made the subset catalog for the team:
      :name: how-we-made-the-subset-catalog-for-the-team

.. container:: cell code

   .. code:: python

      # remove obsolete detections 
      config = HDRconfig('hdr2.1')
      good_det = np.array(pickle.load( open( config.baddetectmask, "rb")), dtype=bool)

.. container:: cell code

   .. code:: python

      sel_chi2 = detects.chi2 <= 1.6

.. container:: cell code

   .. code:: python

      sel_field = (detects.field == b'cosmos') | (detects.field == b'dex-fall') | (detects.field == b'dex-spring') | (detects.field == b'egs') | (detects.field == b'goods-n')

.. container:: cell code

   .. code:: python

      sel_cat = good_det*sel_chi2*sel_field

.. container:: cell code

   .. code:: python

      team_table = detect_table[sel_cat]

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
             [(149.79932 , 1.986114), (149.80261 , 1.991804),
              (149.80466 , 1.994646), ..., ( 36.49977 , 0.405466),
              ( 36.496384, 0.411001), ( 36.49315 , 0.408433)]>

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

         6

.. container:: cell markdown

   .. rubric:: Find a direct line match
      :name: find-a-direct-line-match

.. container:: cell markdown

   If you want to find an exact line match you can use the function
   ``find_match()``

.. container:: cell code

   .. code:: python

      obj_coords = SkyCoord(199.35704 * u.deg, 51.06718 * u.deg, frame='icrs')

.. container:: cell code

   .. code:: python

      wave_obj = 3836.

.. container:: cell code

   .. code:: python

      idx = detects.find_match(obj_coords, wave=wave_obj, radius=5.*u.arcsec, dwave=5 )

.. container:: cell code

   .. code:: python

      detects.detectid[idx]

   .. container:: output execute_result

      ::

         array([2100191119])

.. container:: cell code

   .. code:: python

      detect_table[idx]

   .. container:: output execute_result

      ::

         <Table length=1>
          detectid         fwhm        ... plae_sdss_g_max plae_sdss_g_min
           int64         float64       ...     float32         float32    
         ---------- ------------------ ... --------------- ---------------
         2100191119 1.4780957698822021 ...          1000.0          1000.0

.. container:: cell markdown

   .. rubric:: Check out matched sources in the ElixerWidget
      :name: check-out-matched-sources-in-the-elixerwidget

.. container:: cell markdown

   For this example, we have found 12 detections in this region, we can
   examine these via the ELiXer reports using the ``ElixerWidget()``
   class from ``hetdex_api.elixer_widget_cls.py``. To do so we need to
   save the detectid list to examine in the widget.

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

      elix_widget = ElixerWidget(detectlist = detects_in_region.detectid)
      #elix_widget = ElixerWidget(detectfile='detects_obj.txt')

   .. container:: output display_data

      .. code:: json

         {"model_id":"4ae7cba75e744f61a7ee036b223a8cce","version_major":2,"version_minor":0}

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

         /data/05350/ecooper/hdr2.1/detect/detect_hdr2.1.h5 (File) 'HDR2.1 Detections Database'
         Last modif.: 'Fri Jul 31 12:19:14 2020'
         Object Tree: 
         / (RootGroup) 'HDR2.1 Detections Database'
         /Detections (Table(2591424,)) 'HETDEX Line Detection Catalog'
         /Elixer (Table(2591424,)) 'Elixer Info'
         /Fibers (Table(49570585,)) 'Fiber info for each detection'
         /Spectra (Table(2591424,)) '1D Spectra for each Line Detection'

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

         /Spectra.cols (Cols), 12 columns
           detectid (Column(2591424,), int64)
           wave1d (Column(2591424, 1036), ('<f4', (1036,)))
           spec1d (Column(2591424, 1036), ('<f4', (1036,)))
           spec1d_err (Column(2591424, 1036), ('<f4', (1036,)))
           counts1d (Column(2591424, 1036), ('<f4', (1036,)))
           counts1d_err (Column(2591424, 1036), ('<f4', (1036,)))
           apsum_counts (Column(2591424, 1036), ('<f4', (1036,)))
           apsum_counts_err (Column(2591424, 1036), ('<f4', (1036,)))
           apcor (Column(2591424, 1036), ('<f4', (1036,)))
           flag_pix (Column(2591424, 1036), ('<f4', (1036,)))
           spec1d_nc (Column(2591424, 1036), ('<f4', (1036,)))
           spec1d_nc_err (Column(2591424, 1036), ('<f4', (1036,)))

.. container:: cell markdown

   Flux calibrated, psf-weighted 1D spectra can be retrieved via the API
   for a single detectid through the function ``get_spectrum``:

.. container:: cell code

   .. code:: python

      detectid_nice_lae = 2000202849
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

         /Fibers.cols (Cols), 23 columns
           detectid (Column(15019537,), int64)
           ra (Column(15019537,), float32)
           dec (Column(15019537,), float32)
           multiframe (Column(15019537,), |S20)
           fiber_id (Column(15019537,), |S38)
           x_ifu (Column(15019537,), float32)
           y_ifu (Column(15019537,), float32)
           date (Column(15019537,), int32)
           obsid (Column(15019537,), int32)
           expnum (Column(15019537,), int32)
           distance (Column(15019537,), float32)
           timestamp (Column(15019537,), |S17)
           wavein (Column(15019537,), float32)
           flag (Column(15019537,), int32)
           weight (Column(15019537,), float32)
           ADC (Column(15019537, 5), ('<f4', (5,)))
           amp (Column(15019537,), |S2)
           fibnum (Column(15019537,), int32)
           ifuid (Column(15019537,), |S3)
           ifuslot (Column(15019537,), |S3)
           specid (Column(15019537,), |S3)
           x_raw (Column(15019537,), int32)
           y_raw (Column(15019537,), int32)

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

         <Table length=15>
          detectid      ra       dec         multiframe      ... specid x_raw y_raw
           int64     float32   float32        bytes20        ... bytes3 int32 int32
         ---------- --------- --------- -------------------- ... ------ ----- -----
         2000202849 199.35779 51.066734 multi_025_076_032_RU ...    025   180   123
         2000202849 199.35672   51.0665 multi_025_076_032_RU ...    025   180   132
         2000202849   199.358  51.06743 multi_025_076_032_RU ...    025   178   291
         2000202849 199.35693   51.0672 multi_025_076_032_RU ...    025   178   300
         2000202849 199.35715  51.06789 multi_025_076_032_RU ...    025   176   476
         2000202849 199.35608  51.06766 multi_025_076_032_RU ...    025   176   484
         2000202849  199.3582 51.067028 multi_025_076_032_RU ...    025   181   123
         2000202849 199.35713 51.066795 multi_025_076_032_RU ...    025   180   132
         2000202849 199.35735 51.067493 multi_025_076_032_RU ...    025   178   300
         2000202849 199.35628  51.06726 multi_025_076_032_RU ...    025   178   308
         2000202849 199.35733  51.06643 multi_025_076_032_RU ...    025   180   132
         2000202849 199.35754 51.067123 multi_025_076_032_RU ...    025   178   300
         2000202849 199.35649 51.066895 multi_025_076_032_RU ...    025   178   308
         2000202849 199.35776 51.067818 multi_025_076_032_RU ...    025   176   476
         2000202849 199.35669 51.067585 multi_025_076_032_RU ...    025   176   485

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

      config = HDRconfig(survey='hdr2')
      file_elix = tb.open_file(config.elixerh5)

.. container:: cell code

   .. code:: python

      file_elix.root.Detections

   .. container:: output execute_result

      ::

         /Detections (Table(1098592,)) 'ELiXer Detection Summary Table'
           description := {
           "detectid": Int64Col(shape=(), dflt=0, pos=0),
           "detectname": StringCol(itemsize=64, shape=(), dflt=b'', pos=1),
           "elixer_version": StringCol(itemsize=16, shape=(), dflt=b'', pos=2),
           "elixer_datetime": StringCol(itemsize=21, shape=(), dflt=b'', pos=3),
           "ra": Float32Col(shape=(), dflt=-999.999, pos=4),
           "dec": Float32Col(shape=(), dflt=-999.999, pos=5),
           "wavelength_obs": Float32Col(shape=(), dflt=-999.999, pos=6),
           "wavelength_obs_err": Float32Col(shape=(), dflt=-999.999, pos=7),
           "chi2": Float32Col(shape=(), dflt=-999.999, pos=8),
           "chi2_err": Float32Col(shape=(), dflt=-999.999, pos=9),
           "combined_continuum": Float32Col(shape=(), dflt=-999.999, pos=10),
           "combined_continuum_err": Float32Col(shape=(), dflt=-999.999, pos=11),
           "combined_plae": Float32Col(shape=(), dflt=-999.999, pos=12),
           "combined_plae_err": Float32Col(shape=(), dflt=-999.999, pos=13),
           "continuum_full_spec": Float32Col(shape=(), dflt=-999.999, pos=14),
           "continuum_full_spec_err": Float32Col(shape=(), dflt=-999.999, pos=15),
           "continuum_line": Float32Col(shape=(), dflt=-999.999, pos=16),
           "continuum_line_err": Float32Col(shape=(), dflt=-999.999, pos=17),
           "continuum_sdss_g": Float32Col(shape=(), dflt=-999.999, pos=18),
           "continuum_sdss_g_err": Float32Col(shape=(), dflt=-999.999, pos=19),
           "eqw_rest_lya_full_spec": Float32Col(shape=(), dflt=-999.999, pos=20),
           "eqw_rest_lya_full_spec_err": Float32Col(shape=(), dflt=-999.999, pos=21),
           "eqw_rest_lya_line": Float32Col(shape=(), dflt=-999.999, pos=22),
           "eqw_rest_lya_line_err": Float32Col(shape=(), dflt=-999.999, pos=23),
           "eqw_rest_lya_sdss_g": Float32Col(shape=(), dflt=-999.999, pos=24),
           "eqw_rest_lya_sdss_g_err": Float32Col(shape=(), dflt=-999.999, pos=25),
           "fieldname": StringCol(itemsize=32, shape=(), dflt=b'', pos=26),
           "flux_line": Float32Col(shape=(), dflt=-999.999, pos=27),
           "flux_line_err": Float32Col(shape=(), dflt=-999.999, pos=28),
           "fwhm_line_aa": Float32Col(shape=(), dflt=-999.999, pos=29),
           "fwhm_line_aa_err": Float32Col(shape=(), dflt=-999.999, pos=30),
           "ifuid": StringCol(itemsize=3, shape=(), dflt=b'', pos=31),
           "ifuslot": StringCol(itemsize=3, shape=(), dflt=b'', pos=32),
           "mag_full_spec": Float32Col(shape=(), dflt=-999.999, pos=33),
           "mag_full_spec_err": Float32Col(shape=(), dflt=-999.999, pos=34),
           "mag_sdss_g": Float32Col(shape=(), dflt=-999.999, pos=35),
           "mag_sdss_g_err": Float32Col(shape=(), dflt=-999.999, pos=36),
           "multiline_flag": BoolCol(shape=(), dflt=False, pos=37),
           "multiline_frac_score": Float32Col(shape=(), dflt=-999.999, pos=38),
           "multiline_name": StringCol(itemsize=16, shape=(), dflt=b'', pos=39),
           "multiline_prob": Float32Col(shape=(), dflt=-999.999, pos=40),
           "multiline_raw_score": Float32Col(shape=(), dflt=-999.999, pos=41),
           "multiline_rest_w": Float32Col(shape=(), dflt=-999.999, pos=42),
           "multiline_z": Float32Col(shape=(), dflt=-999.999, pos=43),
           "obsid": Int32Col(shape=(), dflt=0, pos=44),
           "plae_classification": Float32Col(shape=(), dflt=-999.999, pos=45),
           "plae_full_spec": Float32Col(shape=(), dflt=-999.999, pos=46),
           "plae_full_spec_max": Float32Col(shape=(), dflt=-999.999, pos=47),
           "plae_full_spec_min": Float32Col(shape=(), dflt=-999.999, pos=48),
           "plae_line": Float32Col(shape=(), dflt=-999.999, pos=49),
           "plae_line_max": Float32Col(shape=(), dflt=-999.999, pos=50),
           "plae_line_min": Float32Col(shape=(), dflt=-999.999, pos=51),
           "plae_sdss_g": Float32Col(shape=(), dflt=-999.999, pos=52),
           "plae_sdss_g_max": Float32Col(shape=(), dflt=-999.999, pos=53),
           "plae_sdss_g_min": Float32Col(shape=(), dflt=-999.999, pos=54),
           "pseudo_color_blue_flux": Float32Col(shape=(), dflt=-999.999, pos=55),
           "pseudo_color_blue_flux_err": Float32Col(shape=(), dflt=-999.999, pos=56),
           "pseudo_color_flag": Int64Col(shape=(), dflt=0, pos=57),
           "pseudo_color_red_flux": Float32Col(shape=(), dflt=-999.999, pos=58),
           "pseudo_color_red_flux_err": Float32Col(shape=(), dflt=-999.999, pos=59),
           "pseudo_color_rvb_ratio": Float32Col(shape=(), dflt=-999.999, pos=60),
           "pseudo_color_rvb_ratio_err": Float32Col(shape=(), dflt=-999.999, pos=61),
           "response": Float32Col(shape=(), dflt=-999.999, pos=62),
           "seeing_gaussian": Float32Col(shape=(), dflt=-999.999, pos=63),
           "seeing_moffat": Float32Col(shape=(), dflt=-999.999, pos=64),
           "shotid": Int64Col(shape=(), dflt=0, pos=65),
           "sn": Float32Col(shape=(), dflt=-999.999, pos=66),
           "sn_err": Float32Col(shape=(), dflt=-999.999, pos=67),
           "specid": StringCol(itemsize=3, shape=(), dflt=b'', pos=68)}
           byteorder := 'little'
           chunkshape := (159,)
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

      #detects.plae_poii_hetdex

.. container:: cell markdown

   or the nearest neighbour magnitude in an ancillary photometric
   catalog is:

.. container:: cell code

   .. code:: python

      #detects.mag_match

.. container:: cell markdown

   and this comes from the filter:

.. container:: cell code

   .. code:: python

      #detects.cat_filter

.. container:: cell code

   .. code:: python

.. |image0| image:: images/b51469a4acf94ad2b572f96099b221fc147462c1.png
.. |image1| image:: images/4c3a0599c87e4ceb8a4a579e107bdf4f27181b6c.png
.. |image2| image:: images/d11ca64661315436d354fb4827c57046da40b2f5.png
