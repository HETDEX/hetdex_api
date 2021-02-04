``hetdex_api`` release notes
==========================

Release 0.8.7
-------------
- updates for 2.1.3 catalogs
- including apcor as column in catalogs and detections.py
- finalized source_catalog.py command line functionality
- incorporated continuum h5 file into amp_widget.py and elixer_widget_cls.py

Release 0.8.6
-------------
- Updates to SourceCatalog.ipynb

Release 0.8.5
-------------
- Updates to QueryWidget

Release 0.8.4
-------------
- Updates to hetdex_tools/source_catalog.py

Release 0.8.3
-------------

- Release for 2.1.2 catalogs
- updates to make_curated_catalog
- created hetdex_tools/source_catalog.py to make a unique source catalog
- hetdex_tools/interpolate.py developed for data cube and collapsed image
- developed hetdex_api/cube_widget.py to scan a HETDEX cube
- developed hetdex_api/amp_widget.py to interactively explore individual amp observations
- masking tools: hetdex_tools/galmask.py, hetdex_api/mask.py
- organized and made lots of new notebooks

Release 0.7
-----------

- added hetdex_tools/get_contour.py to perform grid line fitting

- updated hetdex_api/elixer_widget_cls.py for new team classifying

- add API for mask development hetdex_api/mask.py

- updated line fitting tools in hetdex_tools/line_fitting.py

- added hetdex_tools/make_curated_catalog.py to document catalog creation

- added tool to grab flim slice and detections hetdex_tools/plot_flim_slice.py

Release 0.6
-----------

- Still release for HDR2.1.1

- using healpy to speed up query speed on FiberIndex

- small changes in docs for cleaner organization

- catalog curation: added 82 new bad shots, a bad fiber and some bad amps, caught 960 charge trap issues.

Release 0.5
-----------

- Release for HDR2.1.1 

- Detections class function refine() now operational for HDR2.1

- Added curated_version option for Detections class to pull stable catalog

- small changes in docs for cleaner organization

- catalog curation: added 82 new bad shots, a bad fiber and some bad amps, caught 960 charge trap issues. 

Development version @ trunk
---------------------------

- The ``f50vals`` array in the ``SensitivityCube`` class has been
  replaced by an array of the 1 sigma noise called ``sigmas``

- When initialising a ``SensitivityCube`` an nsigma parameter
  can be passed that converts the input to signmas via input/nsigma

- You must now pass a ``sncut`` parameter to all the flux limit
  functions, which specifies the S/N cut

- Different flux limit models for different releases are
  now available via the ``flim_model`` option

- The alphas array can be 3D, in which case an alpha can be specified for
each location
