``hetdex_api`` release notes
==========================

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
