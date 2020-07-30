``hetdex_api`` release notes
==========================

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
