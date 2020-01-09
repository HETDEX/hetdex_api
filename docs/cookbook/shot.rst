.. container:: cell markdown

   .. rubric:: 02 - Accessing Processed Data: The Fibers and Images
      Classes
      :name: 02---accessing-processed-data-the-fibers-and-images-classes

.. container:: cell markdown

   All data for a single HETDEX observation, a dither set of three
   exposures, is stored in one HDF5 file. This file contains the
   processed data in the form of 2D image data and individual fiber
   data. It also contains additional branches including information
   about the astrometric and flux-calibration details. Most users will
   be most interested in the flux-calibrated fiber information and 2D
   image data. A user can access all fibers observated at a specific
   RA/DEC and perform analysis on these retrieved fibers. This data can
   be retrieved through the Fibers class API described here. The 2D
   image data represents the intermediary processed data. It provides
   images for each IFU/amp (essential the original CCD image of the
   fiber spectra) along each step of calibration. Please refer to
   Panacea documentation for detailed information on the reduction
   process.

   We begin by first introducing some global shot functions and then
   delve into the Fibers Class API.

.. container:: cell markdown

   .. rubric:: Opening a shot
      :name: opening-a-shot

.. container:: cell markdown

   If you want to just open the shot HDF5 file and explore contents, you
   must know the shot ID: either its integer shotid value (20180124010)
   or the datevobs ('20180124v010'). The last notebook showed how to
   find a list of shots based on querying the Survey Class. A shot HDF5
   file can then be opened using the ``open_shot_file()`` global
   command:

.. container:: cell code

   .. code:: python

      %matplotlib inline
      import matplotlib.pyplot as plt
      from matplotlib.colors import LogNorm
      from astropy.visualization import ZScaleInterval

      import astropy.units as u
      from astropy.coordinates import SkyCoord

      from hetdex_api.shot import *

.. container:: cell code

   .. code:: python

      fileh = open_shot_file('20180124v010')

.. container:: cell markdown

   To see all contents of the shot file, you can do ``print(fileh)``:

.. container:: cell markdown

   .. rubric:: Initiate the Fibers Class
      :name: initiate-the-fibers-class

.. container:: cell markdown

   To access the fibers of a given shot, you can open the HDF5 shot and
   load in the RA/DEC of each fiber as an astropy coordinates object in
   one instance. This will take from 2-10 seconds as it loads in the
   full Fibers dataset.

.. container:: cell code

   .. code:: python

      fibers = Fibers('20180124v010')

.. container:: cell markdown

   This will also initialize an array of rectified wavelengths to be
   used when plotting the 'calfib' column representing the calibrated
   flux for each individual fiber.

.. container:: cell code

   .. code:: python

      fibers.wave_rect

   .. container:: output execute_result

      ::

         array([3470., 3472., 3474., ..., 5536., 5538., 5540.])

.. container:: cell markdown

   The following are functions that act upon the fibers class to aid in
   querying and retrieving fibers. Use ``query_region_idx()`` to get the
   index of all fibers within a defined circular aperture. The center of
   the aperture must be an astropy coords object, for example:

.. container:: cell code

   .. code:: python

      coords = SkyCoord(150.025513 * u.deg, 2.087767 * u.deg, frame='icrs')

.. container:: cell code

   .. code:: python

      idx = fibers.query_region_idx(coords, radius=3./3600.)
      idx

   .. container:: output execute_result

      ::

         array([], dtype=int64)

.. container:: cell raw

   .. raw:: ipynb

      To plot the fiber spectra for each fiber, use `plot_fiber_spectrum()`:

.. container:: cell markdown

   We got plot up all spectra using the plot_fiber_spectrum which takes
   a row index value as an argument and acts upon the Fibers class
   object.

.. container:: cell code

   .. code:: python

      plt.figure(figsize=(8,6))
      for i in idx :
          fibers.plot_fiber_spectrum(i)

   .. container:: output display_data

      ::

         <matplotlib.figure.Figure at 0x2ac4ecfc09e8>

.. container:: cell markdown

   Using the xlim and ylim options, we can vary the axes range:

.. container:: cell code

   .. code:: python

      plt.figure(figsize=(8,6))
      for i in idx:
          fibers.plot_fiber_spectrum(i, xlim=[3680,3740])

   .. container:: output display_data

      ::

         <matplotlib.figure.Figure at 0x2ac4ecfcc630>

.. container:: cell markdown

   Each fiber can be saved to a text file as follows:

.. container:: cell code

   .. code:: python

      for i in idx:
          fibers.save_fiber_spectrum(i, file='spec_' + str(i) + '.dat')

.. container:: cell markdown

   .. rubric:: Some other Fibers class functions
      :name: some-other-fibers-class-functions

.. container:: cell markdown

   To find the closet fiber to a set of coordinates:

.. container:: cell code

   .. code:: python

      idx = fibers.get_closest_fiber(coords)

.. container:: cell markdown

   To find the x,y image value in the 2D images arrays, use get_image_xy
   on the fibers class. A user must provide both a fiber index and a
   wavelength:

.. container:: cell code

   .. code:: python

      x, y = fibers.get_image_xy(idx, 3710)

.. container:: cell markdown

   .. rubric:: Get Image cutouts:
      :name: get-image-cutouts

.. container:: cell markdown

   An image cutout can be extracted for a specific shot either around a
   set of coordinates or for a specific multiframe ID (this is the
   IFU/amp ID). For example:

.. container:: cell code

   .. code:: python

      implot = get_image2D_cutout('20180124v010', coords, 3710)

.. container:: cell code

   .. code:: python

      zscale = ZScaleInterval(contrast=0.5,krej=2.5) 
      vmin, vmax = zscale.get_limits(values=implot)
      plt.imshow(implot,vmin=vmin, vmax=vmax, origin="lower",cmap=plt.get_cmap('gray'),interpolation="none")

   .. container:: output execute_result

      ::

         <matplotlib.image.AxesImage at 0x2ac544703550>

   .. container:: output display_data

      |image0|

.. container:: cell markdown

   Or we can grab an entire amp of interest for a specific shot:

.. container:: cell code

   .. code:: python

      multiframe_obj = 'multi_319_083_023_RL'

.. container:: cell code

   .. code:: python

      im_amp = get_image2D_amp('20180124v010', multiframe_obj)
      zscale = ZScaleInterval(contrast=0.5,krej=2.5) 
      vmin, vmax = zscale.get_limits(values=im_amp)
      plt.imshow(im_amp,vmin=vmin, vmax=vmax, origin="lower",cmap=plt.get_cmap('gray'),interpolation="none")

   .. container:: output execute_result

      ::

         <matplotlib.image.AxesImage at 0x2ac5457d4630>

   .. container:: output display_data

      |image1|

.. |image0| image:: rst/cc6f620da3ea38d3c148e7305d43b6844f8623fd.png
.. |image1| image:: rst/199682cdd638c8603a0b47c78ec4ca652c10f911.png
