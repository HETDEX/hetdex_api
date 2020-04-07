
Extracting Spectra
==================

This notebook demonstrates how to grab 1D aperture summed HETDEX spectra
for an input of ID, RA and DEC using the ``Extract`` Class API from
``hetdex_api``. This can be done interactively using ``get_spectra``
from hte ``hetdex_tools.get_spec`` module. It can also be done in the
command line using the quick entry call ``hetdex_get_spec`` providing
you have hetdex\_api pip installed.

Examples of what you might like to do with the spectra afterwards is
shown later. The output is stored in an astropy table of spectra. For
every HETDEX observation where spectra is found, a spectra is given. It
is up to the user to combine the spectra afterwards.

IMPORTANT NOTE OF CAUTION WITH RUNNING ON TACC!!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because this script involves opening the Fibers class object which
contains all fiber spectra from a 3 dither observation, you will be
pulling in a lot of memory for each shot that is open. **NEVER** run
this script from a login node on TACC. A login node is a node you access
when you ssh in.

You need to request a compute node instead by either

(1) using the idev command : ``idev -t 04:00:00``

(2) using a jupyter notebook

(3) or by submitting the job into the slurm job scheduler (this is
    probably only needed if you are doing more than 1000 extractions per
    shot on more than 1000 shots)

Import all necessary python packages.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are mainly for working within the notebook. The command line tool
already has the necessary preamble built in.

.. code:: ipython3

    %matplotlib inline
    import matplotlib.pyplot as plt
    import numpy as np
    
    import astropy.units as u
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    from astropy.table import Table, join
    
    from hetdex_tools.get_spec import get_spectra

Getting all spectra at a specified RA/DEC. This will search through all shots in HDR2
-------------------------------------------------------------------------------------

If a shotid is not specified the program will search for any shot within
HDR2 that overlaps within an 11 arcmin radius of the input coordinates.
Because of the non-contiguous VIRUS footprint, there is no guarantee the
aperture defined by the input ra/dec/rad will contain enough fibers to
do a measurement. The aperture radius is 3" by default or can be
specified with the --rad argument.

Open a catalog of IDs, RAs, DECs:

.. code:: ipython3

    input_cat = Table.read('/work/05350/ecooper/stampede2/hdr2-tests/hps-muse/highz_emitters.fits')

.. code:: ipython3

    input_cat[0:5]




.. raw:: html

    <i>Table length=5</i>
    <table id="table47889346717848" class="table-striped table-bordered table-condensed">
    <thead><tr><th>ID</th><th>RA</th><th>DEC</th><th>WAVE</th><th>FLUX</th><th>FLUXE_L</th><th>FLUXE_H</th><th>z</th></tr></thead>
    <thead><tr><th></th><th>deg</th><th>deg</th><th>Angstrom</th><th>1e-17 erg / (cm2 s)</th><th>1e-17 erg / (cm2 s)</th><th>1e-17 erg / (cm2 s)</th><th></th></tr></thead>
    <thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>3</td><td>35.30979166666667</td><td>-4.527130555555556</td><td>4973.93</td><td>19.9</td><td>3.1</td><td>4.7</td><td>3.0915</td></tr>
    <tr><td>4</td><td>35.31191666666666</td><td>-4.532388888888889</td><td>5261.37</td><td>42.6</td><td>12.4</td><td>11.2</td><td>1.7561</td></tr>
    <tr><td>5</td><td>35.31308333333333</td><td>-4.531666666666666</td><td>4270.67</td><td>342.1</td><td>14.3</td><td>16.5</td><td>1.757</td></tr>
    <tr><td>6</td><td>35.31816666666666</td><td>-4.4926</td><td>4591.58</td><td>32.7</td><td>3.6</td><td>3.5</td><td>2.777</td></tr>
    <tr><td>11</td><td>35.32691666666666</td><td>-4.459319444444445</td><td>4590.82</td><td>21.2</td><td>4.6</td><td>4.7</td><td>2.7764</td></tr>
    </table>



``get_spectra()`` requires an astropy coordinates object list as an
input.

.. code:: ipython3

    input_coords = SkyCoord(ra=input_cat['RA'], dec=input_cat['DEC'])

.. code:: ipython3

    sources = get_spectra(input_coords, ID=input_cat['ID'])

get\_spectra() options
----------------------

There are a few options to consider when running get\_spectra():

.. code:: ipython3

    help(get_spectra)


.. parsed-literal::

    Help on function get_spectra in module hetdex_tools.get_spec:
    
    get_spectra(coords, ID=None, rad=3.0, multiprocess=True, shotid=None, survey='hdr2', tpmin=0.09, ffsky=False)
        Function to retrieve PSF-weighted, ADR and aperture corrected
        spectral extractions of HETDEX fibers. It will search all shots
        within a specific HETDEX Data Release and return a table of
        spectra for each extraction per shot in which more than 7 fibers
        are found in order to generate an extracted spectrum.
        
        Parameters
        ----------
        coords
            list astropy coordinates
        ID
            list of ID names (must be same length as coords). Will
            generate a running index if no ID is given
        rad
            radius of circular aperture to be extracted in arcsec.
            Default is 3.0
        multiprocess
            boolean flag to use multiprocessing. This will greatly
            speed up its operation as it will extract on 32 shots at
            time. But only use this when on a compute node. Use
            idev, a jupyter notebook, or submit the job as a single
            python slurm job.
        shotid
            list of integer shotids to do extractions on. By default
            it will search the whole survey except for shots located
            in the bad.shotlist file
        survey
            Survey you want to access. User note that HDR1 extractions
            are much slower compared to HDR2.
        tpmin
            Include only shots above tpmin. Default is 0.09.
        ffsky
            Use the full frame 2D sky subtraction model. Default is
            to use the local sky subtracted, flux calibrated fibers.
        
        Returns
        -------
        sources
            an astropy table object of source spectra for all input
            coords/ID that have spectra in the survey shots. There
            is one row per source ID/shotid observation.
    


Reading in the output - astropy FITS files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    sources[0:5]




.. raw:: html

    <i>Table length=5</i>
    <table id="table47890247148600" class="table-striped table-bordered table-condensed">
    <thead><tr><th>ID</th><th>shotid</th><th>wavelength [1036]</th><th>spec [1036]</th><th>spec_err [1036]</th><th>weights [1036]</th></tr></thead>
    <thead><tr><th></th><th></th><th>Angstrom</th><th>1e-17 erg / (Angstrom cm2 s)</th><th>1e-17 erg / (Angstrom cm2 s)</th><th></th></tr></thead>
    <thead><tr><th>int64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>360</td><td>20170130027</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8530224046204286 .. 0.8566766913047349</td></tr>
    <tr><td>360</td><td>20170224006</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8502338248829885 .. 0.8602435340384444</td></tr>
    <tr><td>360</td><td>20170326010</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8750953493900122 .. 0.8787945903307385</td></tr>
    <tr><td>372</td><td>20170126002</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8811780092325096 .. 0.8900289027924926</td></tr>
    <tr><td>372</td><td>20170322016</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8997470730626222 .. 0.9202416808656538</td></tr>
    </table>



Join your input and output table so you can match up any properties you
like

.. code:: ipython3

    output_table = join(input_cat, sources)

.. code:: ipython3

    output_table[0:5]




.. raw:: html

    <i>Table length=5</i>
    <table id="table47890247322648" class="table-striped table-bordered table-condensed">
    <thead><tr><th>ID</th><th>RA</th><th>DEC</th><th>WAVE</th><th>FLUX</th><th>FLUXE_L</th><th>FLUXE_H</th><th>z</th><th>shotid</th><th>wavelength [1036]</th><th>spec [1036]</th><th>spec_err [1036]</th><th>weights [1036]</th></tr></thead>
    <thead><tr><th></th><th>deg</th><th>deg</th><th>Angstrom</th><th>1e-17 erg / (cm2 s)</th><th>1e-17 erg / (cm2 s)</th><th>1e-17 erg / (cm2 s)</th><th></th><th></th><th>Angstrom</th><th>1e-17 erg / (Angstrom cm2 s)</th><th>1e-17 erg / (Angstrom cm2 s)</th><th></th></tr></thead>
    <thead><tr><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>182</td><td>150.05175</td><td>2.2376444444444443</td><td>4174.25</td><td>25.6</td><td>5.2</td><td>5.8</td><td>2.4337</td><td>20170222007</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8779764074573124 .. 0.8919157194039836</td></tr>
    <tr><td>182</td><td>150.05175</td><td>2.2376444444444443</td><td>4174.25</td><td>25.6</td><td>5.2</td><td>5.8</td><td>2.4337</td><td>20170202003</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.9882423520354104 .. 0.9140588857024189</td></tr>
    <tr><td>189</td><td>150.05504166666665</td><td>2.31525</td><td>4195.93</td><td>12.9</td><td>6.7</td><td>8.7</td><td>2.4515</td><td>20170331006</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8342767868348164 .. 0.6895230963441078</td></tr>
    <tr><td>194</td><td>150.05908333333332</td><td>2.2405861111111114</td><td>3997.41</td><td>61.0</td><td>4.3</td><td>4.9</td><td>2.2882</td><td>20181117010</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8559235841764234 .. 0.794577504985113</td></tr>
    <tr><td>194</td><td>150.05908333333332</td><td>2.2405861111111114</td><td>3997.41</td><td>61.0</td><td>4.3</td><td>4.9</td><td>2.2882</td><td>20181114020</td><td>3470.0 .. 5540.0</td><td>nan .. nan</td><td>nan .. nan</td><td>0.8211261682267719 .. 0.730172296780675</td></tr>
    </table>



.. code:: ipython3

    sel = output_table['FLUX'] > 10
    
    for row in output_table[sel][4:8]:
        plt.figure()
        wave_obj = row['WAVE']
        wave = row['wavelength']
        spec = row['spec']
        spec_err = row['spec_err']
        plt.errorbar(wave, spec, yerr=spec_err)
        plt.xlim(wave_obj-50, wave_obj+50)
        plt.xlabel('wave')
        plt.ylabel('spec')
        plt.title(row['ID'])




.. image:: output_22_0.png



.. image:: output_22_1.png



.. image:: output_22_2.png



.. image:: output_22_3.png


Examples of running get\_spec as a command line job:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run these commands from the command line by removing the "!"
command but be sure you are on a compute node by calling ``idev`` first
or submitting these as slurm jobs, one task per line.

.. code:: ipython3

    !hetdex_get_spec --ra 150.02548 --dec 2.087987 --ID cosmos_LAE --outfile cosmos_LAE

Speed things up using multiprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can speed up processes (by up to ~30x) with python multiprocessing
if you are working interactively in a notebook or in an idev session
(**NEVER FROM A LOGIN NODE**). Use the multiprocessing option with the
argument -mp True or --multiprocess True

.. code:: ipython3

    !hetdex_get_spec --multiprocess -ra 150.02548 -dec 2.087987 -id mptest -o mptest

Save output as individual astropy tables for each ID/shot combination:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer to save each spectra to a table, you can do this. We don't
recommend this for large numbers of data, but understand that sometimes
its easy to start with a readable table. The tables will be stored in
the files named ``spec_[ID]_[shotid].tab``

.. code:: ipython3

    !hetdex_get_spec  --multiprocess --single -ra 150.02548 -dec 2.087987 -id cosmos_lae

Getting all spectra at a specified RA/DEC in a specific OBSERVATION/SHOT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Perhaps you only want to focus on a specific shot. Then you can use the
-s argument to put the shotid either as an interger value 'YYYYMMDDOBS'=
20190104008 or as a str '20190104v009'. Note if you don't give an --ID
option the default is 'DEX'

This is a command line routine so remove the "!" if you are running in a
terminal.

.. code:: ipython3

    !hetdex_get_spec  -ra 8.86535 -dec 0.59352  -s 20190104008 -o 20190104008

This is particularly helpful if you plan to submit each shot as a
separate task. For this reason, I suggest changing the default --outfile
option to -o 20190104008 to create the output fits file 20190104008.fits

Work on a list of ID/RA/DECs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can either be a saved astropy table, or an space delimited text
file with 3 columns where the columns are ID, RA, DEC. If you want more
functionality with your input catalog, just talk to Erin. Note that
running this job will take about 30 minutes so only execute if you want
to wait around to explore the output.

.. code:: ipython3

    !cp /work/05350/ecooper/stampede2/3dhst/3dhst_input.cat .

.. code:: ipython3

    !hetdex_get_spec  --multiprocess -i '3dhst_input.cat' -o '3dhst'
