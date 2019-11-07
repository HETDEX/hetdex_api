# HETDEX_API

## Introduction 

This repository contains software developed for the HETDEX Data Release. The team has adopted the HDF5 data format for the release. This standard file format allows for hierarchically structured data to be flexibly delivered in a single container. For a basic briefing of HDF5 files, we suggest checking out the notebook tutorials at [https://github.com/tomkooij/scipy2017](https://github.com/tomkooij/scipy2017).  

This repository contains both the code to create HDF5 container files and also consists of a library of python modules to easily access data within the HDF5 files. To access and contribute to the HETDEX_API, the code library may be installed from github in a working directory on TACC by executing:

```
git clone https://github.com/HETDEX/HETDEX_API.git
```

## Description of Code

The main code to create HDF5 files contained in the DR1 release are:

* create_survey_hdf5.py - this creates the main survey info file which contains general shot and exposure info. A HETDEX shot always has 3 dithered exposures in it, while a parallel shot typically only has one.

* create_shot_hdf5.py - this is the code to all info related to a single shot, including the raw fits image data from each exposure, flux calibrated fibers with their associated astrometric information, response curves, and other calibration information pertaining to an individual shot. It is used along with create_astrometry_hdf5.py and create_cal_hdf5.py which append the astrometry information to an existing shot file.

* create_detect_hdf5.py - this is the code to encode the detections database into an HDF5 container. Included in the HDF5 file are 1d flux-calibrated, PSF-weighted spectra for each source detection.

* hetdex_api/ - contains a number of functions to interact on the Survey, Shot, Fibers, Images and Detections classes. 

* hetdex_api/flux_lmits - deals with creating and accessing HDF5 containers for the datacubes that contain the flux limit. Also computes
the average flux limits across shots and IFUs

* hetdex_api/config.py - configuration file to link to HETDEX DR1 files

* notebooks/ - contains a number of jupyter notebooks to help HETDEX users learn to interact with the data and the API to the data. See readme within this directory on how to start these notebooks.

## Installation of HETDEX_API and flux limit tools

The tools to deal with flux limits need to be installed, to do this, and install the library and binaries in your
home directory, run this command in the top directory (the directory with setup.py)

```
pip install --user .

```

If you want to install them to a root directory, or an Anaconda environment, you can drop the ``user`` flag. If you want to install the library and binaries in a particular location you can use the ``--prefix`` flag e.g.

```
pip install --user --prefix /path/to/install .
```

You will need to add the ``/bin/`` directory created by this command to your path.

## Instructions to contributors

To contribute to a file or add a new file:

```
git add filename

git commit -m "Reason for update or file creation"

git push
```

