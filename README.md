# HETDEX_API

This repository contains software develeloped for the HETDEX Data Releases. The team has adopted the HDF5 data format for the release. This standard file format allows for hierarchically structured data to be flexibly delivered in a single container. For a basic briefing of HDF5 files, we suggest checking out the notebook tutorials at https://github.com/tomkooij/scipy2017.  This repository contains both the code to create HDF5 container files and also consists of a library of python modules to easily access data within the HDF5 files.   

To access and contribute to the HETDEX_API, the code library may installed from github in a working directory on TACC by executing:

git clone https://github.com/grzeimann/HETDEX_API.git


To contribute to a file or add a new file:

git add filename

git commit -m "Reason for update or file creation"

git push

The directory "notebooks" contains a list of jupyter notebooks to familiarize the user with HDF5 file access and querying. It also contains several case uses of API modules developed by the team for scientific analysis.


The main code to create HDF5 files contained in the DR1 release are:

create_survey_hdf5.py - this creates the main survey info file which contains general shot and exposure info. A HETDEX shot always has 3 dithered exposures in it, while a parallel shot typically only has one.

create_hdf5file.py - this is the code to all info related to a single shot, including the raw fits image data from each exposure, flux calibrated fibers
