# HETDEX_API

## Introduction 

This repository contains software developed for the HETDEX Data Release. The team has adopted the HDF5 data format for the release. This standard file format allows for hierarchically structured data to be flexibly delivered in a single container. For a basic briefing of HDF5 files, we suggest checking out the notebook tutorials at [https://github.com/tomkooij/scipy2017](https://github.com/tomkooij/scipy2017).  

The latest and up-to-date documentation is available here [https://hetdex-api.readthedocs.io/en/latest/install.html](https://hetdex-api.readthedocs.io/en/latest/install.html).

This repository contains both the code to create HDF5 container files and also consists of a library of python modules to easily access data within the HDF5 files. To access and contribute to the HETDEX_API, the code library may be installed from github in a working directory on TACC by executing:

```
git clone https://github.com/HETDEX/hetdex_api.git
```

## Description of Code

The main code to create HDF5 files contained in the DR1 release are:

* h5tools/ - contains scripts to generate the H5 container files

* hetdex_tools/ - higher level programs to extract 1D and 2D spectra 

* hetdex_api/ - contains a number of functions to interact on the Survey, Shot, Fibers, Images and Detections classes. 

* hetdex_api/config.py - configuration file to link to HETDEX data release file paths

* notebooks/ - contains a number of jupyter notebooks to help HETDEX users learn to interact with the data and the API to the data. See readme within this directory on how to start these notebooks.

## Installation of HETDEX_API and flux limit tools

The tools to deal with flux limits need to be installed, to do this, and install the library and binaries in your
home directory, run this command in the top directory (the directory with setup.py)

```
pip install -e hetdex_api --user

```

If you want to install them to a root directory, or an Anaconda environment, you can drop the ``user`` flag. If you want to install the library and binaries in a particular location you can use the ``--prefix`` flag e.g.

```
pip install --user --prefix /path/to/install .
```

You will need to add the ``/bin/`` directory created by this command to your path.

If you want to build the documentation yourself, then you will also want to install the packages
for that by adding ``[doc]`` to the name of the package when you install, e.g.

```
pip install -e hetdex_api[doc] --user
``` 

## Instructions to contributors

To contribute to a file or add a new file:

```
git add filename

git commit -m "Reason for update or file creation"

git push
```

