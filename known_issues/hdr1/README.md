## Overview

This directory contains lists of bad pixels, amps, shots, detectids that have been found both prior and after the HDR1 data release. In some cases, the lists (bad pixels/amps) were removed upon ingestion into the release, but as we add more cases, it is important for every user to mask out this bad data.

Note on Detections - the ingested Detections Database was intentionally left to contain a high number of sources to finesse with cuts afterwards. For example, many nearby stars and nearby galaxies (Hbeta) are included, but can easily be excluded with a cut on the continuum levels. We also have lists of known bad shots and detectids that need to be exluded. These will be growing lists. We also very liberal in our S/N and Chi2 cuts in order to not lose any possible sources. It will be an active effort to determine which ones can be removed automatically in software and which ones first need to be evaluated by eye.


Note on Detections - the ingested Detections Database was intentionally left to contain a high number of sources to finesse with cuts afterwards. For example, many nearby stars and nearby galaxies (Hbeta) are included, but can easily be excluded with a cut on the continuum levels. We also have lists of known bad shots and detectids that need to be exluded. These will be growing lists. We also very liberal in our S/N and Chi2 cuts in order to not lose any possible sources. It will be an active effort to determine which ones can be removed automatically in software and which ones first need to be evaluated by ey


## Description of files:

badpix.new  
* copied from /work/03946/hetdex/hdr1/calib
* this contains a list of known bad pixels that are additional to the pixel flats
* format is multiframe X1 X2 Y1 Y2 
* X1, X2, Y1, Y2 are the X and Y pixel values in the ds9 image of the multiframe files
* they need to be flipped and subtracted by 1 to match images in the shot H5 Images class
* moving forward we will likely convert this and work in these dimensions
           
badamps.list
* time dependent list of bad amps
* taken from /work/03946/hetdex/hdr1/reduction/badamps.list
* we expect this to be a growing list
                
baddetects.list
* list of detectids that should be ignored from the detections database
* there was a minor issue in ingesting that lead to the creation of this list
* we expect this list to grow as we will want to put in false positives and
* other artifacts in this list
       
badshots.list
* list of shots that should not be used for analysis despite being found grouped as science frames

## Handling of bad data

If you use the API to access the data we will make every effort to use these working lists to remove bad data. The API for the detections database: detections.py will automatically remove detections in badamps.list, badshots.list and baddetects.list upon initialization.

For this reason you should clone the HETDEX_API github and update it every so often (We will send out updates to hetdex-general when major updates are done). Then you will need to pip install --upgrade --user . within the location of the git clone.

## List of known issues

Please refer to knownissues.md contained within this folder for a list of known issues with the API and Data Release 1


