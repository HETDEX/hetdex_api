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

These engineering shots were accidentally labelled as science shots. Should be ignored for main science.
20180511v015  189.25340   62.20436  324.60000
20180308v002   87.00297   56.10039  303.60000
20180309v001   78.00181   56.09970  303.48000
20180309v002   83.00311   56.09934  303.48000
20180309v003   89.00241   56.09940  303.54000
20180309v004   97.00384   56.09990  303.50000
20181205v012   41.85607   19.36482  117.50000



