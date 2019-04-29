# List of known issues

## Bad shots

A list of badshots that should be removed according to science goals are contained in:

badshots.list

Any detections contained in the shot are removed upon initialization of Detections() API. Users not using this call should remove the detections in their code.

These engineering shots were accidentally labelled as science shots. Should be ignored for main science. 
* 20180511v015 189.25340 62.20436 324.60000 
* 20180308v002 87.00297 56.10039 303.60000 
* 20180309v001 78.00181 56.09970 303.48000 
* 20180309v002 83.00311 56.09934 303.48000 
* 20180309v003 89.00241 56.09940 303.54000 
* 20180309v004 97.00384 56.09990 303.50000 
* 20181205v012 41.85607 19.36482 117.50000

These shots have bad astrometry:
* 20190106008 - was flagged as borderline accepted by RAs due to clouds in 3rd dither, but fails due to astrometry (could be related to the clouds in 3rd dither. Its automatically moved in detections with the TP = 0.0813 so not adding to the bad list for now.

These shots were rejected due to Human Error (ie. any dither has an 'H' ra_flag)

[20170124013, 20170306018, 20170306020, 20170306021, 20170306022,
       20170421016, 20170427006, 20170615009, 20171025009, 20180210009,
       20180421009, 20180519007, 20190208035]
      

These shots were rejected due to Equipment Error ie. RA_flag = E

[20170126003, 20170127007, 20170922025, 20170922028, 20171014007,
       20171016114, 20171029022, 20171222013, 20180213015, 20180310002,
       20180310007, 20180313015, 20180313016, 20180313017, 20180316015,
       20180316016, 20180317012, 20180320009, 20180517015, 20180714007,
       20190104017]
