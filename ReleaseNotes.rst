``hetdex_api`` release notes
==========================

Development version @ trunk
---------------------------

- The f50vals array in the SensitivityCube class has been
  replaced by an array of the 1 sigma noise called sigmas

- The old return_completeness method has been replaced with
  a new method where you must specify the detection 
  significance S/N cut as the final parameter. The old
  version of the method is available as return_completeness_hdr1
  for backward compatibility.

- You may now specify a conversion between signal divided by
  noise and detection significance using the conversion_poly
  parameter when initialising a SensitivityCube 
