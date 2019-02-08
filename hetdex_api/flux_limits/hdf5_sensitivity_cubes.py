"""



"""

import logging
import tables as tb
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

class SensitivityCubeHDF5Container(object):
    """
    Handle accessing and writing the sensitivity
    cubes stored together in an HDF5 container.
    Arguments passed to the init are passed on
    to the open_file function of tables.    

    Parameters
    ----------
    filename : string
        the filename of the HDF5


    Attributes
    ----------
    h5file : tables:File 
        the tables File object

    """
    def __init__(self, filename, **kwargs):

        self.h5file = tb.open_file(filename, **kwargs)
        self.filename = filename

    def add_sensitivity_cube(self, datevshot, ifuslot, scube, flush = False):
        """
        Add a sensitivity cube to the HDF5 file

        Parameters
        ----------
        datevshot : str
            the shot the IFU belongs 
            to
        ifuslot : str
            the IFU identifier
        scube : SensitivityCube
            the sensitivity cube object to
            add
        """

        # Try to get this shot from the database
        # if it doesn't exist add it
        try:
            shot = self.h5file.get_node(self.h5file.root, datevshot)
        except tb.NoSuchNodeError:
            shot =  self.h5file.create_group(self.h5file.root, datevshot)

        # Store this in a compressible array
        array = self.h5file.create_carray(shot, ifuslot, obj=scube.f50vals, title="50% Detection Limits")
       
        # Store the header as an attribute
        array.attrs.header = scube.header
        array.attrs.wavelengths = scube.wavelengths
        array.attrs.alphas = scube.alphas


    def extract_ifu_sensitivity_cube(self, ifuslot, datevshot=None):
        """
        Extract the sensitivity cube
        from IFU (ifuslot). If multiple
        shots are saved in the file specify
        which to use with datevshot

        Parameters
        ----------
        ifuslot : string
            the IFU slot to extract
        datevshot : string (Optional)
            the datevshot if multiple
            shots are stored in the 
            HDF5. If None then test
            that one shot is present 
            and return the IFU
            for that

        Returns
        -------
        scube : hetdex_api.flux_limits.sensitivity_cube:SensitivityCube
            the sensitivity cube
        """ 

        # Use first shot if dateshot not specified
        if not datevshot:            
            shots = self.h5file.list_nodes(self.h5file.root)
            nshots = len(shots)
            if nshots > 1:
                logging.warn("""Datevshot not specified but multiple shots in file!
                                Using first in file.""")

            shot = shots[0]
        else:
            shot = self.h5file.get_node(self.h5file.root, name=datevshot)


        # Now get the desired IFU
        ifu = self.h5file.get_node(shot, name=ifuslot)

        # Extract the data we need for a sensitivity cube
        header = ifu.attrs.header
        wavelengths = ifu.attrs.wavelengths
        alphas = ifu.attrs.alphas
        f50vals = ifu.read() 

        return SensitivityCube(f50vals, header, wavelengths, alphas)

    def flush(self):
        """ Write all alive leaves to disk """
        self.h5file.flush()

    def close(self):
        """ Close the file and destroy the object """
        self.h5file.close()
        logging.info("Closed {:s}".format(self.filename))

    def __enter__(self):
        """ Added to support using the `with` statement """
        return self


    def __exit__(self, type_, value, traceback):
       """ Support tidying up after using the `with` statement """
       self.close()
    

