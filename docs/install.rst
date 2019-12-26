Installation
============

To install hetdex_api you can grab the python package from pypi with a pip3 install. Be sure to use pip3 to install it for use in python3. Use the --user tag when installing on TACC.

.. code-block:: bash

   pip3 install hetdex_api --user

For TACC Users
--------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions and copy over Karl's .bashrc. This will set both your home and work directories to be readable to everyone on TACC.

.. code-block:: bash
   
   cp ~gebhardt/rsetup

Then install all required python packages:

.. code-block:: bash
   
   pip3 install --user tables
   pip3 install --user astropy
   pip3 install --pre --user astroquery 
   pip3 install --user ipywidgets
   pip3 install --user speclite
   
You should now exit the terminal/ssh shell to make sure all your permissions are set.
