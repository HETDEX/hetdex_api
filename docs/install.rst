Installation
============

Pip Install
-----------

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
   pip3 install --user ipyevents
   pip3 install --user astrowidgets
   pip3 install --user jsonschema==3.1.1
   pip3 install --user sep
   pip3 install --user specutils
   pip3 install --user photutils

Copy the git clone repository of hetdex_api 

.. code-block:: bash
		
   git clone https://github.com/HETDEX/hetdex_api.git


Then pip3 install with the -e parameter to update as the repository evolves

.. code-block:: bash
   
   pip3 install -e hetdex_api --user

Be sure to remove any PYTHONPATH to HETDEX_API in your .bashrc. And then add the line

.. code-block:: bash

   export PATH="~/.local/bin:$PATH"

in your .bashrc so the command line quick entry points work.

For Contributors
----------------

To contribute to github

.. code-block:: bash
   
   git add filename
   git commit -m "Reason for update or file creation"
   git push

Please ask to become a member of HETDEX organization on github once you have an account. Please branch your development if you are doing major code work.
