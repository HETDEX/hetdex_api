Installation
============

For TACC Users
--------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions so that both your home and work directories to be readable to everyone on TACC. Use at your own discretion.

.. code-block:: bash

   ssh wrangler or stampede2
   chmod a+rx ../username
   cdw
   chmod a+rx ../*
   chomd a+rx ../../username
   chmod a+rx ../wrangler
   chmod a+rx ../stampede2

Then get the default bash script from TACC by running this script

.. code-block:: bash

   /usr/local/startup_scripts/install_default_scripts

Then uncomment this line:
::

   umask 022

and add in the following module loads/unloads:
::

   module unload python
   module unload python2
   module load python3
   alias python='python3'

and add in the following line to your $PATH:
::

   export PATH=$PATH:$HOME/bin:$HOME/.local/bin

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

One final suggestion is to add a link from your home to your work directory. For example, I would do:

.. code-block:: bash
   
   cd
   ln -s /work/05350/ecooper/ work-stampede2

This will allow you to go to your work directory when you log onto vis.

You can now open up a jupyter notebook and explore some of the notebooks in 
hetdex-api/notebooks or just pop in some of the commands you see throughout this website. 

In your favourite browser goto <vis.tacc.utexas.edu> and log onto stampede2. Choose the 
jupyter notebook option and pick the skx-dev queue. 

For Contributors
----------------

To contribute to github

.. code-block:: bash
   
   git add filename
   git commit -m "Reason for update or file creation"
   git push

Please ask to become a member of HETDEX organization on github once you have an account. Please branch your development if you are doing major code work.
