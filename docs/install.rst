Installation
============

For TACC Users
--------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions so that your /work directories are readable to everyone on TACC. Use at your own discretion, but this will allow you to share classifying work and notebooks with the team.

.. code-block:: bash

   ssh wrangler
   cd $STOCKYARD
   chmod -R a+rX .
   cdw

Then get the default bash script from TACC by running this script

.. code-block:: bash

   /usr/local/startup_scripts/install_default_scripts

Then open your .bashrc and uncomment this line:
::

   umask 022

and add in the following module loads/unloads:
::

   module unload python
   module unload python2
   module load intel/18.0.2
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
   pip3 install --user matplotlib==2.0.0
   pip3 install --user --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex

Copy the git clone repository of hetdex_api 

.. code-block:: bash
		
   git clone https://github.com/HETDEX/hetdex_api.git

Then pip3 install with the -e parameter to update as the repository evolves

.. code-block:: bash
   
   pip3 install -e hetdex_api --user --upgrade

Install Elixer
--------------

We also recommend that you install elixer:

.. code-block:: bash

    git clone https://github.com/HETDEX/elixer.git

.. code-block:: bash

   pip3 install -e elixer --user --upgrade


Compute Nodes on TACC
---------------------

You should not be doing any heavy computing or accessing more than one HDR2 product at a time on a login node. TACC users should use an interactive compute node on a shell by doing:

.. code-block:: bash

    idev -t 04:00:00

This will automatically switch you over to a compute node where you will have access to 48 cores per node and 128 GB of memory. Go nuts there!

Also, it is generally preferred that users store large files on their /data storage drive and any high I/O runs should be done on /tmp.

If you would like to use a jupyter notebook, run this script:

.. code-block:: bash

    ~ecooper/bin/run_jupyter

This will launch from whatever directory you are working in. 
    
If on stampede2 (not relevant for wrangler), one final suggestion is to add a link from your home to your work directory. For example, I would do:

.. code-block:: bash
   
   cd
   ln -s /work/05350/ecooper/ work-stampede2

This will allow you to go to your work directory when you log onto vis.

You can now open up a jupyter notebook and explore some of the notebooks in 
hetdex-api/notebooks or just pop in some of the commands you see throughout this website. 

In your favourite browser goto https://vis.tacc.utexas.edu and log onto stampede2. Choose the 
jupyter notebook option and pick the skx-dev queue. 


For Contributors
----------------

To contribute to github

.. code-block:: bash
   
   git add filename
   git commit -m "Reason for update or file creation"
   git push

Please ask to become a member of HETDEX organization on github once you have an account. Please branch your development if you are doing major code work.

If you want to build the documentation, you can install the necessary packages by adding ``[doc]`` to
the package name when you install, e.g.

.. code-block:: bash
   
   pip3 install -e hetdex_api[doc] --user --upgrade




