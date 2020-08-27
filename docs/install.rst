Installation
============

For TACC Users 
---------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions so that your $WORK, $SCRATCH (on stampede2) and $DATA (on wrangler) directories are readable to everyone on TACC. Use at your own discretion, but this will allow you to share classifying work and notebooks with the team.

For wrangler:

.. code-block:: bash

   ssh username@wrangler.tacc.utexas.edu
   cd $STOCKYARD
   chmod -R a+rX .
   cd $DATA
   chmod -R a+rX .
   cdw

For stampede2:

.. code-block:: bash

   ssh username@stampede2.tacc.utexas.edu
   cd $STOCKYARD
   chmod -R a+rX .
   cd $SCRATCH
   chmod -R a+rX .
   cdw

A note about TACC data drives. $DATA on wrangler, $SCRATCH on stampede2 should host active computing and file creation. It is not subject to a data limit but it is also not backed up. Files untouched may be deleted by the system admin. Your $HOME drive is backed up but has limited storage. $WORK storage is limited to 1 Tb and this is across all computing clusters. For more info please read: 

https://portal.tacc.utexas.edu/tutorials/managingio

Then get the default bash script from TACC by running this script

.. code-block:: bash

   cd $HOME
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

Install Required Packages for hetdex-api
-----------------------------------------

First install the required python modules that are needed
for hetdex-api

.. code-block:: bash

   pip3 install -r /work/05350/ecooper/wrangler/hetdex_api/requirements.txt --user
   pip3 install --user --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex

Install hetdex-api: stable release version
----------------------------------------------

As of HDR2.1 release, a stable release of hetdex-api can now be pip installed from pypi 

.. code-block:: bash

   pip3 install hetdex_api --user --upgrade

Install hetdex-api: latest version
----------------------------------

If you want to be working with the most recent copy of hetdex-api please copy the git 
clone repository of hetdex_api. For anyone not on the core data team, we recommend you 
stick with the release versions.

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

You should not be doing any heavy computing or accessing more than one HDR product at a time on a login node. TACC users should use an interactive compute node on a shell by doing:

.. code-block:: bash

    idev -t 04:00:00

This will automatically switch you over to a compute node where you will have access to 48 cores per node and 128 GB of memory. Go nuts there!

Also, it is generally preferred that users store large files on their $DATA (on wrangler) and $SCRATCH (on stampede2) storage drive and any high I/O runs should be done on /tmp.


Jupyter Notebook Access
-----------------------

Both wrangler and stampede2 are setup for HDR access through hetdex-api. No configuration is needed after install. To access a notebook, in a browser go to:


https://vis.tacc.utexas.edu

Choose the 'all' queue mode under the wrangler or 'skx-dev' under stampede2.

We suggest you add symbolic links from your home to your $WORK and $SCRATCH or $DATA directories 
since a jupyter notebook node will open automatically in your $HOME directory. 

For example, 

.. code-block:: bash

   cd $HOME
   ln -s $WORK work-stampede2
   ln -s $SCRATCH scratch-stampede2

or on wrangler:

.. code-block:: bash

   cd $HOME
   ln -s $WORK work-wrangler
   ln -s $DATA data-wrangler 

This will allow you to go to your work directory when you log onto vis.

You can now open up a jupyter notebook and explore some of the notebooks in
hetdex-api/notebooks or just pop in some of the commands you see throughout this website.

Running a notebook from the command line
----------------------------------------

If accessing a node on https://vis.tacc.utexas.edu fails, you can also run this 
script from a terminal on wrangler:

.. code-block:: bash

    ~ecooper/bin/run_jupyter

This will launch from whatever directory you are working in. 
    
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


