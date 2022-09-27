Installation
============

Recommended Data Access
-----------------------

For new users, we recommend you access HETDEX Data via our JupyterLabs hosted at TACC  and at SciServer (https://sciserver.mpe.mpg.de). JupyterLabs are pre-built computing environments where many common python modules, including HETDEX specific modules are pre-installed. We also provide a number of Jupyter notebook tutorials to allow the user to become familiar with HETDEX-API and the HETDEX data model. Internal users will have access to 1Tb of storage and extra HPC resources.

A TACC account is needed to access the TACC system: https://portal.tacc.utexas.edu/account-request. If you are part of the internal HETDEX collaboration, email Erin and Karl to set up your account. For public users, you can still create an account and you will have access to our public JupyterLab and public data only. No long-term storage will be provided but you will be able to learn how to access the HETDEX Public Catalogs and view some HETDEX visualization tools.

Instructions for the TACC JupyterLab:

1. Goto: https://jupyter.tacc.cloud
2. Internal Users Pick "HETDEX Internal" Spawner Option. Click HPC option. Public users pick the HETDEX Public Option"
3. Internal Users Check out this PDF for more info: http://web.corral.tacc.utexas.edu/hetdex/HETDEX/shared/HETDEX-JupyterLab-Instructions.pdf


Easy Install for stampede2
---------------------------

To get going on stampede2 you need to run a script to install hetdex-api and set up your permissions. You will need to ssh into stampede2. This step will require setting up dual authentication at TACC. Instructions are here: https://portal.tacc.utexas.edu/tutorials/multifactor-authentication. 

.. code-block:: bash

   ssh username@stampede2.tacc.utexas.edu
   ~ecooper/bin/tacc_setup

This script will copy over hetdex-api notebooks to the directory `hetdex-notebook-tutorials` for you to try out some hetdex-api tutorials. It will also make a symlink to your /work directory on stampede2 so when you launch a notebook on vis.tacc.utexas.edu, you will be able to navigate to your work directory.

Once you have executed this script logout and log back in to update your environment.

You can now log onto https://vis.tacc.utexas.edu and run Jupyter notebooks. Choose the options resource ='Stampede2', Session Type='Jupyter Notebook', Queue='skx-dev'.
 

For TACC Users 
---------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions so that your $WORK, $SCRATCH (on stampede2) directories are readable to everyone on TACC. Use at your own discretion, but this will allow you to share classifying work and notebooks with the team.

To do these steps you can either ssh in (for example, with the terminal app on Mac OS) or you can access a terminal from vis.tacc.utexas.edu. For cluster specific connection details, please see `Jupyter Notebook Access`. 

For stampede2:

.. code-block:: bash

   ssh username@stampede2.tacc.utexas.edu
   cd $STOCKYARD
   chmod -R a+rX .
   cd $SCRATCH
   chmod -R a+rX .
   cdw

A note about TACC data drives: $SCRATCH on stampede2 should host active computing and file creation. It is not subject to a data limit but it is also not backed up. Files untouched may be deleted by the system admin. Your $HOME drive is backed up but has limited storage. $WORK storage is limited to 1 Tb and this is across all computing clusters. For more info please read: 

https://portal.tacc.utexas.edu/tutorials/managingio

Copy over Erin's .bashrc script to set up your environment.

.. code-block:: bash

    cp ~ecooper/.bashrc $HOME/.bashrc

Now you need to exit the terminal session by completely closing your ssh connection or your vis.tacc.utexas.edu notebook session. Then log back in to finish the install.

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

This will allow you to go to your work directory when you log onto vis.

You can now open up a jupyter notebook and explore some of the notebooks in
hetdex-api/notebooks or just pop in some of the commands you see throughout this website.
We recommend you copy over the notebook tutorials to explore in your local directory.

.. code-block:: bash

    cp -r /work/05350/ecooper/wrangler/hetdex_api/notebooks $WORK/

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


