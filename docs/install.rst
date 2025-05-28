Installation
============

JupyterHub Access
-----------------

For new users, we recommend you access HETDEX Data via our JupyterLabs hosted at TACC and at SciServer (https://sciserver.mpe.mpg.de). JupyterLabs are pre-built computing environments where many common python modules, including HETDEX specific modules are pre-installed. We also provide a number of Jupyter notebook tutorials to allow the user to become familiar with HETDEX-API and the HETDEX data model. Internal users will have access to 1Tb of storage and extra HPC resources.

A TACC account is needed to access the TACC system: https://portal.tacc.utexas.edu/account-request. If you are part of the internal HETDEX collaboration, email Erin and Karl to set up your account. For public users, you can still create an account and you will have access to our public JupyterLab and public data only. No long-term storage will be provided but you will be able to learn how to access the HETDEX Public Catalogs and view some HETDEX visualization tools.

Instructions for the TACC JupyterLab:

1. Goto: https://jupyter.tacc.cloud
2. Internal Users Pick "HETDEX Internal" Spawner Option. Public users pick the HETDEX Public Option"
3. Internal Users Check out this PDF for more info: http://web.corral.tacc.utexas.edu/hetdex/HETDEX/shared/HETDEX-JupyterLab-Instructions.pdf

All software is pre-built so you do not need to do anything further. Note that all directories in the JupyterLab environment are temporary aside from the /work directory which is where your persistent storage is stored (for internal users only). Public users should download data as they generate it in the notebook.
   
For Internal TACC Users
-----------------------

If this is your first time on a TACC cluster we recommend a few setup steps. First set your permissions so that your $WORK, $SCRATCH (on stampede2) directories are readable to everyone on TACC. Use at your own discretion, but this will allow you to share classifying work and notebooks with the team.

For ls6:

.. code-block:: bash

   ssh username@ls6.tacc.utexas.edu
   cd $STOCKYARD
   chmod -R a+rX .
   cd $SCRATCH
   chmod -R a+rX .
   cdw

A note about TACC data drives: $SCRATCH on ls6 should host active computing and file creation. It is not subject to a data limit but it is also not backed up. Files untouched may be deleted by the system admin. Your $HOME drive is backed up but has limited storage. $WORK storage is limited to 1 Tb and this is across all computing clusters. For more info please read this to understand best practices at TACC.

https://docs.tacc.utexas.edu/basics/conduct/

Not necessary anymore but you can consider copying over Erin's .bashrc script to set up your environment. It will uncomment the line "unmask 022" to allow for shared read-access of your files.

.. code-block:: bash
    cp ~ecooper/.bashrc $HOME/.bashrc

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


