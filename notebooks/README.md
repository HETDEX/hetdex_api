## First, if you don't have an account on TACC:

To get one, go to:

https://portal.tacc.utexas.edu/

Then, click:

Create a TACC account

Finally, send Karl your username so he can add you to the HETDEX group. FYI, you may need to wait up to 24 hours for your password and group authentication to get approved.

## Instructions to get the notebooks working on stampede2:

```
ssh username@stampede2.tacc.utexas.edu
```

In your home (technically this is ``home1/``)​

1. remove any paths to Karl’s anaconda in your .bashrc

2. create a symbolic link to connect to your work directory:
```
ln -s /work/magicnumber/username/stampede2/ stampede2
```

3. pip install HETDEX_API

Option 1: pip install HETDEX_API. You can just install from the HETDEX_API installed in HDR1 (this will likely be out of date soon).

```
pip install --user /work/03946/hetdex/hdr1/software/HETDEX_API
```

Option 2: clone the GitHub HETDEX_API to somewhere on ``/work`` (Not necessary if you are pip installing from HDR1/software/HETDEX_API)

```
git clone https://github.com/grzeimann/HETDEX_API.git
```

Go into the HETDEX_API directory and install with pip. This will install the API & dependencies NOT including the ELiXer API.

```
cd HETDEX_API
pip install --user .
```

As we update HETDEX_API, this latter option will be needed. For updates, 

```
cd HETDEX_API
git pull
pip install --user --upgrade .
```

4 To use the ELiXer API, you will need to install or reference ELiXer. Please see /work/03946/hetdex/hdr1/software/elixer/docs/Elixer_readme.pdf (Installation section on page 3)

```
pip install photutils --user
```

5. Also install other dependent python modules

```
pip install --user tables
pip install --user astropy
pip install --user pickle
pip install --user ipywidgets
pip install --user --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex
```

You don't need to install vdrp unless you plan to use the reduction scripts located within it.

Note: do not worry about pip version warnings. I have had issues after upgrading pip so I recommend you remain with version 9.0.3 for now. If you do upgrade and have issues, email Erin (erin@astro) and she can help you revert back.

6. Copy over the notebooks to a working directory on /work. These create some example output files so you probably want these somewhere you can delete and don't make your git clone too messy.

```
cdw (alias to cd into your working directory)
mkdir notebookplay
cp /work/03946/hetdex/hdr1/software/HETDEX_API/notebooks/* notebookplay/
```

7. goto https://vis.tacc.utexas.edu/ 

8. log onto stampede2 ipython jupyter notebook (development queue is fastest but will only be active for 2 hours. Normal queue will give you longer)

9. cd stampede2/notebookplay
