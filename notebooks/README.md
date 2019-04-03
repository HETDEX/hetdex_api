## First, if you don't have an account on TACC:

To get one, go to:

https://portal.tacc.utexas.edu/

Then, click:

Create a TACC account

Finally, send Karl your username so he can add you to the HETDEX group.

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

3. clone the GitHub HETDEX_API to somewhere on ``/work`` (Not necessary if you are pip installing from HDR1/software/HETDEX_API)

```
git clone https://github.com/grzeimann/HETDEX_API.git
```

4. Go into the HETDEX_API directory and install with pip. This will install the API & dependencies NOT including the ELiXer API.

```
cd HETDEX_API
pip install --user .
```
or you can just install from the HETDEX_API installed in HDR1 (this will likely be out of date soon).

```
pip install --user /work/03946/hetdex/hdr1/software/HETDEX_API
```

4.1 To use the ELiXer API, you will need to install or reference ELiXer. Please see /work/03946/hetdex/hdr1/software/elixer/docs/Elixer_readme.pdf (Installation section on page 3)

5. Also install other dependent python modules

```
pip install --user tables
pip install --user astropy
pip install --user --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex
pip install --user /work/03946/hetdex/hdr1/software/vdrp
```

6. goto https://vis.tacc.utexas.edu/ 

7. log onto stampede2 ipython jupyter notebook (normal queue worked for me)

8. cd stampede2/path_to_HETDEX_API/notebooks
