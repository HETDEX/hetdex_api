## Instructions to get the notebooks working on stampede2:

```
ssh stampede2
```

In your home (technically this is ``home1/``)​

1. remove any paths to Karl’s anaconda in your .bashrc

2. create a symbolic link to connect to your work directory:
```
ln -s /work/magicnumber/username/stampede2/ stampede2
```

3. clone the GitHub HETDEX_API to somewhere on ``/work``

```
git clone https://github.com/grzeimann/HETDEX_API.git
```

4. Go into the HETDEX_API directory and install with pip. This will install the API & dependencies

```
cd HETDEX_API
pip install --user .
```

5. goto https://vis.tacc.utexas.edu/ 

6. log onto stampede2 ipython jupyter notebook (normal queue worked for me)

7. cd stampede2/path_to_HETDEX_API/notebooks
