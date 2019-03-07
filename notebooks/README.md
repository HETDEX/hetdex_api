Instructions to get the notebooks working on stampede2:

ssh stampede2

In your home (technically this is home1/)​

1. remove any paths to Karl’s anaconda in your .bashrc

2. Pip install tables and astropy 

pip install --user tables
pip install --user astropy

4. create a symbolic link to connect to your work directory:

ln -s /work/magicnumber/username/stampede2/ stampede2

5. clone the GitHub HETDEX_API to somewhere on /work

git clone https://github.com/grzeimann/HETDEX_API.git

6. goto https://vis.tacc.utexas.edu/ 

7. log onto stampede2 ipython jupyter notebook (normal queue worked for me)

8. cd stampede2/path_to_HETDEX_API/notebooks
