# Bash Scripts used to Generate the files

If you're not on the data reduction team this might not be 
very interesting. It's just any BASH scripts we used to set
commands running in the reduction.


################################################################

SHOT FILE CREATION ON *****STAMPEDE2******

Instructions to generate all shot HDF5 files using the lists
hdr1/reduction/hdr1.scilist
hdr1/reduction/hdr1.callist 

In hdr1/reduction/data, execute script contained in 
/work/03946/hetdex/hdr1/software/HETDEX_API/hdr1_bash_scripts

cd /work/03946/hetdex/hdr1/reductions/data/
/work/03946/hetdex/hdr1/software/HETDEX_API/make_shot_run.sh

(This requires the scripts jobsplitter_stampede2 and mkhdf5_stampede2.sh)
This will create many run_shot_XX.run and run_shot_XX.slurm files, 
to be executed as:

sbatch run_shot_1.slurm
sbatch run_shot_2.slurm
etc...

These take some time to enter the queue, but once they are in,
it will take about 5-20 minutes to finish.

NOTE: to convert on wrangler please use build_shot_slurmfiles.py

