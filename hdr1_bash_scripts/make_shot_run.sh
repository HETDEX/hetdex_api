hdr1_path='/work/03946/hetdex/hdr1'
awk '{print "/work/03946/hetdex/hdr1/software/HETDEX_API/hdr1_bash_scripts/mkhdf5_stampede2.sh",$1,$2}' $hdr1_path/reduction/hdr1.scilist > run_shot
awk '{print "/work/03946/hetdex/hdr1/software/HETDEX_API/hdr1_bash_scripts/mkhdf5_stampede2.sh",$1,$2}' $hdr1_path/reduction/hdr1.callist >> run_shot
$hdr1_path/software/HETDEX_API/hdr1_bash_scripts/jobsplitter_stampede2 run_shot 192 1 '00:30:00'
