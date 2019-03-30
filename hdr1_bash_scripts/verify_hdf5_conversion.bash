#
#  Verify all flux limit
#  HDF5s are present

# List of dates and shots to verify
INLIST=$1

echo "If nothing is printed then everything worked: "

while read line
do

    date=`echo $line | awk 'BEGIN {FS=" "}{print $1}'`
    shot=`echo $line | awk 'BEGIN {FS=" "}{print $2}'`
    dateshot="${date}v${shot}"
    fn="${date}v${shot}_sensitivity_cube.h5"

    # Test folder exists
    if [ ! -f $fn ]; then
        echo "Missing! $fn"
        echo "bash ../hdr1_code/HETDEX_API/hdr1_bash_scripts/generate_sensitivity_cube_hdf5s.bash ${dateshot}" >> hdf5_repeat.run
        continue 
    fi
   
done < $INLIST

