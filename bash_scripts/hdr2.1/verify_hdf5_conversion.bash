#
#  Verify all flux limit
#  HDF5s are present

# List of dates and shots to verify
INLIST=/data/00115/gebhardt/getcen/fulllist
BADLIST=/data/05350/ecooper/hdr2.1/survey/badshots.list

echo "If nothing is printed then everything worked: "

while read line
do

    date=`echo $line | awk 'BEGIN {FS=" "}{print $1}'`
    shot=`echo $line | awk 'BEGIN {FS=" "}{print $2}'`
    dateshot="${date}v${shot}"
    fn="${date}v${shot}_sensitivity_cube.h5"

    # Test HDF exists
    if [ ! -f $fn ]; then
 
        # test if a bad shot
        bad=`grep "${date}${shot}" $BADLIST`

        if [ "$bad" == "" ] ; then
            echo "Missing! $fn"
            echo $dateshot >> missing
        fi
        continue 
    fi
   
done < $INLIST

