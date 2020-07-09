#
#
#
#  Verify the flux limit cubes ran and
#  no left over uncompleted temporary 
#  directories exist
#

# List of dates and shots to verify
INLIST=$1

echo "If nothing is printed then everything worked: "

while read line
do

    date=`echo $line | awk 'BEGIN {FS=" "}{print $1}'`
    shot=`echo $line | awk 'BEGIN {FS=" "}{print $2}'`
    folder="${date}v${shot}"

    # Test folder exists
    if [ ! -d $folder ]; then
        echo "Non-existant folder! $folder"
        continue 
    fi
   
    # Test stuff in the folder 
    test=$(ls ${folder}/flim/*.fits 2>/dev/null | wc -l) 
    if [ "$test" == 0 ]; then 
        echo "No results! "$folder
    fi

    # Test for left over directories from finished runs
    failures=`ls -l "${folder}/flim/" | grep "^d" | awk 'BEGIN {FS = " "}{print $9}'`

    if [ "$failures" != "" ]; then
        for fail in $failures
        do
            echo "Unfinished file! $fail"
            if [[ $fail == *"multi"* ]]; then
                echo "Suggested solution"
                ifuslot=`echo $fail | awk 'BEGIN {FS="_"}{print "ifu"$4}'`
                echo "vdrp_setup_flim --commit --cores_per_job 2 --ifu_slots $ifuslot -- $date $shot"
            fi
        done
    fi
done < $INLIST


# Command to submit lots
# for dir in `cat failures | awk 'BEGIN {FS=" "}{print $3}'`; do cd $dir/flim/; ls *.slurm; cd -; done 
