#
#
# Replace failed IFUs
# and take note of the HDFs
# that need to be regenerated

datevshot=""

for file in `ls ifu_redo/*flim.fits`
do
    ifu=`basename $file`

    # find original location of file (even if it's been moved already)
    original=`ls run*/${ifu} 2>/dev/null` 
    if [ "$original" == "" ]
    then
        original=`ls run*/${ifu}.old | sed 's/.old//g'`
    fi

    # save datevshots with issues
    tdatevshot="${ifu:0:12}"
    if [ "$datevshot" != "$tdatevshot" ]
    then
        dir_=`dirname $original`
        datevshot=$tdatevshot
        echo $dir_
        echo $datevshot >> $dir_/hdfs_to_redo
    fi
 
    # save original, copy new version
    mv $original "${original}.old3"
    cp $file $original
done
