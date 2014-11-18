#!/bin/bash

exit_stat=0
for fn in maccelerator/testing/test_*
do
    echo "$fn"
    python $fn $1
    exit_stat=$(($exit_stat + $?))
    echo "Cumulative exit status: $exit_stat"
done

exit $exit_stat
