qstat | grep -E tmat-\([0-9]+\).job -o | sort -V | uniq -c
#qstat | grep -E tmat-\([0-9]+\).job -o | sort -n -k6,7 | uniq -c
