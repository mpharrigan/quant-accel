qstat | grep -E tmat-\([0-9]+\).job -o | sort -V | uniq -c
