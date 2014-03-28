#!/bin/bash

if [ -z "$1" ]; then
	echo usage: $0 dirindex
	echo NOTE: Runcopy is one less than the number specified
	exit
fi

RUNCOPY=$1
DIRNAME=k-run-$1

echo dirname: $DIRNAME runcopy: $RUNCOPY
pwd

mkdir $DIRNAME
cd ./$DIRNAME
pwd
python -m quantaccel.make_tmat_jobs $RUNCOPY --adapt weights
chmod +x submitter.sh

./submitter.sh

