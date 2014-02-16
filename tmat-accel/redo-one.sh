#!/bin/bash

if [ -z "$1" ]; then
	echo usage: $0 dirindex
	echo NOTE: Runcopy is one less than the number specified
	exit
fi

RUNCOPY=$1
DIRNAME=h-run-$1

echo dirname: $DIRNAME runcopy: $RUNCOPY
pwd

cd ./$DIRNAME
pwd
python ../../src/quantaccel/make_tmat_jobs.py $RUNCOPY
python ../../src/quantaccel/submit_tmat_jobs.py $RUNCOPY
chmod +x submitter.sh
./submitter.sh

cd ../
