#!/bin/bash

for f in /home/pgranger/weightor/build/app/hd_AV*.root;
do
	fname=$(basename $f .root)
	ofname="${fname}_splines.root"
	logname="${fname}.log"
	echo "Processing $f"
	./build/bin/make_xsec_response_1d -m $f -w $f -o $ofname -selec numu > $logname 2>&1
done

