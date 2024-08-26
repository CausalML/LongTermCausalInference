#!/bin/bash

for i in `seq 0 24`
do
	echo $i
	nohup r --no-save --args $i < surrogate.r &
done

nohup R --no-save --args $i < ridge.R &
