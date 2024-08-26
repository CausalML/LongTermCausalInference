#!/bin/bash

for eta in 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6
do
	nohup R --no-save --args $eta < pipeline.R &
	nohup R --no-save --args $eta < pipeline_ridge.R &
done
