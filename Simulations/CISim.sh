#!/bin/bash
fin=Configs/config10.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	studies=$(echo ${line} | awk "{print \$1}")
	rate=$(echo ${line} | awk "{print \$2}")
	alpha=$(echo ${line} | awk "{print \$3}")
	beta=$(echo ${line} | awk "{print \$4}")

	# Run simulation.
	Rscript Rscripts/ci_sim.R --studies ${studies} --rate ${rate} --alpha ${alpha} --beta ${beta} --reps 10 --mc 10 --out "Results/";
done