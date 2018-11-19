#!/bin/bash -x

for i in {1,2} #,2,3,4,8}
do
	qsub -v MPISIZE=$i,IMG="./data/input/edgenew768x768.pgm" -o "./data/qsub/$i.txt"  scripts/image.pbs  	

done
