#!/bin/bash -x

#for i in {1,2,3,4,8}
for i in {1,2,3,4,8,16,24,32,64}
do
	qsub -v MPISIZE=$i,INPUT_IMG="./data/input/edgenew192x128.pgm",OUTPUT_IMG="./data/qsub/small/imagenew192x128_$i.pgm" -o "./data/qsub/small/$i.txt"  scripts/image.pbs  	
	qsub -v MPISIZE=$i,INPUT_IMG="./data/input/edgenew512x384.pgm",OUTPUT_IMG="./data/qsub/medium/imagenew512x384_$i.pgm" -o "./data/qsub/medium/$i.txt"  scripts/image.pbs  	
	qsub -v MPISIZE=$i,INPUT_IMG="./data/input/edgenew768x768.pgm",OUTPUT_IMG="./data/qsub/large/imagenew768x768_$i.pgm" -o "./data/qsub/large/$i.txt"  scripts/image.pbs  	
done

