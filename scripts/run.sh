#!/bin/bash 

nodes=1
for i in {1,2,3,4,8,16,32,64}
do
	
	if [ "$i" -eq "2" ]
	then
		nodes=2
		
	fi

	if [ "$i" -eq "4" ]
	then
		nodes=4
		
	fi
	

	
	qsub -q R380254 -l select=$nodes:ncpus=36 -v NPROC=$i,INPUT_IMG="./data/input/edgenew192x128.pgm",OUTPUT_IMG="./data/qsub/small/imagenew192x128_$i.pgm" -o "./data/qsub/small/$i.txt"  scripts/image.pbs  	
	qsub -q R380254 -l select=$nodes:ncpus=36 -v NPROC=$i,INPUT_IMG="./data/input/edgenew512x384.pgm",OUTPUT_IMG="./data/qsub/medium/imagenew512x384_$i.pgm" -o "./data/qsub/medium/$i.txt"  scripts/image.pbs  	
	qsub -q R380254 -l select=$nodes:ncpus=36 -v NPROC=$i,INPUT_IMG="./data/input/edgenew768x768.pgm",OUTPUT_IMG="./data/qsub/large/imagenew768x768_$i.pgm" -o "./data/qsub/large/$i.txt"  scripts/image.pbs  	
done
