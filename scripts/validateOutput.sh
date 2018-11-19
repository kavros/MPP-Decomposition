#!/bin/bash

# run the application using various number of processes
# and compare the output images with the image from the serial execution.(delta is disabled)
for i in {1,2,3,4,8}
do
	mpirun -n $i ./build/image -o "./data/output/imagenew192x128_$i.pgm" 
	diff  "data/output/imagenew192x128_$i.pgm" data/output/imagenew192x128_expectedOutput.pgm
done

mpirun -n 1 ./build/image -d 1 -o ./data/output/imagenew192x128_d_1.pgm 
for i in {2,3,4,8}
do
	mpirun -n $i ./build/image -d 1 -o "./data/output/imagenew192x128_d_$i.pgm" 
	diff ./data/output/imagenew192x128_d_1.pgm  "./data/output/imagenew192x128_d_$i.pgm" 
done

