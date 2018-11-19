## Build & Run
* Before running any script or the project make sure that your directory is MPP-Decomposition.
* Build project using: ```make```
* Run at the frontend of cirrus using 4 processes and delta disabled: ```make run n=4```
* Validate output image without delta: ```make validation```
* Run at the frontend of cirrus using 4 processes and delta enabled: ```make run_delta n=4```
* Run at the backend of cirrus using 4 processes: ```make qsub n=4```

## Run without Makefile
* We can run the project in 2 ways:
* without command line arguments using :```mpirun -n numberOfprocesses ./build/image ```
* and wih command line arguments using: ```mpirun -n numberOfprocesses ./build/image -i inputImagePath -e outputImagePath -t totalNumber -d 0or1 ``` 
  * inputImagePath: is the edge image path. 
  * outputImagePath: is the output image after the calculations
  * totalNumber: is the total number that we are going to print average value of the pixels
  * 0or1: we can use 0 or 1 to disable or enable delta correspondingly.By default is disabled.

## Scripts Description
* run.sh: 
* image.pbs:
* get_mpich.sh:
* image.pbs:
* validateOutput.sh:


## Scripts Instructions
* Before running any script or the project make sure that your current directory is MPP-Decomposition.
* Run experiments at the backend of cirrus and generate graphs(located at ./data/graphs using the following commands:
```
./scripts/run.sh 
python  scripts/generateGraphs.py 
```

## Testing
* For testing I used continuous integration(travis) and output validation using scripts.

## Acknowledgement
* I use the following tutorial in order to setup [travis-mpi](https://d-meiser.github.io/2016/01/10/mpi-travis.html)
