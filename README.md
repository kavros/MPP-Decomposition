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
* The following scripts are inside the scripts directory:
  * run.sh: run my experiments on qsub
  * image.pbs: qsub script
  * get_mpich.sh: install mpi on travis
  * runAndValidateLocal.sh: run 18 tests (with and without delta enabled) and validate that their output images are correct
  * generateGraphs.py : generates speedup and total time graphs based on qsub results
  * qsubResultsValidation.py: validates that the results of qsub (images + logfiles) are correct.

## Scripts Instructions
* Before running any script or the project make sure that your current directory is MPP-Decomposition.
* Run experiments at the backend of cirrus, validate their results and generate graphs(located at ./data/graphs using the following commands:
```
./scripts/run.sh 
python scripts/qsubResultsValidation.py
python  scripts/generateGraphs.py 
```
* Run experiments and validate the results local using:
```
./scripts/runAndValidateLocal.sh
```

## Testing
* For testing I used continuous integration(travis) and scripts.(see Scripts Description)

## References
* I used the following tutorial in order to setup [travis-mpi](https://d-meiser.github.io/2016/01/10/mpi-travis.html)
* I used an external library for command line parsing called [Argtable](https://github.com/argtable/argtable3)
