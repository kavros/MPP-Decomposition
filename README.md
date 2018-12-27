## Build & Run [![Build Status](https://travis-ci.com/kavros/MPP-Decomposition.svg?token=FpiEbZxSSgCpKyAubrTW&branch=master)](https://travis-ci.com/kavros/MPP-Decomposition)
* Before running any script or the project make sure that your directory is MPP-Decomposition.
* Build project using: ```make```
* Run at the frontend of cirrus using 4 processes and delta disabled: ```make run n=4```
* Validate output image without delta: ```make validation```
* Run at the frontend of cirrus using 4 processes and delta enabled: ```make run_delta n=4```
* Run at the backend of cirrus using 4 processes: ```make qsub n=4```

## Run Description
* We can run the project in 2 ways:
  * without command line arguments using :```mpirun -n numberOfprocesses ./build/image ```
  * and wih command line arguments using: ```mpirun -n numberOfprocesses ./build/image -i inputImagePath -e outputImagePath -t totalNumber -d 0or1 ``` 
* Using the following command line arguments we can change the input and the output. Also, we can enable/disable termination condition and determine the total prints of the average value of pixels.

Argument | Description
---      | ---
-i | determines the input image path and name. 
-e | determines the output image path and name. 
-t | determines the total number that we are going to print average value of the pixels.
-d | we can use 0 or 1 to disable or enable delta correspondingly.By default is disabled.
-h | provide a help message.

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
* I used the [matplotlib](https://matplotlib.org/) which is a python library in order to plot my graphs.
