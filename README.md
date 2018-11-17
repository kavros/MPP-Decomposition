## Build & Run
* Build project using: ```make```
* Before running any script or the project make sure that your directory is MPP-Decomposition.
* Run at the frontend of cirrus using: ```make run n=4```
* Run at the backend of cirrus using: ```make qsub```

# Testing
* Unit Tests
* Continuous integration 
* Output validation usin scripts

# Validation
* Validate output image and expected image ( delta need to be disabled) using: ```make validation```

# Acknowledgement
* I use the following tutorial in order to setup [travis-mpi](https://d-meiser.github.io/2016/01/10/mpi-travis.html)
