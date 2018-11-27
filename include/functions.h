//B145772
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h> 
#include <assert.h>
#include <string.h>
#include "globalVariables.h"


/*
 * Each processes has an instance of this struct.
 * The function initTopology(...) initialize this struct and
 * each process use this struct for number of reasons such as
 * to communicate with neighbors. 
 */
typedef struct topology
{
    int Np;
    int Mp;
    int rank;
    int right;
    int left;
    int up;
    int down;
    int coords[2];
}topology;

/**
 * Returns the fixed sawtooth value.
 * @param i
 * @param m
 * @return 
 */
double boundaryval(int i, int m);

/**
 * Implements the main algorithm for the image reconstruction.
 * Many things has been modified compared to the serial implementation.
 * However the most important changes are the function calls that implement:
 * 1) horizontal and vertical halo swaps,
 * 2) delta termination condition, 
 * 3) computation of sawtooth,
 * 4) print average value of the pixels
 * @param topo
 * @param edge
 * @param buf
 * @param old
 * @param new
 * @param dims
 */
void imageRecontruction( topology topo,double** edge,double** buf,double** old,double** new,int* dims);

/**
 * Returns true if the number is prime and false otherwise
 * @param num
 * @return 
 */
bool isNumberPrime(int num);

/**
 * Computes the boundary conditions.
 * I changed the serial implementation in order to compute correct
 * the boundary conditions.
 * @param topo
 * @param dims
 * @param old
 */// determins how ofter 
void computeBoundaryConditions(topology topo,int *dims,double **old);

/**
 * Print average value of the pixels if the target iteration 
 * equals with the current iteration. 
 * Every process send their sum of pixels using 
 * MPI_Reduce and process 0 calculates and prints the average value of pixels.
 * @param iter
 * @param targetIter
 * @param cnt
 * @param old
 * @param topo
 */
void printAverages(int iter,int targetIter,int* cnt ,double **old,topology topo );

/**
 * Exit the loop  across all the processes if the value of delta below 0.1.
 * Delta is calculated every 200 iterations but user also can change it using the
 * appropriate cmd line flags. I used MPI_Allreduce for the calculation of max delta across all processes.
 * @param topo
 * @param maxDelta
 * @param iter
 * @return 
 */
bool isTheLastIteration(topology topo,double maxDelta,int iter);


/**
 * Allocates arrays everything using provided function arralloc.
 * @param topo
 * @param masterbuf
 * @param buf
 * @param old
 * @param new
 * @param edge
 */
void allocations(topology topo,double ***masterbuf,double ***buf,double ***old,double ***new, double ***edge);

/**
 * Deallocates arrays.
 * @param masterbuf
 * @param buf
 * @param old
 * @param new
 * @param edge
 */
void deallocations(double **masterbuf,double **buf,double **old,double **new, double **edge);

/**
 * The user can optional set any of the following flags:
 * t: total prints of average pixel value.
 * i: input image
 * e: output image
 * d: 0 or 1 to enable or disable termination condition delta.
 * @param argc
 * @param argv
 */
void parseCmdLine(int argc, char *argv[]);

/**
 * Calculate max delta, if termination condition is enabled.
 * @param i
 * @param j
 * @param new
 * @param old
 * @param maxDelta
 */
void calculateMaxDelta(int i,int j,double** new,double** old,double* maxDelta);

#endif /* FUNCTIONS_H */

