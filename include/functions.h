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
 * 
 * @param topo
 * @param edge
 * @param buf
 * @param old
 * @param new
 * @param dims
 */
void imageRecontruction(topology topo,double** edge,double** buf,double** old,double** new,int* dims);

/**
 * Returns true if the number is prime and false otherwise
 * @param num
 * @return 
 */
bool isNumberPrime(int num);




void computeBoundaryConditions(topology topo,int *dims,double **old);
void printAverages(int iter,int targetIter,int* cnt ,double **old,topology topo );
bool isTheLastIteration(topology topo,double maxDelta,int iter);



void allocations(topology topo,double ***masterbuf,double ***buf,double ***old,double ***new, double ***edge);
void deallocations(double **masterbuf,double **buf,double **old,double **new, double **edge);
void parseCmdLine(int argc, char *argv[]);
void calculateMaxDelta(int i,int j,double** new,double** old,double* maxDelta);



#endif /* FUNCTIONS_H */

