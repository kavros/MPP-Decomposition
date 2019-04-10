//B145772
#ifndef INITIALIZATIONS_H
#define INITIALIZATIONS_H
#include "../include/functions.h"

/**
 * Initialize the global variables.
 * By default delta is disabled but user can enabled by passing -d 1 as command line flag.
 * Also total prints for the average pixel by default is 1 but user can change it
 * command line flag -t.
 */
void initGlobalVariables();

/**
 * Initialize an array which holds the iterations that we print the average pixel value.
 * @param printAvgAtIter
 */
void initAvgPrints(int* printAvgAtIter);

/**
 * Initialize topology values using Cartesian grid topology.
 * @param topo
 * @param worldSize
 * @param comm2d
 * @param dims
 */
void initTopology(topology* topo,int worldSize,MPI_Comm* comm2d,int* dims);

/**
 * Initialize the main(global) derived data types of the application.
 * @param topo
 */
void initDataTypes(topology topo);

/**
 * Rank 0 loads images using this function.
 * @param topo
 * @param masterbuf
 * @param input
 */
void loadImage(topology topo,double** masterbuf,char* input);

void loadImageInParallel(topology topo,double** buf,char* input,MPI_Comm comm2d);
#endif /* INITIALIZATIONS_H */

