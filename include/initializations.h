
#ifndef INITIALIZATIONS_H
#define INITIALIZATIONS_H
#include "../include/functions.h"

void initGlobalVariables();
void initAvgPrints(int* printAvgAtIter);
void initTopology(topology* topo,int worldSize,MPI_Comm* comm2d,int* dims);
void initDataTypes(topology topo);
void loadImage(topology topo,double** masterbuf,char* input);

#endif /* INITIALIZATIONS_H */

