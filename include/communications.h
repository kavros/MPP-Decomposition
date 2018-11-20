
#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H
#include "functions.h"
#include "pgmio.h"

void scatter(double** masterbuf,double** buf,topology topo,int worldSize,MPI_Comm comm2d);
void gather(topology topo,double** masterbuf,double** buf,MPI_Comm comm2d,char* output,int worldSize);
void halloSwapsHorizontal(double** old,topology topo);
void halloSwapsVertical(double** old,topology topo);



#endif /* COMMUNICATIONS_H */

