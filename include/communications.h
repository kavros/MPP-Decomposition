
#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H
#include "functions.h"
#include "pgmio.h"

/**
 * Using this function rank 0 send to each processes their part of the image.
 * All processes except rank 0 are receiving their part of the image.
 * @param masterbuf
 * @param buf
 * @param topo
 * @param worldSize
 * @param comm2d
 */
void scatter(double** masterbuf,double** buf,topology topo,int worldSize,MPI_Comm comm2d);

/**
 * All processes are sending their part of the image after they finish with the calculations to rank 0.
 * Rank 0 compose all the parts in to masterbuf array.
 * @param topo
 * @param masterbuf
 * @param buf
 * @param comm2d
 * @param output
 * @param worldSize
 */
void gather(topology topo,double** masterbuf,double** buf,MPI_Comm comm2d,char* output,int worldSize);

/**
 * Implements the horizontal halo swaps using derived data types.
 * @param old
 * @param topo
 */
void halloSwapsHorizontal(double** old,topology topo);

/**
 * Implements vertical halo swaps without derived data types.
 * @param old
 * @param topo
 */
void halloSwapsVertical(double** old,topology topo);

/**
 * Rank 0 writes image to the file.
 * The user of the api must first call gather function.
 * @param topo
 * @param masterbuf
 */
void saveImage(topology topo, double** masterbuf);



#endif /* COMMUNICATIONS_H */

