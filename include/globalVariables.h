

#ifndef GLOBALVARIABLES_H
#define GLOBALVARIABLES_H
#include <mpi.h>
#define MAXITER   1500
#define PRINTFREQ  200                      
#define DELTAFREQ 50                        //deterimines how ofter our program calculates delta




bool isDeltaActivated;                      //determines if delta is activated or not
                                            //by default is deactivated

int totalAveragePrints;                     //determines the total 
                                            //print of the average value of the pixels

char *input;                                //input image file path and name
char *output;                               //output image file path and name

MPI_Datatype vectorMpxNp;                   //declaration of global data types for halo swaps
MPI_Datatype vectorMpxNP_N;    
MPI_Datatype vectorMpx1;

int M,N;                                    //M=width, N=height of the image


#endif /* GLOBALVARIABLES_H */

