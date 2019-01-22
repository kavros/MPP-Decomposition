//B145772
#include <stdio.h>
#include <stdlib.h>
#include "../include/functions.h"
#include "../include/pgmio.h"
#include "../include/ioutils.h"
#include "../include/initializations.h"
#include "../include/communications.h"

int main (int argc, char *argv[])
{
    MPI_Comm comm2d;
    int worldSize;
    topology topo;

    double **old,**new,**edge,**buf,**masterbuf;
    int dims[2];
    double start=0,end=0;
    
    initGlobalVariables();                      //initialize global variables 
                                                //with default values
    
    parseCmdLine(argc,argv);                    //parsing command line values
    
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    iosize(input,&M,&N);
    //printf("M=%d,N=%d\n",M,N);
    
    dims[0]=0;dims[1]=0;    
    initTopology(&topo,worldSize,&comm2d,dims);     //initialization of virtual topology 
                                                    

    allocations( topo,&masterbuf,&buf,&old,&new,&edge); //allocation of buffers
    
    initDataTypes(topo);                        //initialization of derived data types.
    
    loadImage(topo,masterbuf,input);            //process zero load image from file
    
    if(topo.rank == 0)
        start = MPI_Wtime();                    //starts timer
    
    scatterUsingSubArray( masterbuf, buf, topo, worldSize, comm2d);  //process 0 sends 
                                                        // part of the image to 
                                                        //to other process
                                                        // and the others are        
                                                        //receive their parts.
    
    imageRecontruction( topo, edge, buf, old, new,dims);    //every process run 
                                                            //the reconstruction
                                                            //algorithm
    
    
    gatherUsingSubArray(masterbuf,buf,topo,worldSize,comm2d);//process 0
                                                            //receive parts of the
                                                            //reconstructed image
                                                            //and compose it.
    
    if(topo.rank == 0)
    {
        end = MPI_Wtime();                                  //stops timer.
        fprintf(stdout,"total time is %f (sec), woldSize is %d \n",end-start,worldSize);
    }

    saveImage(topo,masterbuf);                          //save recontructed image 
                                                    // to the file.
    
    deallocations(masterbuf,old,new,edge,buf);          //dealocation of buffers
    
    MPI_Finalize();
    
    

    return 0;
} 
