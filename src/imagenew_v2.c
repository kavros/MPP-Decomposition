#include <stdio.h>
#include <stdlib.h>
#include "../include/functions.h"
#include "../include/pgmio.h"
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
    
    initGlobalVariables();   
    parseCmdLine(argc,argv);
    
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    pgmsize(input,&M,&N);
    
    dims[0]=0;dims[1]=0;    
    initTopology(&topo,worldSize,&comm2d,dims);

    allocations( topo,&masterbuf,&buf,&old,&new,&edge);
    
    initDataTypes(topo);
    
    loadImage(topo,masterbuf,input);
    
    if(topo.rank == 0)
        start = MPI_Wtime();
    
    scatter( masterbuf, buf, topo, worldSize, comm2d);
    
    imageRecontruction( topo, edge, buf, old, new,dims);
    
    gather( topo, masterbuf,buf,comm2d,output,worldSize);
    
    
    if(topo.rank == 0)
    {
        end = MPI_Wtime();
        fprintf(stdout,"total time is %f (sec), woldSize is %d \n",end-start,worldSize);
    }
    
    MPI_Finalize();
    
    deallocations(masterbuf,old,new,edge,buf);

    return 0;
} 
