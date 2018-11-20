#include <stdbool.h>
#include "../include/initializations.h"
#include "../include/globalVariables.h"
#include "../include/pgmio.h"
#include <assert.h>
#include <mpi.h>

void initGlobalVariables()
{
    isDeltaActivated = false;
    totalAveragePrints=1;
    input = "./data/input/edgenew192x128.pgm";
    output="./data/output/imagenew192x128.pgm";
}

void initAvgPrints(int* printAvgAtIter)
{
    int i;
    if(totalAveragePrints <= 0)
    {
        assert(0);
    }
    //int printAvgAtIter[totalAveragePrints];
    
    for(i=0; i<totalAveragePrints; i++)
    {
        if(i==0)
        {
            printAvgAtIter[i] = 0;
        }
        printAvgAtIter[i] =i*(MAXITER/totalAveragePrints);
        //printf("%d, ",printAvgAtIter[i]);
        
    }
    //printf("\n");  
}

void initTopology(topology* topo,int worldSize,MPI_Comm* comm2d,int* dims)
{
    
    
    int disp=1,reorder=false,ndims =2;
    int period[ndims];
    MPI_Comm comm;
    comm  = MPI_COMM_WORLD;   
    
    dims[0]=0;
    dims[1]=0;
    if(isNumberPrime(worldSize))
    {
        if(M%worldSize ==0)
        {
            dims[1]=1;
        }
        else if(N%worldSize ==0)
        {
            dims[0]=1;
        }
    }
    
    MPI_Dims_create(worldSize,ndims,dims);
    //printf("(%d x %d) \n",dims[0],dims[1]);
    bool isSerialDecomp = (dims[0] == 1) && (dims[1] == 1 );
    bool isVerticalDecomp = (dims[0] > 1) && (dims[1] == 1 );
    bool isHorizontalDecomp = (dims[0] == 1) && (dims[1] > 1 );
    bool is2D               = (dims[0] > 1) && (dims[1] > 1 );
    if(isSerialDecomp || isVerticalDecomp)
    {
        period[0] =  false;
        period[1] =  false;
    }
    else if(is2D || isHorizontalDecomp)
    {
        period[0] = false;
        period[1] = true;   //cyclic on y axis
    }
    else
    {
        
        MPI_Finalize();
        assert(0);
    }
    
    MPI_Cart_create(comm,ndims,dims,period,reorder,comm2d);
    MPI_Cart_shift(*comm2d,0,disp,&(topo->left),&(topo->right));
    MPI_Cart_shift(*comm2d,1,disp,&(topo->down),&(topo->up));
    MPI_Comm_rank(*comm2d,&(topo->rank));    
    MPI_Cart_coords(*comm2d, topo->rank, 2, topo->coords);
    
    
    //check if columns and rows are divisible
    bool isTopologyDivisible = ( (N%dims[1]) == 0) && ( (M % dims[0]) ==0 );
    if(isTopologyDivisible == false)
    {
        if(topo->rank==0)
            printf("Error: Number of processes is not supported\n");
        MPI_Finalize();
        exit(-1);
    }
    assert(dims[0]!=0);
    assert(dims[1]!=0);
    int Mp = (M/dims[0]);
    int Np = (N/dims[1]);
    
    topo->Mp = Mp;
    topo->Np = Np;

}


void initDataTypes(topology topo)
{
    int count =topo.Mp;          
    int blocklength=topo.Np;
    int stride = topo.Np;
    MPI_Type_vector(count, blocklength, stride, MPI_DOUBLE, &vectorMpxNp);
    MPI_Type_commit(&vectorMpxNp);
    
    
    int count_1 =topo.Mp;          
    int blocklength_1=topo.Np;
    int stride_1 = N;
    MPI_Type_vector(count_1, blocklength_1, stride_1, MPI_DOUBLE, &vectorMpxNP_N);
    MPI_Type_commit(&vectorMpxNP_N);
    
    
    count = topo.Mp;
    blocklength = 1;
    stride = topo.Np+2;
    MPI_Type_vector(count, blocklength, stride, MPI_DOUBLE, &vectorMpx1);
    MPI_Type_commit(&vectorMpx1);
    
}


void loadImage(topology topo,double** masterbuf,char* input)
{
    if(topo.rank == 0)
    {
        //printf("---master\n");
        //printf("Processing %d x %d image\n", M, N);
        //printf("Number of iterations = %d\n", MAXITER);
        printf("\nReading <%s>\n", input);
        pgmread(input, &masterbuf[0][0], M, N);
    }
}