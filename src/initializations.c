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
    
    
    int disp;
    const int ndims =2;
    
    int period[ndims];
    int reorder;
    MPI_Comm comm;

    comm  = MPI_COMM_WORLD;          
    reorder = false;
    disp = 1;  
    if(worldSize == 1)
    {
        topo->Mp = M;
        topo->Np = N;
        dims[0] = 0;
        dims[1] = 0;
        period[0] =  false;
        period[1] =  false;
        
        MPI_Dims_create(worldSize,ndims,dims);
        MPI_Cart_create(comm,ndims,dims,period,reorder,comm2d);
        MPI_Cart_shift(*comm2d,0,disp,&(topo->left),&(topo->right));
        MPI_Cart_shift(*comm2d,1,disp,&(topo->down),&(topo->up));
        MPI_Comm_rank(*comm2d,&(topo->rank));
    }
    else if( isNumberPrime(worldSize) == true )
    {   
        //finds a 1D decomposition
        
        //printf("1D\n");
        if( (M % worldSize == 0) )
        {
            //vertical decomposition            
            //printf("vertical\n");
            topo->Mp = M/worldSize;
            topo->Np = N;
            dims[0] = 0;
            dims[1] = 1;        //one dimension on y axis
            period[0] = false;
            period[1] = false;
        }
        else if( (N % worldSize == 0) )
        {
            //printf("horizontal\n");
            //horizontal decomposition
            topo->Mp = M;
            topo->Np = N/worldSize;
            dims[0] = 1;        //one dimension on x axis
            dims[1] = 0;
            period[0] = false;
            period[1] = true;   //cyclic on y axis
        }
        else
        {
            //printf("Combination of image and thread number is not supported.\n");
            MPI_Finalize();
            exit(1);
        }
        MPI_Dims_create(worldSize,ndims,dims);
        MPI_Cart_create(comm,ndims,dims,period,reorder,comm2d);
        MPI_Cart_shift(*comm2d,0,disp,&(topo->left),&(topo->right));
        MPI_Cart_shift(*comm2d,1,disp,&(topo->down),&(topo->up));
        MPI_Comm_rank(*comm2d,&(topo->rank));
        
    }
    else
    {
        //find 2D decomposition
        
        dims[0]=0;
        dims[1]=0;
        period[0]=false;
        period[1]=true;
        
        //TODO:find a smart way for assignment of dims[0] and dims[1]
        
        MPI_Dims_create(worldSize,ndims,dims);
        MPI_Cart_create(comm,ndims,dims,period,reorder,comm2d);
        MPI_Cart_shift(*comm2d,0,disp,&(topo->left),&(topo->right));
        MPI_Cart_shift(*comm2d,1,disp,&(topo->down),&(topo->up));
        MPI_Comm_rank(*comm2d,&(topo->rank));
        //printf("2d (%dx%d) \n",dims[0],dims[1]);
        bool isTopologyDivisible = ( (N%dims[0]) == 0) && ( (M % dims[1]) ==0 );
        if(  isTopologyDivisible   )
        {
            topo->Mp=M/dims[0];
            topo->Np=N/dims[1];
        }
        else
        {
            //printf("Combination of image and thread number is not supported.\n");
            MPI_Finalize();
            exit(1);
        } 
        
    }
    MPI_Cart_coords(*comm2d, topo->rank, 2, topo->coords);
    
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