//B145772
#include <stdbool.h>
#include "../include/initializations.h"
#include "../include/globalVariables.h"
#include "../include/ioutils.h"
#include <assert.h>
#include <mpi.h>

void initGlobalVariables()
{
    isDeltaActivated = false;
    totalAveragePrints=1;
    
//    output = "./data/output/imagenew192x128.pgm";
    output = "./data/output/imagenew0192x0128.dat";
    input = "./data/input/cedgefiles/cedgenew0192x0128.dat";
}

void initAvgPrints(int* printAvgAtIter)
{
    int i;
    if(totalAveragePrints <= 0)
    {
        assert(0);                        //validate that totalAveragePrints
                                           //cannot be negative or zero.
    }
    
    //initialize an array that holds the iterations that we 
    //will print the average pixel values.
    for(i=0; i<totalAveragePrints; i++)
    {
        if(i==0)
        {
            printAvgAtIter[i] = 0;
        }
        printAvgAtIter[i] =i*(MAXITER/totalAveragePrints);
        
    }
}

void initTopology(topology* topo,int worldSize,MPI_Comm* comm2d,int* dims)
{
    
    
    int disp=1,reorder=false,ndims =2;
    int period[ndims];
    MPI_Comm comm;
    comm  = MPI_COMM_WORLD;   
    
    dims[0]=0;
    dims[1]=0;
    
    //if the number of processes is prime then the decomposition can only be 1D.
    //Thus we check if horizontal or vertical decomposition is suitable.
    if(isNumberPrime(worldSize))
    {
        if(M%worldSize ==0)         //check if columns are divisible with the num of processes.
        {
            dims[1]=1;
        }
        else if(N%worldSize ==0)
        {
            dims[0]=1;
        }
    }
    
    MPI_Dims_create(worldSize,ndims,dims);

    //determines if cyclic allowed.This is useful for finding neighboors correct.
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
        //printf("\nReading <%s>\n", input);
        //pgmread(input, &masterbuf[0][0], M, N);
        ioread(input, &masterbuf[0][0], M*N);
    }
}

void loadImageInParallel(topology topo,double** buf,char* input,MPI_Comm comm2d)
{
    MPI_File fh;
    MPI_Status status;
    MPI_Datatype subarray;
    int sizes[2]    = {M,N};  
    int subsizes[2] = {topo.Mp,topo.Np};
    int my_coords[2];
    MPI_Cart_coords(comm2d, topo.rank, 2,my_coords );
    int startY = my_coords[1]*topo.Np;
    int startX = my_coords[0]*topo.Mp;
    int starts[2]   = {startX,startY};
    printf("startX = %d , startY = %d, subsizesX = %d,subsizesY=%d input = %s\n",
            starts[0],
            starts[1],
            subsizes[0],
            subsizes[1],
            input);

    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subarray);
    MPI_Type_commit(&subarray);
            
    //createfilename(input, "cinput", M, N, -1);
    if(MPI_File_open(comm2d,input,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh) != MPI_SUCCESS)
    {
        printf("Open error on rank %d\n", topo.rank);
        return;
    }
    
    if(MPI_File_set_view(fh,0,MPI_DOUBLE,subarray,"native",MPI_INFO_NULL) != MPI_SUCCESS)
    {
      printf("View error on rank %d\n", topo.rank);
      return;
    }
    
    if(MPI_File_read_all(fh, &buf[0][0],topo.Mp*topo.Np,MPI_DOUBLE,&status) != MPI_SUCCESS)
    {
      printf("Read error on rank %d\n", topo.rank);
      return;
    }
    //printf("%f %f %f\n",buf[0][0],buf[0][1],buf[0][2]);
    if (MPI_File_close(&fh) != MPI_SUCCESS)
    {
      printf("Close error on rank %d\n", topo.rank);
    }

}
