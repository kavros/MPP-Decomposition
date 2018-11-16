/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/arralloc.h"
#include "../include/pgmio.h"
#include <mpi.h>
#include <stdbool.h> 
#include <assert.h>
#define MAXITER   1500
#define PRINTFREQ  200

typedef struct topology
{
    int Np;
    int Mp;
    int rank;
    int right;
    int left;
    int up;
    int down;
    int coords[2];
}topology;

double boundaryval(int i, int m);
void imageRecontruction(topology topo,double** edge,double** buf,double** old,double** new,int N,int* dims);
bool isNumberPrime(int num);
void initialization(topology* topo,int worldSize,int M,int N,MPI_Comm* comm2d,int* dims);
void scatter(double** masterbuf,double** buf,topology topo,int worldSize,char* input,MPI_Comm comm2d,MPI_Datatype masterbufType,MPI_Datatype vectorMpxNp,int M,int N);
void gather(topology topo,double** masterbuf,double** buf,MPI_Datatype vectorMpxNp,MPI_Datatype masterbufType,int M,int N,MPI_Comm comm2d,char* output,int worldSize);
void createDataTypes(topology topo,MPI_Datatype* vectorMpxNp,MPI_Datatype* masterbufType,int N);
void computeBoundaryConditions(topology topo,int *dims,double **old,int N);
void halloSwapsHorizontal(double** old,topology topo);

int main (void)
{
    
    MPI_Comm comm2d;
    int worldSize;
    topology topo;
    int M,N;
    char *input,*output;
    input = "./images/input/edgenew192x128.pgm";
    output="./images/output/imagenew192x128.pgm";
    
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    pgmsize(input,&M,&N);
    int dims[2];
    dims[0]=0;dims[1]=0;
    
    initialization(&topo,worldSize,M,N,&comm2d,dims);
    //printf("p%d left=%d,right=%d up=%d down=%d \n",topo.rank,topo.left,topo.right,topo.up,topo.down);
    //printf("p%d Mp=%d, Np=%d\n",topo.rank,topo.Mp,topo.Np);
    //printf("M=%d,N=%d\n",M,N);
    //printf("p%d (%d,%d)\n",topo.rank,topo.coords[0],topo.coords[1]);
    
    double **old,**new,**edge;
    double **buf,**masterbuf;
    masterbuf = (double**) arralloc(sizeof(double),2,M,N);
    old = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    new = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    edge = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    buf = (double**) arralloc(sizeof(double),2,topo.Mp,topo.Np);
    
    MPI_Datatype vectorMpxNp;    
    MPI_Datatype masterbufType;    
    
    
    createDataTypes(topo,&vectorMpxNp,&masterbufType, N);
    
    scatter( masterbuf, buf, topo, worldSize, input, comm2d, masterbufType,vectorMpxNp ,M, N);
    
    imageRecontruction( topo, edge, buf, old, new,N,dims);
    
    gather( topo, masterbuf,buf,vectorMpxNp,masterbufType, M, N, comm2d,output,worldSize);
    
    MPI_Finalize();
    
    free(masterbuf);
    free(old);
    free(new);
    free(edge);
    free(buf);
    return 0;
} 

void createDataTypes(topology topo,MPI_Datatype* vectorMpxNp,MPI_Datatype* masterbufType,int N)
{
    int count =topo.Mp;          
    int blocklength=topo.Np;
    int stride = topo.Np;
    MPI_Type_vector(count, blocklength, stride, MPI_DOUBLE, vectorMpxNp);
    MPI_Type_commit(vectorMpxNp);
    
    
    int count_1 =topo.Mp;          
    int blocklength_1=topo.Np;
    int stride_1 = N;
    MPI_Type_vector(count_1, blocklength_1, stride_1, MPI_DOUBLE, masterbufType);
    MPI_Type_commit(masterbufType);
    
}

void gather(topology topo,double** masterbuf,double** buf,MPI_Datatype vectorMpxNp,MPI_Datatype masterbufType,int M,int N,MPI_Comm comm2d,char* output,int worldSize)
{
    MPI_Request request;
    MPI_Status status;
    int i,j;
     if(topo.rank!=0)
    {
        
        MPI_Isend(&buf[0][0],1,vectorMpxNp,0,0,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
    }
    
    
    if(topo.rank==0)
    {
        for(i=0;i<M;i++)
            for(j=0;j<N;j++)
                masterbuf[i][j]=0;                  //clean masterbuf for testing
        
        int coords[2];
        for(i=1;i<worldSize;i++)
        {
            MPI_Cart_coords(comm2d, i, 2, coords);
            int startY = coords[1]*topo.Np;
            int startX = coords[0]*topo.Mp;
            //printf("p %d, coords[0]=%d,coords[1]=%d ,startY=%d, startX=%d\n",i,coords[0],coords[1],startY,startX);                             
            
            MPI_Irecv(&masterbuf[startX][startY],1,masterbufType,i,0,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,&status);
        }
        
        for(i=0;i<topo.Mp;i++)
            for(j=0;j<topo.Np;j++)
                masterbuf[i][j] = buf[i][j];
        
        pgmwrite(output, &masterbuf[0][0], M, N);
    }
}


void scatter(double** masterbuf,double** buf,topology topo,int worldSize,char* input,MPI_Comm comm2d,MPI_Datatype masterbufType,MPI_Datatype vectorMpxNp,int M,int N)
{
    int i,j;
    MPI_Request request;
    MPI_Status status;
    if(topo.rank == 0)
    {
        //printf("---master\n");
        //printf("Processing %d x %d image\n", M, N);
        //printf("Number of iterations = %d\n", MAXITER);
        
        printf("\nReading <%s>\n", input);
        pgmread(input, &masterbuf[0][0], M, N);
        
        
        int coords[2];
        
        for(i=1; i < worldSize; i++)
        {
            MPI_Cart_coords(comm2d, i, 2, coords);
            int startY = coords[1]*topo.Np;
            int startX = coords[0]*topo.Mp;
            //printf("p %d, startY=%d, startX=%d\n",i,startY,startX);
            
            MPI_Isend(&masterbuf[ startX][startY], 1, masterbufType, i, 0, MPI_COMM_WORLD,&request);
            
            MPI_Wait(&request,&status);
        }
        for(i=0;i<topo.Mp;i++)
            for(j=0;j<topo.Np;j++)
                buf[i][j] = masterbuf[i][j];
    }
    
    
    
    if(topo.rank !=0)
    {
        MPI_Irecv(&buf[0][0], 1, vectorMpxNp, 0, 0, MPI_COMM_WORLD, &request); 
        MPI_Wait(&request,&status);
        
    }
    
}

bool isNumberPrime(int num)
{
    assert(num > 0);
    int i;
    if(num == 1)
    {
        return false;
    }
    for(i=2; i<num; i++)
    {
        if(num%i == 0)
        {
            return false;
        }
    }
    return true;
}

void initialization(topology* topo,int worldSize,int M,int N,MPI_Comm* comm2d,int* dims)
{
    
    int disp;
    const int ndims =2;
    
    int period[ndims];
    int reorder;
    MPI_Comm comm;

    comm  = MPI_COMM_WORLD;          
    reorder = false;
    disp = 1;                   
    if( isNumberPrime(worldSize) == true )
    {   
        //finds a 1D decomposition
        
        printf("1D\n");
        if( (M % worldSize == 0) )
        {
            //vertical decomposition            
            printf("vertical\n");
            topo->Mp = M/worldSize;
            topo->Np = N;
            dims[0] = 0;
            dims[1] = 1;        //one dimension on y axis
            period[0] = false;
            period[1] = false;
        }
        else if( (N % worldSize == 0) )
        {
            printf("horizontal\n");
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
            printf("Combination of image and thread number is not supported.\n");
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
        printf("2d (%dx%d) \n",dims[0],dims[1]);
        bool isTopologyDivisible = ( (N%dims[0]) == 0) && ( (M % dims[1]) ==0 );
        if(  isTopologyDivisible   )
        {
            topo->Mp=M/dims[0];
            topo->Np=N/dims[1];
        }
        else
        {
            printf("Combination of image and thread number is not supported.\n");
            MPI_Finalize();
            exit(1);
        } 
        
    }
    MPI_Cart_coords(*comm2d, topo->rank, 2, topo->coords);
    
}


double boundaryval(int i, int m)
{
    double val;
    
    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) val = 2.0-val;
    
    return val;
}

void imageRecontruction(topology topo,double** edge,double** buf,double** old,double** new,int N,int* dims)
{
    int i,j,iter;
    int Np = topo.Np;
    int Mp = topo.Mp;
    MPI_Request request,request2,request3,request4;
    MPI_Status status;
    bool isSerialExecution = (dims[0] == 1)&& (dims[1]==1);
    bool isVerticalDecomposition = (dims[1] ==1 );
    
    for (i=1;i<Mp+1;i++)
    {
        for (j=1;j<Np+1;j++)
	{
            edge[i][j]=buf[i-1][j-1];
	}
    }
    
    for (i=0; i<Mp+2;i++)
    {
        for (j=0;j<Np+2;j++)
	{
            old[i][j]=255.0;
	}
    }
    
    /* Set fixed boundary conditions on the left and right sides */
    computeBoundaryConditions(topo,dims,old,N);
    
    
    for (iter=1;iter<=MAXITER; iter++)
    {
        if(iter%PRINTFREQ==0)
	{
            printf("Iteration %d\n", iter);
	}
        
        /* Implement periodic boundary conditions on bottom and top sides */
        
        //vertical
        MPI_Isend(&old[Mp][1],Np,MPI_DOUBLE,topo.right,0,MPI_COMM_WORLD,&request);
        MPI_Isend(&old[1][1],Np,MPI_DOUBLE,topo.left,0,MPI_COMM_WORLD,&request2);
        MPI_Irecv(&old[0][1],Np,MPI_DOUBLE,topo.left,0,MPI_COMM_WORLD,&request3);
        MPI_Irecv(&old[Mp+1][1],Np,MPI_DOUBLE,topo.right,0,MPI_COMM_WORLD,&request4);
      
        MPI_Wait(&request,&status);
        MPI_Wait(&request2,&status);
        MPI_Wait(&request3,&status);
        MPI_Wait(&request4,&status);
        
        
        //horizontal
        halloSwapsHorizontal(old,topo);
        
        
        if(isSerialExecution || isVerticalDecomposition)
        {
            
            for (i=1; i < Mp+1; i++)
            {
                old[i][0]   = old[i][Np];
                old[i][Np+1] = old[i][1];
            }
        }
        for (i=1;i<Mp+1;i++)
	{
            for (j=1;j<Np+1;j++)
	    {
                new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
                        - edge[i][j]);
	    }
	}
	
        for (i=1;i<Mp+1;i++)
	{
            for (j=1;j<Np+1;j++)
	    {
                old[i][j]=new[i][j];
	    }
	}
    }
    
    printf("\nFinished %d iterations\n", iter-1);
    
    for (i=1;i<Mp+1;i++)
    {
        for (j=1;j<Np+1;j++)
	{
            buf[i-1][j-1]=old[i][j];
	}
    }
    
    
}

void halloSwapsHorizontal(double** old,topology topo)
{
    
    int Np = topo.Np;
    int Mp = topo.Mp;
    
    int i;
    double lastRow[Mp];
    double firstRow[Mp];

    
    
    for(i=0; i < Mp; i++)
    {
        lastRow[i] = old[i+1][Np];
    }
    
    for(i=0; i < Mp; i++)
    {
        firstRow[i] = old[i+1][1];
    }
    
    MPI_Request request,request2,request3,request4;
    MPI_Status status;
    MPI_Isend(&firstRow[0],Mp,MPI_DOUBLE,topo.down,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&lastRow[0],Mp,MPI_DOUBLE,topo.up,0,MPI_COMM_WORLD,&request2);
    MPI_Wait(&request,&status);
    MPI_Wait(&request2,&status);
    
    
    MPI_Irecv(&firstRow[0],Mp,MPI_DOUBLE,topo.up,0,MPI_COMM_WORLD,&request4);
    MPI_Irecv(&lastRow[0],Mp,MPI_DOUBLE,topo.down,0,MPI_COMM_WORLD,&request3);
    
    MPI_Wait(&request3,&status);
    MPI_Wait(&request4,&status);
    
    for(i=0; i < Mp; i++)
    {
        old[i+1][Np+1] = firstRow[i];
    }
    
    
    for(i=0; i < Mp; i++)
    {
        old[i+1][0] = lastRow[i];
    }
    
    
}


void haloSwaps()
{
    
}

void computeBoundaryConditions(topology topo,int *dims,double **old,int N)
{
    int j;
    double val;
    int Np = topo.Np;
    int Mp = topo.Mp;
    for (j=1; j < Np+1; j++)
    {
        val = boundaryval(j+(topo.coords[1]*topo.Np ), N);
        int myXCoord = topo.coords[0];
        int lastXCoord = (dims[0]-1);
        
        if(myXCoord == 0)
        {
            old[0][j]   = (int)(255.0*(1.0-val));
        }
        
        if(myXCoord == lastXCoord )
        {
            old[Mp+1][j] = (int)(255.0*val);
        }
    }
}
