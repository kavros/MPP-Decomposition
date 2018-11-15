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
}topology;

double boundaryval(int i, int m);
void imageRecontruction(int Mp,int Np,double** edge,double** buf,double** old,double** new);
bool isNumberPrime(int num);
void initialization(topology* topo,int worldSize,int M,int N);

int main (void)
{
    
    int worldSize;
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    
    int M,N;
    char *filename;
    filename = "./images/input/edgenew192x128.pgm";
    pgmsize(filename,&M,&N);
    
    topology topo;
    
    int Mp =M;
    int Np=N;
    initialization(&topo,worldSize,M,N);
    printf("p%d left=%d,right=%d up=%d down=%d \n",topo.rank,topo.left,topo.right,topo.up,topo.down);
    printf("p%d Mp=%d, Np=%d\n",topo.rank,topo.Mp,topo.Np);
    printf("M=%d,N=%d\n",M,N);
    
    
    double **old,**new,**edge;
    double **buf,**masterbuf;
    masterbuf = (double**) arralloc(sizeof(double),2,M,N);
    old = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    new = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    edge = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    buf = (double**) arralloc(sizeof(double),2,Mp,Np);
    
    
   /* 
    printf("Processing %d x %d image\n", M, N);
    printf("Number of iterations = %d\n", MAXITER);
    
    printf("\nReading <%s>\n", filename);
    pgmread(filename, &buf[0][0], M, N);
    printf("\n");
    
    
    imageRecontruction(Mp,Np,edge,buf,old,new);
    
    int i;
    for(i=1; i < 50;i++)
    {
        if(isNumberPrime(i))
        {
            printf("%d,",i);
        }
    }
    
    filename="./images/output/imagenew192x128.pgm";
    printf("\nWriting <%s>\n", filename); 
    pgmwrite(filename, &buf[0][0], M, N);
    */
    MPI_Finalize();
    
    free(masterbuf);
    free(old);
    free(new);
    free(edge);
    free(buf);
    return 0;
    
    
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

void initialization(topology* topo,int worldSize,int M,int N)
{
    
    int disp;
    const int ndims =2;
    int dims[ndims];
    int period[ndims];
    int reorder;
    MPI_Comm comm;
    MPI_Comm comm2d;
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
        MPI_Cart_create(comm,ndims,dims,period,reorder,&comm2d);
        MPI_Cart_shift(comm2d,0,disp,&(topo->left),&(topo->right));
        MPI_Cart_shift(comm2d,1,disp,&(topo->down),&(topo->up));
        MPI_Comm_rank(comm2d,&(topo->rank));
        
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
        MPI_Cart_create(comm,ndims,dims,period,reorder,&comm2d);
        MPI_Cart_shift(comm2d,0,disp,&(topo->left),&(topo->right));
        MPI_Cart_shift(comm2d,1,disp,&(topo->down),&(topo->up));
        MPI_Comm_rank(comm2d,&(topo->rank));
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
}


double boundaryval(int i, int m)
{
    double val;
    
    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) val = 2.0-val;
    
    return val;
}

void imageRecontruction(int Mp,int Np,double** edge,double** buf,double** old,double** new)
{
    int i,j,iter;
    double val;
    
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
    
    for (j=1; j < Np+1; j++)
    {
        
        val = boundaryval(j, Np);
        
        old[0][j]   = (int)(255.0*(1.0-val));
        old[Mp+1][j] = (int)(255.0*val);
    }
    
    for (iter=1;iter<=MAXITER; iter++)
    {
        if(iter%PRINTFREQ==0)
	{
            printf("Iteration %d\n", iter);
	}
        
        /* Implement periodic boundary conditions on bottom and top sides */
        
        for (i=1; i < Mp+1; i++)
	{
            old[i][0]   = old[i][Np];
            old[i][Np+1] = old[i][1];
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