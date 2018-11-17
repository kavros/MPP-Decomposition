#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/arralloc.h"
#include "../include/pgmio.h"
#include "../include/argtable3.h"
#include <mpi.h>
#include <stdbool.h> 
#include <assert.h>
#define MAXITER   1500
#define PRINTFREQ  200

//#define IGNORE_PRINTF
#ifdef IGNORE_PRINTF
#define printf(fmt, ...) 
#endif


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
void scatter(double** masterbuf,double** buf,topology topo,int worldSize,MPI_Comm comm2d);
void gather(topology topo,double** masterbuf,double** buf,int M,int N,MPI_Comm comm2d,char* output,int worldSize);
void createDataTypes(topology topo,int N);
void computeBoundaryConditions(topology topo,int *dims,double **old,int N);
void halloSwapsHorizontal(double** old,topology topo);
void halloSwapsVertical(double** old,topology topo);
void printAverages(int iter,double **old,topology topo,int targetIter );
bool isTheLastIteration(topology topo,double maxDelta,int iter);
void loadImage(topology topo,double** masterbuf,char* input,int M,int N);
void allocations(topology topo,int M,int N,double ***masterbuf,double ***buf,double ***old,double ***new, double ***edge);
void deallocations(double **masterbuf,double **buf,double **old,double **new, double **edge);
void cmdLineParser(int argc, char *argv[]);

MPI_Datatype vectorMpxNp;    
MPI_Datatype vectorMpxNP_N;    
MPI_Datatype vectorMpx1;


struct arg_lit *help;
struct arg_int *delta, *targetIter;
struct arg_file *input, *output;
struct arg_end *end;

int main (int argc, char *argv[])
{
    cmdLineParser(argc,argv);

    MPI_Comm comm2d;
    int worldSize;
    topology topo;
    int M,N;
    char *input,*output;
    double **old,**new,**edge,**buf,**masterbuf;
    int dims[2];
    double start,end;
    
    input = "./data/input/edgenew192x128.pgm";
    output="./data/output/imagenew192x128.pgm";
    
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    pgmsize(input,&M,&N);
    
    dims[0]=0;dims[1]=0;    
    initialization(&topo,worldSize,M,N,&comm2d,dims);
    //printf("p%d left=%d,right=%d up=%d down=%d \n",topo.rank,topo.left,topo.right,topo.up,topo.down);
    //printf("p%d Mp=%d, Np=%d\n",topo.rank,topo.Mp,topo.Np);
    //printf("M=%d,N=%d\n",M,N);
    printf("p%d (%d,%d)\n",topo.rank,topo.coords[0],topo.coords[1]);

   
    allocations( topo, M, N,&masterbuf,&buf,&old,&new,&edge);
    
    createDataTypes(topo, N);
    
    loadImage(topo,masterbuf,input, M, N);
    
    start = MPI_Wtime();
    scatter( masterbuf, buf, topo, worldSize, comm2d);
    
    imageRecontruction( topo, edge, buf, old, new,N,dims);
    
    gather( topo, masterbuf,buf, M, N, comm2d,output,worldSize);
    end = MPI_Wtime();
    if(topo.rank == 0)
        fprintf(stdout,"total time is %f(sec), woldSize=%d \n",end-start,worldSize);
    
    MPI_Finalize();
    
    deallocations(masterbuf,old,new,edge,buf);

    return 0;
} 

void cmdLineParser(int argc, char *argv[])
{
    void *argtable[] = {
        help    = arg_litn("h", "help", 0, 1, "display this help and exit"),
        delta   = arg_intn("d","delta","0 or 1",0,1,"activates delta"),
        targetIter = arg_intn("t","target"," must be less than 1500",0,1,"select when to print average value"),
        input   = arg_filen("i", NULL, "<file>", 0, 100, "input file"),
        output  = arg_filen("o", NULL, "<file>", 0, 100, "output file"),
        end     = arg_end(20),
    };
     

    char progname[] = "coursework";
    
    int nerrors;
    nerrors = arg_parse(argc,argv,argtable);
    
    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
    {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        printf("Demonstrate command-line parsing in argtable3.\n\n");
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");

        goto exit;
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout, end, progname);
        printf("Try '%s --help' for more information.\n", progname);

        goto exit;
    }
    if(input->count > 0 )
    {
        printf("input= %s\n",input->filename[0]);
    }
    if (output->count > 0 )
    {
        printf("output= %s\n",output->filename[0]);
    }
    if (targetIter->count > 0 )
    {
        printf("targetIter= %d\n",*(targetIter->ival));
    }
    if(delta->count > 0 )
    {
        printf("delta= %d\n",*(delta->ival));
    }
    
    
exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return;
    
}


void allocations(topology topo,int M,int N,double ***masterbuf,double ***buf,double ***old,double ***new, double ***edge)
{
    *masterbuf = (double**) arralloc(sizeof(double),2,M,N);
    *old = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    *new = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    *edge = (double**) arralloc(sizeof(double),2,topo.Mp+2,topo.Np+2);
    *buf = (double**) arralloc(sizeof(double),2,topo.Mp,topo.Np);
    
}

void deallocations(double **masterbuf,double **buf,double **old,double **new, double **edge)
{
    free(masterbuf);
    free(old);
    free(new);
    free(edge);
    free(buf);
}

void createDataTypes(topology topo,int N)
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

void gather(topology topo,double** masterbuf,double** buf,int M,int N,MPI_Comm comm2d,char* output,int worldSize)
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
            
            MPI_Irecv(&masterbuf[startX][startY],1,vectorMpxNP_N,i,0,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,&status);
        }
        
        for(i=0;i<topo.Mp;i++)
            for(j=0;j<topo.Np;j++)
                masterbuf[i][j] = buf[i][j];
        
        pgmwrite(output, &masterbuf[0][0], M, N);
    }
}


void scatter(double** masterbuf,double** buf,topology topo,int worldSize,MPI_Comm comm2d)
{
    MPI_Request request;
    MPI_Status status;
    if(topo.rank == 0)
    {
        int coords[2];
        int i,j;

        for(i=1; i < worldSize; i++)
        {
            MPI_Cart_coords(comm2d, i, 2, coords);
            int startY = coords[1]*topo.Np;
            int startX = coords[0]*topo.Mp;
            //printf("p %d, startY=%d, startX=%d\n",i,startY,startX);
            
            MPI_Isend(&masterbuf[ startX][startY], 1, vectorMpxNP_N, i, 0, MPI_COMM_WORLD,&request);
            
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
void loadImage(topology topo,double** masterbuf,char* input,int M,int N)
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
    int targetIter=1000;
    double delta,globalMaxDelta,maxDelta;
    delta=globalMaxDelta=maxDelta=-1;
    
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
            //printf("Iteration %d\n", iter);
	}
        
        /* Implement periodic boundary conditions on bottom and top sides */
        //vertical
        halloSwapsVertical(old,topo);
        
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
            delta = abs(new[i][j] - old[i][j]);
            if(maxDelta < delta)
            {
                maxDelta  = delta;
            }

	}
	
        for (i=1;i<Mp+1;i++)
	{
            for (j=1;j<Np+1;j++)
	    {
                old[i][j]=new[i][j];
	    }
	}
        
        //print averages at specified iteration number
        printAverages(iter,old,topo,targetIter );
        
        //terminate based on delta
        /*if( isTheLastIteration(topo,maxDelta,iter) ) 
        {
            //break;
        }*/
            
       
    }
    
    //printf("\nFinished %d iterations\n", iter-1);
    
    for (i=1;i<Mp+1;i++)
    {
        for (j=1;j<Np+1;j++)
	{
            buf[i-1][j-1]=old[i][j];
	}
    }
    
    
}

bool isTheLastIteration(topology topo,double maxDelta,int iter)
{
    double globalMaxDelta = -1.0f;
    MPI_Allreduce(&maxDelta,&globalMaxDelta,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    /*if(topo.rank==0)
    {
        printf("p%d  maxDelta=%f globalMaxDelta=%f \n",topo.rank,maxDelta,globalMaxDelta);
    }*/
    
    if(globalMaxDelta < 0.1 )
    {
        printf(" my_rank = %d, globalMaxDelta= %f, iterations = %d \n",topo.rank,globalMaxDelta,iter );
        return true;
    }
    maxDelta =-10;
    return false;
}

void printAverages(int iter,double **old,topology topo,int targetIter )
{
    int i,j;
    if(targetIter == iter)
        {
            double sum = 0;
            double totalSum;
            for (i=1;i<topo.Mp+1;i++)
            {
                for (j=1;j<topo.Np+1;j++)
                {
                    sum=+old[i][j];
                }
            }
            //printf("p%d sum=%f\n",topo.rank,sum);
            MPI_Reduce(&sum,&totalSum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            if(topo.rank == 0)
            {
                int worldSize;
                MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
                printf("average value of pixels is =%f\n",(totalSum/(double)worldSize));
            }
        }
}
 
void halloSwapsVertical(double** old,topology topo)
{
    MPI_Request request,request2,request3,request4;
    MPI_Status status;
    int Np = topo.Np;
    int Mp = topo.Mp;
    
    MPI_Isend(&old[topo.Mp][1],Np,MPI_DOUBLE,topo.right,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&old[1][1],Np,MPI_DOUBLE,topo.left,0,MPI_COMM_WORLD,&request2);
    MPI_Irecv(&old[0][1],Np,MPI_DOUBLE,topo.left,0,MPI_COMM_WORLD,&request3);
    MPI_Irecv(&old[Mp+1][1],Np,MPI_DOUBLE,topo.right,0,MPI_COMM_WORLD,&request4);
    
    MPI_Wait(&request,&status);
    MPI_Wait(&request2,&status);
    MPI_Wait(&request3,&status);
    MPI_Wait(&request4,&status);
}

void halloSwapsHorizontal(double** old,topology topo)
{
    
    int Np = topo.Np;
   
    MPI_Request request,request2,request3,request4;
    MPI_Status status;
    
    MPI_Isend(&old[1][1],1,vectorMpx1,topo.down,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&old[1][Np],1,vectorMpx1,topo.up,0,MPI_COMM_WORLD,&request2);
    MPI_Irecv(&old[1][Np+1],1,vectorMpx1,topo.up,0,MPI_COMM_WORLD,&request4);
    MPI_Irecv(&old[1][0],1,vectorMpx1,topo.down,0,MPI_COMM_WORLD,&request3);
    
    MPI_Wait(&request,&status);
    MPI_Wait(&request2,&status);    
    MPI_Wait(&request3,&status);
    MPI_Wait(&request4,&status);
    
    
    
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
