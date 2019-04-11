//B145772
#include <mpi.h>
#include <stdbool.h> 
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/arralloc.h"
#include "../include/pgmio.h"
#include "../include/argtable3.h"
#include "../include/functions.h"
#include "../include/initializations.h"
#include "../include/communications.h"
//#define IGNORE_PRINTF
#ifdef IGNORE_PRINTF
#define printf(fmt, ...) 
#endif
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

void parseCmdLine(int argc, char *argv[])
{
    struct arg_lit *help;
    struct arg_int *deltaArg, *totalAvgPrintsArg;
    struct arg_file *inputArg, *outputArg;
    struct arg_end *end;

   void* argtable[] = {
        help    = arg_litn("h", "help", 0, 1, "display this help and exit"),
        deltaArg   = arg_intn("d","delta"," please set -d 1 to activate delta termination",0,1,"activates delta"),
        totalAvgPrintsArg = arg_intn("t","target"," must be positive number less than 1500 ",0,1,"select when to print average value"),
        inputArg   = arg_filen("i", NULL, "<file>", 0, 100, "input file"),
        outputArg  = arg_filen("e", NULL, "<file>", 0, 100, "output file"),
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

        
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        exit(-1);
        
    }
    if(inputArg->count > 0 )
    {
        
        input = (char*)(inputArg->filename[0]);
        
    }
    if (outputArg->count > 0 )
    {
        
        output =(char*) (outputArg->filename[0]);
        
       //printf("output= %s\n",outputArg->filename[0]);
    }
    if (totalAvgPrintsArg->count > 0 )
    {
        totalAveragePrints = *(totalAvgPrintsArg->ival);
        if(totalAveragePrints <= 0 )
        {
            totalAveragePrints =1;
            printf("Negative or zero values are not accepted for the total prints of averages.\n");
            printf("total prints of averages is now 1\n");
        }
        //printf("totalAvgPrintsArg= %d\n",totalAveragePrints);
    }
    if(deltaArg->count > 0 )
    {
        isDeltaActivated = *(deltaArg->ival);
        //printf("delta= %d\n",*(deltaArg->ival));
    }
    
}


void allocations(topology topo,double ***masterbuf,double ***buf,double ***old,double ***new, double ***edge)
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

bool isNumberPrime(int num)
{
    assert(num > 0);       //validates that the program does not give negative number as input.
    int i;
    if(num == 1)            //we consider 1 as prime because if the worldsize is 1 then
                            //the decomposition is 1D.
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

double boundaryval(int i, int m)
{
    double val;
    
    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) val = 2.0-val;
    
    return val;
}

void imageRecontruction(topology topo,double** edge,double** buf,double** old,double** new,int* dims)
{
    int i,j,iter;
    int Np = topo.Np;
    int Mp = topo.Mp;
    int cntAvgPrints=0;
    int printAvgAtIter[totalAveragePrints];
    initAvgPrints( printAvgAtIter);
    

    
    double maxDelta=-1;     //set maxDelta to -1 in order to take the value of delta at the first iteration 
    
    
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
    computeBoundaryConditions(topo,dims,old);
    
    
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

        #pragma omp parallel for default(none) shared(new,old,maxDelta,Mp,Np,edge) private(i,j)
        for (i=1;i<Mp+1;i++)
        {
            for (j=1;j<Np+1;j++)
            {
                new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]
                        - edge[i][j]);
                calculateMaxDelta(i,j,new,old,&maxDelta);
            }
        }
    
        
        #pragma omp for schedule(static)
        for (i=1;i<Mp+1;i++)
	    {
            for (j=1;j<Np+1;j++)
            {
                old[i][j]=new[i][j];
            }
	    }
        
        //print averages at specified iteration number
        printAverages((iter-1),printAvgAtIter[cntAvgPrints],&cntAvgPrints,old,topo);
        
        //terminate based on delta
        if( isTheLastIteration(topo,maxDelta,iter) ) 
        {
            break;
        }
            
        maxDelta = -1; //set maxDelta to -1 in order to take the value of delta at the first iteration 
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

void calculateMaxDelta(int i,int j,double** new,double** old,double* maxDelta)
{   
    if(isDeltaActivated)
    {
        
        double delta = fabs( (new[i][j] - old[i][j]) );
        if( (*maxDelta) < delta)
        {
            ///printf("delta = %f\n",*maxDelta);
            *maxDelta  = delta;
        }
    }
}

bool isTheLastIteration(topology topo,double maxDelta,int iter)
{
    if (isDeltaActivated == false) //continue iterations if delta is disabled
    {
        return false;
    }
    
    if(iter%DELTAFREQ != 0 ) return false;  //continue iterations 
                                            //if current iteration 
                                            //is not devisible by DELTAFREQ
    
    //printf("iter = %d\n",iter);
    double globalMaxDelta = -1.0f;
    MPI_Allreduce(&maxDelta,&globalMaxDelta,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    if(globalMaxDelta < 0.1  )
    {
        if(topo.rank == 0)
            printf("\n globalMaxDelta = %f at iterations = %d \n",globalMaxDelta,iter );
        return true;
    }
    return false;
    
}

void printAverages(int iter,int targetIter,int* cnt ,double **old,topology topo )
{
    
    if( (*cnt) == totalAveragePrints ) return;      //if we have print all the
                                                    //averages we return.
    
    //printf("targetIter = %d\n",targetIter);
    int i,j;
    if(targetIter == iter)
    {
        *cnt=*cnt+1;
        double sum = 0;
        double totalSum=0;
        for (i=1;i<topo.Mp+1;i++)
        {
            for (j=1;j<topo.Np+1;j++)
            {
                sum=sum+old[i][j];
            }
        }
        //printf("p%d sum=%f\n",topo.rank,sum);
        MPI_Reduce(&sum,&totalSum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        if(topo.rank == 0)
        {
            int totalPixels =N*M;
            printf("avg number of pixels=%f, iteration=%d \n",(totalSum/(double)totalPixels),targetIter);
        }
    }
}
 
void computeBoundaryConditions(topology topo,int *dims,double **old)
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
        
        if(myXCoord == 0)                           //compute boundary condition
                                                    //at the most left column
        {
            old[0][j]   = (int)(255.0*(1.0-val));
        }
        
        if(myXCoord == lastXCoord )                 //compute boundary condition
                                                    //conditions at the most right column
        {
            old[Mp+1][j] = (int)(255.0*val);
        }
    }
}

