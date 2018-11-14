/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "../include/pgmio.h"
#include<assert.h>
#include "../include/arralloc.h"
//#define M 192
//#define N 128
#define MAXITER   1500
#define PRINTFREQ  200

#define HORIZONTAL 1
#define VERTICAL 2
#define TWO_DIM      3
#define SERIAL  4
typedef struct decompositionType
{
    int type;
}decompositionType;

double boundaryval(int i, int m);


void validation(double** expectedImage,double** masterbuf);
void horizontalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf);
void horizontalComposition(double** masterbuf2,double** buf1,double** buf,int worldSize,int Mp,int Np);
void verticalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf);
void verticalComposition(double** masterbuf2,double** buf,int worldSize,int Mp,int Np);
void composition(double** buf,double** masterbuf, int worldRank,decompositionType decomp,int worldSize,int Mp,int Np);
void decomposition(decompositionType decomp,int worldRank,int Np,int Mp,double** buf,double** masterbuf);
void twoDimensionalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf);
void twoDimensionalComposition(int worldSize,double** buf,double** masterbuf2,int Np,int Mp);
void imageRecontruction(decompositionType decomp, int worldRank,int worldSize,int Mp,int Np,double** edge,double** buf,double** old,double** new);
void computeSawtooth(decompositionType decomp,double** old,int Np,int Mp);
void computeSawtoothVertical(double** old,int Np,int Mp);
void initDecompositionType(decompositionType* decomp,int* Np,int* Mp);
void retriveImage(char* filename,double** masterbuf,double** buf,int worldRank,int worldSize);
void computeSawtoothSerial(double** old, int Np,int Mp);
void computeSawtoothTwoDim(double** old,int Np,int Mp);
void computeSawtoothHorizontal(double** old,int Np,int Mp);

void halloSwapsVertical(double** old,int Np,int Mp,int worldRank,int worldSize);
void halloSwaps(decompositionType decomp,double** old,int Np,int Mp,int worldRank,int worldSize);
void halloSwapsSerial(double** old,int Np,int Mp,int worldRank,int worldSize);
void halloSwapsHorizontal(double** old,int Np,int Mp,int worldRank,int worldSize);
int M,N;
int main (void)
{
    
    MPI_Status status;
    int worldSize,worldRank;
    char *filename;
    
    //initialization of width and height
    filename = "./images/input/edgenew192x128.pgm";
    pgmsize(filename,&M,&N);
    printf("M=%d,N=%d\n",M,N);
    
    //initialization of MPI variables
    MPI_Request request;
    MPI_Init(NULL,NULL);        
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    
    //initialization of Mp,Np
    int Np,Mp;
    decompositionType decomp;
    initDecompositionType(&decomp,&Np,&Mp);
    printf("Mp=%d, Np=%d\n",Mp,Np);
    
    //memory allocations
    double** masterbuf,**old,**new,**edge,**buf;    
    masterbuf = (double**) arralloc(sizeof(double),2,M,N);    
    old = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    new = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    edge = (double**) arralloc(sizeof(double),2,Mp+2,Np+2);
    buf = (double**) arralloc(sizeof(double),2,Mp,Np);
    

    retriveImage(filename,masterbuf,buf,worldRank,worldSize);
    
    decomposition( decomp, worldRank, Np, Mp, buf, masterbuf);
    
    imageRecontruction(decomp,worldRank,worldSize,Mp,Np,edge,buf,old,new);
        
    if(worldRank !=0)
    {
        //printf("--->p%d buf[0][0]=%f\n",worldRank,buf[0][0]);
        MPI_Isend(&buf[0][0],Mp*Np,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
    }
    
    composition(buf,masterbuf,worldRank,decomp,worldSize,Mp,Np);
    
    
    MPI_Finalize();
    free(masterbuf);
    free(buf);
    free(old);
    free(edge);
    free(new);
    return 0;
} 

void retriveImage(char* filename,double** masterbuf,double** buf,int worldRank,int worldSize)
{
    if(worldRank == 0)
    {
        
        printf("Processing %d x %d image\n", M, N);
        printf("Number of iterations = %d\n", MAXITER);
        
        printf("\nReading <%s>\n", filename);
        if(worldSize == 1)
        {
            pgmread(filename, &buf[0][0], M, N);
        }
        else
        {
            pgmread(filename, &masterbuf[0][0], M, N);
        }
        printf("\n");
    }    
    MPI_Bcast(&masterbuf[0][0],(M*N),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void initDecompositionType(decompositionType* decomp,int* Np,int* Mp)
{
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    
    
    if( (worldSize == 4) && ((N%2)==0)  && ((M%2)==0)  )
    {
        decomp->type=TWO_DIM;
        *Np = N/2;
        *Mp = M/2;
        printf("2D\n");
    }
    else if(worldSize == 1)
    {
        decomp->type =SERIAL;
        *Np=N;
        *Mp =M;
        printf("Serial\n");
    }
    else if( N%worldSize == 0)
    {
        *Mp = M;
        *Np = N/worldSize;
        decomp->type=HORIZONTAL;
        printf("horizontal\n");
    }
    else if(M%worldSize == 0)
    {
        *Np = N;
        *Mp = M/worldSize;
        decomp->type= VERTICAL;
        printf("vertical\n");
    }

    

}

void computeSawtooth(decompositionType decomp,double** old,int Np,int Mp)
{

    if(decomp.type == VERTICAL)
    {
        computeSawtoothVertical(old,Np,Mp);
    }
    else if(decomp.type == SERIAL)
    {
        computeSawtoothSerial(old,Np,Mp);
    }
    else if(decomp.type == HORIZONTAL)
    {
        computeSawtoothHorizontal(old,Np,Mp);
    }
    else if(decomp.type == TWO_DIM)
    {
        computeSawtoothTwoDim(old,Np,Mp);
    }
}

void computeSawtoothHorizontal(double** old,int Np,int Mp)
{
    double val;
    int j,worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    for (j=1; j < Np+1; j++)
    {
        // compute sawtooth value    
        val = boundaryval(j+(worldRank*Np), N);
        
        old[0][j]   = (int)(255.0*(1.0-val));
        old[Mp+1][j] = (int)(255.0*val);
    }
}

void computeSawtoothTwoDim(double** old,int Np,int Mp)
{
    
}


void computeSawtoothVertical(double** old,int Np,int Mp)
{
    
    int worldRank,worldSize,j;
    double val;
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    if(worldRank == 0 )
    {
        for (j=1; j < Np+1; j++)
        {
            // compute sawtooth value    
            val = boundaryval(j, Np);
            
            old[0][j]   = (int)(255.0*(1.0-val));
        }
    }
    if(worldRank == worldSize-1 )
    {
        for (j=1; j < Np+1; j++)
        {
            // compute sawtooth value    
            val = boundaryval(j, Np);
            
            
            old[Mp+1][j] = (int)(255.0*val);
        }
    }
}

void computeSawtoothSerial(double** old, int Np,int Mp)
{
    int j;
    double val;
    for (j=1; j < Np+1; j++)
    {
        // compute sawtooth value    
        val = boundaryval(j, Np);
        
        old[0][j]   = (int)(255.0*(1.0-val));
        old[Mp+1][j] = (int)(255.0*val);
    }
}



void imageRecontruction(decompositionType decomp, int worldRank,int worldSize,int Mp,int Np,double** edge,double** buf,double** old,double** new)
{
    int i,j,iter;

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
 
    
    computeSawtooth(decomp,old,Np,Mp);
    
    for (iter=1;iter<=MAXITER; iter++)
    {
        if(iter%PRINTFREQ==0)
        {
            //printf("Iteration %d\n", iter);
        }
        
        // Implement periodic boundary conditions on bottom and top sides 
        halloSwaps(decomp,old,Np,Mp,worldRank,worldSize);
      
        //if(decomp.type== SERIAL)
        //{
            
        //}
        
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

void halloSwaps(decompositionType decomp,double** old,int Np,int Mp,int worldRank,int worldSize)
{
    if(decomp.type == VERTICAL)
    {
        halloSwapsVertical(old,Np,Mp,worldRank,worldSize);
    }
    else if(decomp.type == HORIZONTAL)
    {
        halloSwapsHorizontal(old, Np,Mp, worldRank, worldSize);
    }
    else if(decomp.type == TWO_DIM)
    {
        
    }
    else if(decomp.type==SERIAL)
    {
        halloSwapsSerial(old,Np,Mp,worldRank,worldSize);
    }
}

void halloSwapsSerial(double** old,int Np,int Mp,int worldRank,int worldSize)
{
    int i;
    for (i=1; i < Mp+1; i++)
    {
        old[i][0]   = old[i][Np];
        old[i][Np+1] = old[i][1];
    }
}

void halloSwapsVertical(double** old,int Np,int Mp,int worldRank,int worldSize)
{
    
    MPI_Request request,request2,request3,request4;
    MPI_Status status;
    
    
    int right = worldRank+1;
    int left  = worldRank-1;
    if(right == worldSize) right=MPI_PROC_NULL;
    if(left == -1) left =MPI_PROC_NULL;
    
    MPI_Isend(&old[Mp][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&old[1][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request2);

    MPI_Wait(&request,&status);
    MPI_Wait(&request2,&status);

    MPI_Irecv(&old[0][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request3);
    MPI_Irecv(&old[Mp+1][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request4);

    MPI_Wait(&request3,&status);
    MPI_Wait(&request4,&status);
    halloSwapsSerial(old,Np,Mp,worldRank,worldSize);
}

void halloSwapsHorizontal(double** old,int Np,int Mp,int worldRank,int worldSize)
{
    
    int i;
    double lastRow[Mp];
    double firstRow[Mp];

    
    int down = worldRank-1;
    int up = worldRank+1;
    
    if(up == worldSize) up = 0;
    if(down == -1) down=worldSize-1;
    
    //printf("p%d up=%d down=%d \n",worldRank,up,down);
    
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
    MPI_Isend(&firstRow[0],Mp,MPI_DOUBLE,down,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&lastRow[0],Mp,MPI_DOUBLE,up,0,MPI_COMM_WORLD,&request2);
    MPI_Wait(&request,&status);
    MPI_Wait(&request2,&status);
    
    
    MPI_Irecv(&firstRow[0],Mp,MPI_DOUBLE,up,0,MPI_COMM_WORLD,&request4);
    MPI_Irecv(&lastRow[0],Mp,MPI_DOUBLE,down,0,MPI_COMM_WORLD,&request3);
    
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

void decomposition(decompositionType decomp,int worldRank,int Np,int Mp,double** buf,double** masterbuf)
{
    if(decomp.type==HORIZONTAL)
    {
        horizontalDecomposition(worldRank,Np,Mp,buf,masterbuf);
    }
    else if(decomp.type == VERTICAL)
    {
        verticalDecomposition(worldRank,Np,Mp,buf,masterbuf);
    }    
    else if(decomp.type == TWO_DIM)
    {
        twoDimensionalDecomposition( worldRank, Np, Mp, buf,masterbuf);
    }
    else if(decomp.type == SERIAL)
    {
        //nothing to do here
        
    }
}



void twoDimensionalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf)
{
    int r,c;
    for( r =0; r <Mp; r++)
    {    
        for( c=0; c < Np; c++)
        {       
            buf[r][c] = masterbuf[r+(Mp*(worldRank/2))][c+(Np*(worldRank%2))];   
        } 
        
    }
}

void twoDimensionalComposition(int worldSize,double** buf,double** masterbuf2,int Np,int Mp)
{
    int r,c;
    double buf1[Mp][Np];
    double buf2[Mp][Np];
    double buf3[Mp][Np];
    MPI_Request request01;
    MPI_Request request02;
    MPI_Request request03;
    MPI_Status status;
    for( r = 1 ; r < worldSize; r++)
    {
        //int col = (r/2)*Mp;
        //int row = (r%2)*Np;
        if(  r ==1)
        {
            MPI_Irecv(&buf1[0][0],Mp*Np,MPI_DOUBLE,r,0,MPI_COMM_WORLD,&request01);
        }
        else if( r==2)
        {
            MPI_Irecv(&buf2[0][0],Mp*Np,MPI_DOUBLE,r,0,MPI_COMM_WORLD,&request02);
        }
        else if(r==3)
        {
            MPI_Irecv(&buf3[0][0],Mp*Np,MPI_DOUBLE,r,0,MPI_COMM_WORLD,&request03);
        }
        
        
        /*    printf("---p%d (%d,%d)\n",i,col,row);
         //printf("p%d masterbuf2[0][0] = %f \n",i,masterbuf2[0][0]);
         
         MPI_Irecv(&masterbuf2[col][row],Mp*Np,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request2);
         MPI_Wait(&request2,&status);
         */   
        //printf("p%d masterbuf2 = %f \n",i,masterbuf2[col][row]);
        
    }
    MPI_Wait(&request01,&status);
    MPI_Wait(&request02,&status);
    MPI_Wait(&request03,&status);
    
    //int r,c;
    for( r=0; r <M;r++ )
    {
        for( c=0; c <N;c++ )
        {
            if(r >=0 && r <M/2  )
            {
                if(c<N/2)
                {
                    masterbuf2[r][c] = buf[r][c];
                }
                else
                {
                    masterbuf2[r][c] = buf1[r][c-Np];
                }
            }
            else if(r >= M/2)
            {
                if(c<N/2)
                {
                    
                    assert((r-Mp)>= 0 );
                    assert((r-Mp)< Mp  );
                    
                    masterbuf2[r][c] = buf2[r-Mp][c];
                }
                else
                {
                    assert((r-Mp)>= 0 );
                    assert((r-Mp)< Mp  );
                    assert((c-Np) >= 0);
                    assert((c-Np) < Np);
                    
                    
                    masterbuf2[r][c] = buf3[r-Mp][c-Np];
                }
            }
            
        }
    }
}


void composition(double** buf,double** masterbuf,int worldRank,decompositionType decomp,int worldSize,int Mp,int Np)
{
    if(worldRank == 0)
    {

        double** buf1=NULL;
        double** masterbuf2=NULL;

        if(decomp.type != SERIAL)
        {
            buf1 = (double**) arralloc(sizeof(double),2,Mp,Np);
            masterbuf2 = (double**) arralloc(sizeof(double),2,M,N);        
        }
        
        
        if(decomp.type == HORIZONTAL)
        {
            horizontalComposition( masterbuf2, buf1, buf, worldSize,Mp, Np);
        }
        else if(decomp.type == VERTICAL)
        {
            verticalComposition(masterbuf2,buf,worldSize,Mp,Np);
        }
        else if(decomp.type == TWO_DIM )
        {
            twoDimensionalComposition(worldSize,buf,masterbuf2,Np,Mp);
        }else if( decomp.type == SERIAL )
        {
            printf("serial composition\n");
            masterbuf2 = buf;
        }
        
        //validation(masterbuf,masterbuf2);
        
        char* filename="./images/output/imagenew192x128.pgm";
        printf("\nWriting <%s>\n", filename); 
        pgmwrite(filename, &masterbuf2[0][0], M, N);
        
        if(decomp.type != SERIAL)
        {
            free(masterbuf2); 
            free(buf1);
        }
    }
   
}

void horizontalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf)
{
    int i,j;
    for(i=0; i < Mp;i++)
    {
        for(j=0; j < Np; j++)
        {
            if(worldRank == 0)
            {
                buf[i][j] = masterbuf[i][j];
            }
            else
            {
                buf[i][j] = masterbuf[i][j+(Np*worldRank)];
            }
        }
    }
}

void verticalDecomposition(int worldRank,int Np,int Mp,double** buf,double** masterbuf)
{
    int i,j;
    for( i =0; i <Mp; i++)
    {
        for( j=0; j < Np; j++)
        {
            buf[i][j] = masterbuf[i+(Mp*worldRank)][j];
            
        }       
    }

}

void horizontalComposition(double** masterbuf2,double** buf1,double** buf,int worldSize,int Mp,int Np)
{
    int i,j;
    MPI_Request request;
    MPI_Status status;
    

    int p=0;
    for(p =1; p < worldSize; p++)
    {
        MPI_Irecv(&buf1[0][0],Mp*Np,MPI_DOUBLE,p,0,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);

        for(i =0; i < Mp; i++)
        {
            for(j=0; j < Np; j++)
            {
                masterbuf2[i][j+(p*Np)] = buf1[i][j];
            }
        }
    }

    for(i =0; i < Mp; i++)
    {
        for(j=0; j < Np; j++)
        { 
            masterbuf2[i][j] = buf[i][j];
        }

    }
    
}

void verticalComposition(double** masterbuf2,double** buf,int worldSize,int Mp,int Np)
{
    int i,j;
    MPI_Request request;
    MPI_Status status;
    for( i = 1 ; i < worldSize; i++)
    { 
        MPI_Irecv(&masterbuf2[i*Mp][0],Mp*Np,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
    }
    
    for( i=0; i <Mp;i++ )
    {
        for( j=0; j <Np;j++ )
        {
            masterbuf2[i][j] = buf[i][j];
        }
    }
}


void validation(double** expectedImage,double** masterbuf)
{

    int i,j,cnt=0;
    
    for(i=0;i<M;i++)
    {
        for(j=0; j < N; j++)
        {
            if(expectedImage[i][j] != masterbuf[i][j])
            {
                cnt++;
            }
        }
    }
    int totalPixels  = N*M;
    double percentage = ((double)cnt)/((double)totalPixels);
    printf("difference from expected image = %f",percentage);
    printf("cnt = %d\n",cnt);
    
}

double boundaryval(int i, int m)
{
    double val;
    
    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) val = 2.0-val;
    
    return val;
}

void communication()
{
      //vertical
        /*
         * 
         *     int right = worldRank+1;
    int left  = worldRank-1;
    if(right == worldSize) right=MPI_PROC_NULL;
    if(left == -1) left =MPI_PROC_NULL;
         * 
         * 
         MPI_Isend(&old[Mp][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request);
        MPI_Isend(&old[1][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request2);
        
        MPI_Wait(&request,&status);
        MPI_Wait(&request2,&status);
        
        MPI_Irecv(&old[0][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request3);
        MPI_Irecv(&old[Mp+1][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request4);
        
        MPI_Wait(&request3,&status);
        MPI_Wait(&request4,&status);
       */
        /*
        for(k=1; k < Mp; k++)
        {
            send1[k-1]=old[k][1];  //first row
        }
        
        for(k=1; k < Mp; k++)
        {
            send2[k-1]=old[k][Np];  //last row
        }

        
     
        //horizontal
        MPI_Isend(&send1[0],Mp,MPI_DOUBLE,up,0,MPI_COMM_WORLD,&request);
        MPI_Isend(&send2[0],Mp,MPI_DOUBLE,down,0,MPI_COMM_WORLD,&request2);
        MPI_Wait(&request,&status);
        MPI_Wait(&request2,&status);
        
        
        MPI_Irecv(&old[1][Np+1],Mp,MPI_DOUBLE,up,0,MPI_COMM_WORLD,&request4);
        MPI_Irecv(&old[1][0],Mp,MPI_DOUBLE,down,0,MPI_COMM_WORLD,&request3);
        
        MPI_Wait(&request3,&status);
        MPI_Wait(&request4,&status);
       
        */
}


    /*if(worldRank == 0)
    {
        double masterbuf2[M][N];
        double buf1[Mp][Np];
        double buf2[Mp][Np];
        double buf3[Mp][Np];
        MPI_Request request01;
        MPI_Request request02;
        MPI_Request request03;
        
        
        for( i = 1 ; i < worldSize; i++)
        {
            //int col = (r/2)*Mp;
            //int row = (r%2)*Np;
            if(  i ==1)
            {
                MPI_Irecv(&buf1[0][0],Mp*Np,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request01);
            }
            else if( i==2)
            {
                MPI_Irecv(&buf2[0][0],Mp*Np,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request02);
            }
            else if(i==3)
            {
                MPI_Irecv(&buf3[0][0],Mp*Np,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request03);
            }
            
        }
        MPI_Wait(&request01,&status);
        MPI_Wait(&request02,&status);
        MPI_Wait(&request03,&status);

        //int r,c;
        for( i=0; i <M;i++ )
        {
            for( j=0; j <N;j++ )
            {
                if(i >=0 && i <M/2  )
                {
                    if(j<N/2)
                    {
                        masterbuf2[i][j] = buf[i][j];
                    }
                    else
                    {
                        masterbuf2[i][j] = buf1[i][j-Np];
                    }
                }
                else if(i >= M/2)
                {
                    if(j<N/2)
                    {
                        
                        assert((i-Mp)>= 0 );
                        assert((i-Mp)< Mp  );
                        
                        masterbuf2[i][j] = buf2[i-Mp][j];
                    }
                    else
                    {
                        assert((i-Mp)>= 0 );
                        assert((i-Mp)< Mp  );
                        assert((j-Np) >= 0);
                        assert((j-Np) < Np);
                        
                        
                        masterbuf2[i][j] = buf3[i-Mp][j-Np];
                    }
                }
                
            }
        }
        
//        int cnt=0;
//        for(i =0; i < M; i++){
//            for(j=0; j< N; j++)
//            {
//                if((masterbuf2[i][j] != masterbuf[i][j]))
//                {
//                    cnt++;  
//                }
//            }
//        }
//        printf("--cnt=%d\n",cnt);
        
        
        //validation("imagenew192x128.pgm",masterbuf2);
        
        filename="imagenew192x128.pgm";
        printf("\nWriting <%s>\n", filename); 
        pgmwrite(filename, masterbuf2, M, N);
     
    }   */


 
    
    /*if(worldRank == 0)
    {
        pgmwrite("x0.pgm", &buf[0][0], Mp, Np);
    }
    if(worldRank == 1)
    {
        pgmwrite("x1.pgm", &buf[0][0], Mp, Np);
    }
    
   */
    
    /*int right,left,up,down;
   
    if(worldRank == 0)
    {
        right = 2;
        up = 1;
        down=1;
        left=MPI_PROC_NULL;
    }
    
    if(worldRank==1)
    {
        right = 3;
        up = 0;
        down=0;
        left=MPI_PROC_NULL;
    }
    
    
    if(worldRank==2)
    {
        right = MPI_PROC_NULL;
        up = 3;
        down=3;
        left=0;
    }
    if(worldRank==3)
    {
        right = MPI_PROC_NULL;
        up = 2;
        down=2;
        left=1;
    }
    //printf("p%d left=%d right=%d up = %d, down=%d \n",worldRank,left,right,up,down);
    */
