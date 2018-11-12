/*
 * A simple serial solution to the Case Study exercise from the MPP
 * course.  Note that this uses the alternative boundary conditions
 * that are appropriate for the assessed coursework.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "pgmio.h"
#include<assert.h>
#define M 192
#define N 128

#define MAXITER   1500
#define PRINTFREQ  200

double boundaryval(int i, int m);


int main (void)
{
    
    MPI_Status status;
    int worldSize,worldRank;
    
    MPI_Request request;
    MPI_Request request2;    
    
    MPI_Init(NULL,NULL);
    
    
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    
    int Np = N;
    int Mp = M/worldSize;
    
    if(worldSize%2 == 0)
    {
        Np = N/2;
        Mp = (M/2);
        printf("Np = %d, Mp=%d\n",Np,Mp);
    }
    
    double masterbuf[M][N];
    double old[Mp+2][Np+2], new[Mp+2][Np+2], edge[Mp+2][Np+2];   
    double buf[Mp][Np];
    
    int r, c, iter, maxiter;
    char *filename;
    double val;
    
    if(worldRank == 0)
    {
        //printf("master process\n");
        
        printf("Processing %d x %d image\n", M, N);
        printf("Number of iterations = %d\n", MAXITER);
        
        filename = "edgenew192x128.pgm";
        
        printf("\nReading <%s>\n", filename);
        pgmread(filename, &masterbuf, M, N);
        printf("\n");
        
    }
    
    //printf("before p%d buf[0][0] = %f\n",worldRank,masterbuf[0][0]);
    MPI_Bcast(&masterbuf,(M*N),MPI_DOUBLE,0,MPI_COMM_WORLD);
    //printf("after p%d buf[0][0] = %f\n",worldRank,masterbuf[0][0]);
    int x =0;
    printf("++++p%d,(%d,%d) ",worldRank,(Mp*(worldRank/2)),Np*(worldRank%2));
    for( r =0; r <Mp; r++)
    {    
        for( c=0; c < Np; c++)
        {       
            buf[r][c] = masterbuf[r+(Mp*(worldRank/2))][c+(Np*(worldRank%2))];   
        } 
        
    }
    
    if(worldRank == 0)
    {
        pgmwrite("x0.pgm", buf, Mp, Np);
    }
    if(worldRank == 1)
    {
        pgmwrite("x1.pgm", buf, Mp, Np);
    }
    if(worldRank == 2)
    {
        pgmwrite("x2.pgm", buf, Mp, Np);
    }
    if(worldRank == 3)
    {
        pgmwrite("x3.pgm", buf, Mp, Np);
    }
   /* 
    
    int right,left;
    right = worldRank+1;
    left = worldRank-1;
    if(right == worldSize) right=MPI_PROC_NULL;
    if(left == -1) left =MPI_PROC_NULL;
    //printf("p%d left=%d right=%d \n",worldRank,left,right);
    //printf("Mp = %d\n",Mp);
    
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
    
    
    // Set fixed boundary conditions on the left and right sides 
    
    for (j=1; j < Np+1; j++)
    {
        // compute sawtooth value
        
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
        
        // Implement periodic boundary conditions on bottom and top sides 
        MPI_Isend(&old[Mp][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request);
        MPI_Irecv(&old[0][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request);
        
        //printf("p%d send to %d \n",worldRank,left);
        MPI_Isend(&old[1][1],Np,MPI_DOUBLE,left,0,MPI_COMM_WORLD,&request2);
        MPI_Irecv(&old[Mp+1][1],Np,MPI_DOUBLE,right,0,MPI_COMM_WORLD,&request2);
        
        MPI_Wait(&request,&status);
        MPI_Wait(&request2,&status);
        
        
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
    */
    
    if(worldRank !=0)
    {
        //printf("--->p%d buf[0][0]=%f\n",worldRank,buf[0][0]);
        MPI_Isend(&buf[0][0],Mp*Np,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request);
    }
    
    
    if(worldRank == 0)
    {
        double masterbuf2[M][N];
        double buf1[Mp][Np];
        double buf2[Mp][Np];
        double buf3[Mp][Np];
        MPI_Request request01;
        MPI_Request request02;
        MPI_Request request03;
        
        
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
        
        int cnt;
        for(r =0; r < M; r++){
            for(c=0; c< N; c++)
            {
                if((masterbuf2[r][c] != masterbuf[r][c]))
                {
                    cnt++;  
                }
            }
        }
        printf("cnt=%d\n",cnt);
        
        filename="imagenew192x128.pgm";
        printf("\nWriting <%s>\n", filename); 
        pgmwrite(filename, masterbuf2, M, N);
        
    }
    
    MPI_Finalize();
    return 0;
} 

double boundaryval(int i, int m)
{
    double val;
    
    val = 2.0*((double)(i-1))/((double)(m-1));
    if (i >= m/2+1) val = 2.0-val;
    
    return val;
}
