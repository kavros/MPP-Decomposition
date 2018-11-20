#include "../include/communications.h"


void gather(topology topo,double** masterbuf,double** buf,MPI_Comm comm2d,char* output,int worldSize)
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