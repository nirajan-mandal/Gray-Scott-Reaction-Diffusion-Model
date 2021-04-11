// Modeling Gray-Scott Reaction-Diffusion equation
//FTCS numerical method is used in the modeling process
//the boundry wraps such that the recation occours on a torus surface
//Nirajan Mandal, April, 2009
//original code: B. Schumacher -- 13 April 2009

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "mpi.h"
using namespace std;

double d=0.01, h=0.5;	                       //delta, del x = del y = h
double Du=0.16, Dv=0.08, f=0.046, k=0.061;	   //Gray-Scot Model constants: D=difusion constant ration
double t=0, tmax=1000;	                       //time parameters
int x=10;                                      //outputs jpg file and data every x time unit
const int XSIZE =100;	                       //size of grid (i had issues if I increased grid size to more than 150
                                               //this must have something to do with the memory 
int filenumber=100000;  //partof filenumber for jpg graphs, it is assumed that no more than 100000 files will be written

int oldt=0,newt=1, u=0, v=1;                   //old grid, updated grid, u-grid and v-grid  
double grid[2][2][XSIZE+2][XSIZE+2];  //old-new, u-v, x-size, y-sixe, hence 4 dimensional
double sendedgeleft[XSIZE+2],sendedgeright[XSIZE+2], sendtop[XSIZE+2], sendbottom[XSIZE+2];  //used for boundry condition
double recvedgeleft[XSIZE+2],recvedgeright[XSIZE+2], recvetop[XSIZE+2],recvebottom[XSIZE+2]; //used for boundry condition
int me, myleft, myright, size;                  //keeps track of who is processors

ofstream outd, outfu, outfv, gnuf;              //writing file streams

double mysumu=0, mysumv=0;                    //used for density calculations
double allsumu=0, allsumv=0;
double densityu=0, densityv=0;

int i, j, m, n;                             //used for index and other loops

void initialize(void) // initialize grid;
{
   srand((unsigned)time(0)+me);    //seed for random numbers
   for (m=1;m<XSIZE+1;m++)
      for (n=1;n<XSIZE+1;n++)
      {
         grid[oldt][u][m][n] =1;
	     grid[oldt][v][m][n] =0;
      }

  for (m=10*(1+rand()%3);m<50+10*(rand()%4);m++)    //random sizes are created
      for (n=10*(1+rand()%3);n<50+10*(rand()%3);n++)
      {
         grid[oldt][u][m][n] =0.5;
	    grid[oldt][v][m][n] =0.25;
      }
}//end initialize

void boundary(void)              // handle the boundary by passing messages
{
   MPI_Status stat;
   int i, s;
 for(s=0; s<2; s++)  //twice for U and V
 {
   for (i=0;i<XSIZE+2;i++) //stores the four boundry data 
   {
      sendedgeleft[i] = grid[oldt][s][i][1];
      sendedgeright[i] = grid[oldt][s][i][XSIZE];  
      sendtop[i]=grid[oldt][s][1][i]; 
      sendbottom[i]=grid[oldt][s][XSIZE][i];   
   }
            //sends all the four boundry data to the respective process
   MPI_Send(sendedgeleft,XSIZE+2,MPI_DOUBLE,myleft,1,MPI_COMM_WORLD);
   MPI_Send(sendedgeright,XSIZE+2,MPI_DOUBLE,myright,1,MPI_COMM_WORLD);
   MPI_Send(sendtop,XSIZE+2,MPI_DOUBLE,me,1,MPI_COMM_WORLD);
   MPI_Send(sendbottom,XSIZE+2,MPI_DOUBLE,me,1,MPI_COMM_WORLD);
            //receives all the four boundry data to the respective process
   MPI_Recv(recvedgeleft,XSIZE+2,MPI_DOUBLE,myleft,1,MPI_COMM_WORLD,&stat);
   MPI_Recv(recvedgeright,XSIZE+2,MPI_DOUBLE,myright,1,MPI_COMM_WORLD,&stat);
   MPI_Recv(recvebottom,XSIZE+2,MPI_DOUBLE,me,1,MPI_COMM_WORLD,&stat);
   MPI_Recv(recvetop,XSIZE+2,MPI_DOUBLE,me,1,MPI_COMM_WORLD,&stat);

   for (i=0;i<XSIZE+2;i++)  //writes the four boundry data
   {
      grid[oldt][s][i][0] = recvedgeleft[i];
      grid[oldt][s][i][XSIZE+1] = recvedgeright[i];
      grid[oldt][s][0][i] = recvetop[i];
      grid[oldt][s][XSIZE+1][i] = recvebottom[i];
   }

 }//end for two grids
}//end boundry

void updatecells(void)            // apply rule to old grid to get new grid using FTCS integration method
{
   double u1, u2, v1, v2;
   int i,j;
   boundary();                    // take care of boundary stuff
   for (i=1;i<XSIZE+1;i++)
      for (j=1;j<XSIZE+1;j++)
      {
	u1=Du*(grid[oldt][u][i+1][j]+grid[oldt][u][i-1][j]+grid[oldt][u][i][j+1]+grid[oldt][u][i][j-1]-4*grid[oldt][u][i][j])*d/(h*h);    //del square part
	u2=d*(f*(1-grid[oldt][u][i][j])-grid[oldt][u][i][j]*grid[oldt][v][i][j]*grid[oldt][v][i][j]); 
	grid[newt][u][i][j]=grid[oldt][u][i][j]+u1+u2;  //updates u-grid
	
	v1=Dv*(grid[oldt][v][i+1][j]+grid[oldt][v][i-1][j]+grid[oldt][v][i][j+1]+grid[oldt][v][i][j-1]-4*grid[oldt][v][i][j])*d/(h*h);    //del square part
	v2=d*(grid[oldt][u][i][j]*grid[oldt][v][i][j]*grid[oldt][v][i][j]-(f+k)*grid[oldt][v][i][j]); 
	grid[newt][v][i][j]=grid[oldt][v][i][j]+v1+v2;  //updates v-grid   
        
      }
   oldt = (oldt+1) % 2;           // new becomes old, old becomes new
   newt = (newt+1) % 2;
}//end updatecells

void get_density(void)              // getting the average density
{
   mysumu=0, mysumv=0;              //resetting all variables to zero
   allsumu=0, allsumv=0;
   densityu=0, densityv=0;    
   
  for (i=1; i<=XSIZE; i++)      //gets the sum of all the grid
      for (j=1; j<=XSIZE; j++)
      {
         mysumu = mysumu+grid[u][0][i][j];
         mysumv = mysumv+grid[v][0][i][j];
      }
            //collects the sum
   MPI_Reduce(&mysumu,&allsumu,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&mysumv,&allsumv,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   
   allsumu = mysumu;
   allsumv = mysumv;

   if (me == 0) //writes teh density
   {
      densityu = allsumu/(size*XSIZE*XSIZE);
      densityv = allsumv/(size*XSIZE*XSIZE);
      outd<<t<<" "<<densityu<<" "<<densityv<<endl;
   }
   
}//end get_density

void graphit(void)
{
    MPI_Status stat;
    if (me == 0)
   {
        double paper[size][XSIZE+2][XSIZE+2]; //scrath grid to receive the grid data from all process 
      
        outfu.open("/home/mandal/gray_scott/grayscott_data_u.txt");  //u-grid data
        outfv.open("/home/mandal/gray_scott/grayscott_data_v.txt");	//v-grid data
        
        MPI_Send(grid[0][u],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,0,1,MPI_COMM_WORLD);  //send u-grid to itself
    
       for (m=0;m<size;m++)
       {
            MPI_Recv(paper[m],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,m,1,MPI_COMM_WORLD,&stat);	// receive grid from all
        }
   	   for(i=1; i<XSIZE+1; i++)
   	   {
          
          for(m=0; m<size; m++)
          {
                for(j=1; j<XSIZE+1; j++)
                {
                 outfu<<paper[m][i][j]<<" ";  //writes the u-grid data for each position
                }
             }
        outfu<<endl;                           //moves to the next line 
      } //finished writing u-grid
      outfu.flush();                        //clears the writing buffer
      outfu.close();                        //closes u-grid data file
      MPI_Barrier(MPI_COMM_WORLD);
      
      //starts writing V-grid 
      MPI_Send(grid[0][v],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,0,1,MPI_COMM_WORLD);  //send v-grid to itself
      for (m=0;m<size;m++)
       {
            MPI_Recv(paper[m],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,m,1,MPI_COMM_WORLD,&stat);	// receive grid from all
        }
   	   for(i=1; i<XSIZE+1; i++)
   	   {
          for(m=0; m<size; m++)
          {
                for(j=1; j<XSIZE+1; j++)
                {
                 outfv<<paper[m][i][j]<<" ";    //writes the v-grid data for each position
                }
             }
        outfv<<endl;                            //moves to the next line
      } //finished writing v-grid
      outfv.flush();                            //clears the writing buffer
      outfv.close();                            //closes u-grid data file
      MPI_Barrier(MPI_COMM_WORLD);
   }//end if
   
   else  //if process=0, then they all send their data to process=0
   {
      MPI_Send(grid[0][u],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,0,1,MPI_COMM_WORLD);  
      MPI_Barrier(MPI_COMM_WORLD);   
      MPI_Send(grid[0][v],(XSIZE+2)*(XSIZE+2),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);   
    }//end else 
    
    
    if(me==0)   //processor zero writes the jpeg grapg file and the inital conditions that ran the code
    {
        //these lines the gnu command to create different image file name
 
        gnuf.open("/home/mandal/gray_scott/grayscott_gnu_command.txt"); 
        gnuf<<"reset"<<endl;
        gnuf<<"set terminal jpeg small font times size 800,600"<<endl;
        gnuf<<"set output '/home/mandal/gray_scott/video_files/graph"<<filenumber<<".jpg'"<<endl;
        gnuf<<"set multiplot"<<endl;
        gnuf<<"unset key "<<endl;
        gnuf<<"set view map"<<endl;
        gnuf<<"set size ratio "<<1.0/size<<endl;  //this is the height to width raio
        gnuf<<"set lmargin at screen 0.06"<<endl;
        gnuf<<"set rmargin at screen 0.85"<<endl;
        gnuf<<"set bmargin at screen 0.55"<<endl;
        gnuf<<"set tmargin at screen 0.95"<<endl;
        gnuf<<"set title ' grid u : t="<<t<<" / "<<tmax<<"'"<<endl;
        gnuf<<"splot \"/home/mandal/gray_scott/grayscott_data_u.txt\" matrix with image"<<endl;
        
        gnuf<<"set title ' grid v : t="<<t<<" / "<<tmax<<"'"<<endl;   
        gnuf<<"set lmargin at screen 0.06"<<endl;
        gnuf<<"set rmargin at screen 0.85"<<endl;
        gnuf<<"set bmargin at screen 0.05"<<endl;
        gnuf<<"set tmargin at screen 0.55"<<endl;  
        gnuf<<"splot \"/home/mandal/gray_scott/grayscott_data_v.txt\" matrix with image"<<endl;   
        gnuf<<"unset multiplot"<<endl;
        gnuf<<"set output"<<endl;
        gnuf<<"exit"<<endl;
        gnuf.flush();
        gnuf.close();
                //runs GNUplot and reads the commad file to create the jpg file
        system("/usr/local/bin/gnuplot /home/mandal/gray_scott/grayscott_gnu_command.txt");  
        
        filenumber++; //increases the file number so that it is in sequence to make video
     }//end if

}//end graphit    

main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);                // Initialize MPI
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   
   myleft = (me + size - 1) % size;
   myright = (me + 1) % size;
   if(me==0)
      outd.open("/home/mandal/gray_scott/density.txt"); //only process=0 writes the density file
   
   initialize();                 // get ready
   graphit();                   //graph the initial condition
   get_density();               //write teh initial density
   
   for (t=d; t<tmax; 1)        // updating loop
   {
      updatecells();
      if(int(t*(1.0/d)+1)%(x*int((1.0/d)+1))==0)  //outputs jpg file and density every x time unit
      {
             graphit();  
             get_density();
       }
      MPI_Barrier(MPI_COMM_WORLD);                //stops until every process reaches this point  
      t+=d;                         //increase time
   }
   
   boundary(); //useful to see the final result
      
   if(me==0)
   {
   outfv.flush();   //clears writing buffer
   outfv.close();   //close density file
    }  
    
  MPI_Finalize();               // all MPI done!
  
}//end main
