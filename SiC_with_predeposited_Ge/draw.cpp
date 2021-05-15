#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <graphics.h>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"
#include "draw.h"

int dx, dy;

void DrawMesh( void )
{
  int i,j;
  int cell;
  int k, l;
  
  
  
  dx = 800 / xcell;
  dy = 800 / ycell; 
   
  
  setcolor( 15 );
  setfillstyle(0,100);
  //pieslice( 250, 250, 0, 0, 10 );
  //fillellipse( 100, 100, 5, 5);

  for(i = 1; i <=xcell; i++)
    for(j = 1; j <= ycell; j++)
     {
       fillellipse( 10 + dx * (i-1), 790 - (dy * (j-1)), 5, 5);
       fillellipse( 10 + (dx * (i-1) + dx/2), 790 - (dy * (j-1) + dy/2), 5, 5);
     }
  
  //testing of the mesh  
  //line(100,750,100,100);
  //line(100,750,750,750);
  //return 0;
    
}

void DrawEvolution ( int number )
{
  int i,j;
  
  
  
  setfillstyle(1,15);
  //NeighborList(number);
  
  for(i = 0; i <= xcell*ycell*2; i++)
   { 
      if(adsite[i] == 1 && i % 2 == 1) 
      {
        NeighborList(i);
        fillellipse( 10 + dx * atomcell[i][0], 790 - dy * atomcell[i][1], 5, 5);
      }
      if(adsite[i] == 1 && i % 2 == 0) 
      {
        NeighborList(i);
        fillellipse( 10 + dx * atomcell[i][0] + dx/2, 790 - dy * atomcell[i][1]+dy/2, 5, 5);
      }
   }
  
}
