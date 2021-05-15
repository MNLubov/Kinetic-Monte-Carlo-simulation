#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>

#include "array.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"


//building list of neighbours function for the diamond/zb structure
void NeighbourListZB2D ( void )
{
  int i, j, k, l, n; 
  int NNx, NNz;  // number of atoms along X,Z directions in the simulation box in function  'NeighbourList()'
  int NSimulationBoxTotal;    // total number of the atoms in the simulation box

  NeighborListMemroyAllocation();
  //printf("NeighborListMemroyAllocation\n");
  //getch();
  
  NNx = BoxSizeX * UnitCellAtomsNumber;
  NNz = BoxSizeZ * UnitCellAtomsNumber;
  NSimulationBoxTotal = BoxSizeX * BoxSizeZ * UnitCellAtomsNumber;
  
  //diamond/zb structure
  for(i = 1; i <= NSimulationBoxTotal; i++) 
  {
    j = i % NNx;           //current number of the atom's column (column width is 1 atom)
    if( j == 0)
      j = NNx;
    if((i % 2) == 1)
    {
      n = (int)(i / NNx);
      k = 2 * n + 1;  //current number of the atom's row            
    }
    if((i % 2) == 0)
    {        
      n = (int)(i / NNx) + 1;
      if( i % NNx == 0)
        n = i / NNx;
      k = 2 * n;              
    }  
    l = i * 4;                // current number in the NearestNeighbor and NextNearestNeighbor lists i.e. (NearestNeighbornumber = NextNearestNeighbornumber)         

    //printf("i = %d  j = %d  k = %d\n",i, j, k);
    //getch();
    //traversal order of the neighbors is clockwise
    // calculation of the neighbors 
    if(j != 1 && j != NNx && k != 1 && k != NNz)
    {
      //odd numbers
      if(i % UnitCellAtomsNumber == 1)
      {
        NearestNeighborList[l - 3] = (i - NNx) + 1;
        NearestNeighborList[l - 2] = (i - NNx) - 1;
        NearestNeighborList[l - 1] =  i - 1;//
        NearestNeighborList[l] = i + 1;//
        /*
        NearestNeighborList[l - 3] = i + 1;
        NearestNeighborList[l - 2] = (i - NNx) + 1;
        NearestNeighborList[l - 1] = (i - NNx) - 1;
        NearestNeighborList[l] = i - 1;
        */
      }  
      //even numbers
      if(i % UnitCellAtomsNumber == 0)
      {
        NearestNeighborList[l - 3] = i + 1;
        NearestNeighborList[l - 2] = i - 1;
        NearestNeighborList[l - 1] = (i + NNx) - 1;//
        NearestNeighborList[l] = (i + NNx) + 1;//(i + NNx) - 1;
        /*
        NearestNeighborList[l - 3] = (i + NNx) + 1;
        NearestNeighborList[l - 2] = i + 1;
        NearestNeighborList[l - 1] = i - 1;
        NearestNeighborList[l] = (i + NNx) - 1;
        */
      }  
      // calculation of the next nearest neughbors
      if(j == (NNx - 1))
        NextNearestNeighborList[l - 3] = (i - NNx) + 2;
      if(j != (NNx - 1))      
        NextNearestNeighborList[l - 3] = i + 2;     
        
      if(k == 2)
        NextNearestNeighborList[l - 2] = (NSimulationBoxTotal - NNx) + i;
      if(k != 2)  
        NextNearestNeighborList[l - 2] = i - NNx;
      
      if(j == 2)  
        NextNearestNeighborList[l - 1] = (i + NNx) - 2;     
      if(j != 2)
        NextNearestNeighborList[l - 1] = i - 2;

      if(k == (NNz - 1))
        NextNearestNeighborList[l] = i - (NSimulationBoxTotal - NNx); 
      if(k != (NNz - 1))
        NextNearestNeighborList[l] = i + NNx;
    }
    if(j == 1 && k != 1)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = (i - NNx) + 1;
      NearestNeighborList[l - 2] = i - 1;
      NearestNeighborList[l - 1] = (i + NNx) - 1; //
      NearestNeighborList[l] = i + 1;
      /*
      NearestNeighborList[l - 3] = i + 1;
      NearestNeighborList[l - 2] = (i - NNx) + 1;
      NearestNeighborList[l - 1] = i - 1;
      NearestNeighborList[l] = (i + NNx) - 1;
      
      */
      
      //printf("i = %d  j = %d  k = %d\n",i, j, k);
      //getch();
      // calculation of the next nearest neughbors
      NextNearestNeighborList[l - 3] = i + 2;
      NextNearestNeighborList[l - 2] = i - NNx;
      NextNearestNeighborList[l - 1] = (i + NNx) - 2;
      if(k == (NNz - 1))
        NextNearestNeighborList[l] = 1;
      if(k != (NNz - 1))
        NextNearestNeighborList[l] = i + NNx;      
      }
    if(k == 1 && j != 1)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = i + (NSimulationBoxTotal - NNx) + 1;
      NearestNeighborList[l - 2] = i + (NSimulationBoxTotal - NNx) - 1;
      NearestNeighborList[l - 1] = i - 1;
      NearestNeighborList[l] = i + 1;
      /*
      NearestNeighborList[l - 3] = i + 1;
      NearestNeighborList[l - 2] = i + (NSimulationBoxTotal - NNx) + 1;
      NearestNeighborList[l - 1] = i + (NSimulationBoxTotal - NNx) - 1;
      NearestNeighborList[l] = i - 1;
      */
      // calculation of the next nearest neughbors
      if(j == (NNx - 1))
        NextNearestNeighborList[l - 3] = 1;
      if(j != (NNx - 1))
        NextNearestNeighborList[l - 3] = i + 2;
      NextNearestNeighborList[l - 2] = i + (NSimulationBoxTotal - NNx);
      NextNearestNeighborList[l - 1] = i - 2;
      NextNearestNeighborList[l] = i + NNx;         
    }
    if(k == NNz && j != NNx)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = i + 1;
      NearestNeighborList[l - 2] = i - 1;
      NearestNeighborList[l - 1] = i - (NSimulationBoxTotal - NNx) - 1;
      NearestNeighborList[l] = i - (NSimulationBoxTotal - NNx) + 1;
      /*
      NearestNeighborList[l - 3] = i - (NSimulationBoxTotal - NNx) + 1;
      NearestNeighborList[l - 2] = i + 1;
      NearestNeighborList[l - 1] = i - 1;
      NearestNeighborList[l] = i - (NSimulationBoxTotal - NNx) - 1;
      */
      // calculation of the next nearest neughbors
      NextNearestNeighborList[l - 3] = i + 2;
      NextNearestNeighborList[l - 2] = i - NNx;
      if(j == 2)
        NextNearestNeighborList[l - 1] = NSimulationBoxTotal;         
      if(j != 2)
        NextNearestNeighborList[l - 1] = i - 2;   
      NextNearestNeighborList[l] = i - (NSimulationBoxTotal - NNx);          
     }
    if(j == NNx && k != NNz)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = (i - NNx) + 1;
      NearestNeighborList[l - 2] = i - 1;
      NearestNeighborList[l - 1] = (i + NNx) - 1;
      NearestNeighborList[l] = i + 1;
      /*
      NearestNeighborList[l - 3] = i + 1;
      NearestNeighborList[l - 2] = (i - NNx) + 1;
      NearestNeighborList[l - 1] = i - 1;
      NearestNeighborList[l] = (i + NNx) - 1;
      */
      //printf("i = %d  j = %d  k = %d\n",i, j, k);
      //getch();
      // calculation of the next nearest neughbors
      NextNearestNeighborList[l - 3] = (i - NNx) + 2;
      if(k == 2)
        NextNearestNeighborList[l - 2] = NSimulationBoxTotal;
      if(k != 2)
        NextNearestNeighborList[l - 2] = i - NNx;
      NextNearestNeighborList[l - 1] = i - 2;
      NextNearestNeighborList[l] = i + NNx;         
    }
    if(j == 1 && k == 1)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = (NSimulationBoxTotal - NNx) + 2;
      NearestNeighborList[l - 2] = NSimulationBoxTotal;
      NearestNeighborList[l - 1] = NNx;
      NearestNeighborList[l] = 2;
      /*
      NearestNeighborList[l - 3] = 2;
      NearestNeighborList[l - 2] = (NSimulationBoxTotal - NNx) + 2;
      NearestNeighborList[l - 1] = NSimulationBoxTotal;
      NearestNeighborList[l] = NNx;
      */
      // calculation of the next nearest neughbors
      NextNearestNeighborList[l - 3] = i + 2;
      NextNearestNeighborList[l - 2] = (NSimulationBoxTotal - NNx) + 1;
      NextNearestNeighborList[l - 1] = NNx - 1;
      NextNearestNeighborList[l] = i + NNx;     
    }
    if(j == NNx && k == NNz)
    {
      // calculation of the nearest neighbors 
      NearestNeighborList[l - 3] = (NSimulationBoxTotal - NNx) + 1;
      NearestNeighborList[l - 2] = i - 1;
      NearestNeighborList[l - 1] = NNx - 1;
      NearestNeighborList[l] = 1;
      /*
      NearestNeighborList[l - 3] = 1;
      NearestNeighborList[l - 2] = (NSimulationBoxTotal - NNx) + 1;
      NearestNeighborList[l - 1] = i - 1;
      NearestNeighborList[l] = NNx - 1;
      */
      // calculation of the next nearest neughbors
      NextNearestNeighborList[l - 3] = (NSimulationBoxTotal - NNx) + 2;
      NextNearestNeighborList[l - 2] = NSimulationBoxTotal - NNx;
      NextNearestNeighborList[l - 1] = i - 2;
      NextNearestNeighborList[l] = NNx;     
    }
    //printf("%d  %d  %d  %d  %d\n",i, NearestNeighborList[l - 3], NearestNeighborList[l - 2], NearestNeighborList[l - 1], NearestNeighborList[l]);
    //getch();
    //printf("%d  %d  %d  %d  %d\n",i, NextNearestNeighborList[l - 3], NextNearestNeighborList[l - 2], NextNearestNeighborList[l - 1], NextNearestNeighborList[l]);
    //getch();
  }
    
}

int NearestNeighborCalc( int SiteNumber )
{
  int i, j, k, l;
  int NearestSum;  // number of the nearest neighbor for the atoms from the FreeSite  array  
  
  NearestSum = 0;
  for(i = (NearestNeighbornumber - 1); i >=0; i--)
  {
    k = SiteNumber * 4 - i;	
    l = NearestNeighborList[k];       // nearest neighbors of the site
    if(SimulationBox[l] >= 1)
    {
      NearestSum = NearestSum + 1;
    }     
  }
  
  return NearestSum;
}

int NearestNeighborBeforeProcess( int SiteNumber )
{
	int i, j, k, l;
	int NearestSum, AtomType;
	
	NearestSum = 0;
	
	//setting site from which atom moved to occupied 
	//AtomType = SimulationBox[ProcessedAtomNumber];
	//SimulationBox[ProcessedAtomNumber] = 1;
	for(i = (NearestNeighbornumber - 1); i >=0; i--)
    {
      k = SiteNumber * 4 - i;	
      l = NearestNeighborList[k];       // nearest neighbors of the site
      if(SimulationBox[l] >= 1 || AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] > 0)
      {
        NearestSum = NearestSum + 1;
      }     
    }
    
    return NearestSum;
}

int NextNearestNeighborBeforeProcess( int SiteNumber )
{
	int i, j, k, l;
	int NextNearestSum, AtomType;
	
	NextNearestSum = 0;
	
	//setting site from which atom moved to occupied 
	//AtomType = SimulationBox[ProcessedAtomNumber];
	//SimulationBox[ProcessedAtomNumber] = 1;
	for(i = (NearestNeighbornumber - 1); i >=0; i--)
    {
      k = SiteNumber * 4 - i;	
      l = NextNearestNeighborList[k];       // nearest neighbors of the site
      if(SimulationBox[l] >= 1 || AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] > 0)
      {
        NextNearestSum = NextNearestSum + 1;
      }     
    }
        
    return NextNearestSum;
}



void NeighborListCalc ( int DepositedStructureIndex )
{
  //DepositedStructureIndex = 1 means crystal structures of the substrate and deposited material are the same
  if(SubstrateType == 5 && DepositedStructureIndex == 1)
  {
  	NeighbourListZB2D();
  }
   
}

void SubstrateNeighbors( void )
{
	int i, j, k;
}


