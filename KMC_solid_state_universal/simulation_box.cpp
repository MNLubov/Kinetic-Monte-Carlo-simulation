#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"
#include "deposition.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"

using namespace std;

//function for generation simulation box
void SimulationBoxGeneration( void )
{
  int i, j, k, l, m, n;
  int Nx, Ny, Nz;  // number of atoms along X,Y,Z directions in the simulation box
  int NSubsZ; // number of atoms along Z directions in the substrate box
  double Lx, Ly, Lz; // length of the simulation box along X,Y,Z directions
  
  //setting unit cell parameters
  CrystalStructureType = SubstrateType;
  UnitCellTermination = SubstrateTermination;
  
  
  if( Dimension == 3)
    UnitCellGeneration3D(CrystalStructureType);       
  
  if( Dimension == 2)
  { 
    SimulationBoxMemoryAllocation();
    //printf("SimulationBoxGeneration\n");

    Nx = BoxSizeX * UnitCellAtomsNumber;
    Nz = BoxSizeZ * UnitCellAtomsNumber;
    Lx = (BoxSizeX - 1) * UnitCellSizeX + UnitCellLengthX;
    Lz = (BoxSizeZ - 1) * UnitCellSizeX + UnitCellLengthZ;
    
    //setting coordinates of the substrate sites (occupied) in the simulation box
    for(i = 1; i <= SubstrateSizeZ; i++) 
      for(j = 1; j <= Nx; j++)
      {
        //setting atoms type in global array of free and occupied states
        l = j + (i - 1) * Nx;                   // current number in the SimulationBox[] array
        k = (l - 1) % UnitCellAtomsNumber + 1;  // current number of the atom in the unit cell
        m = l * 3;                              // current number in the SimulationBoxCoordinates[] array
        n = (int)((j - 1) / UnitCellAtomsNumber); // current number of the (row - 1)  in the simulatioon box
        
        SimulationBox[l] = UnitCellAtomsType[k];
        //printf("i=%d  j = %d   k = %d   l = %d   m = %d  n = %d\n", i,j,k,l,m, n);
        SimulationBoxCoordinates[m - 2] = UnitCellCoordinates[k][1] + UnitCellSizeX * n;
        SimulationBoxCoordinates[m - 1] = UnitCellCoordinates[k][2];
        SimulationBoxCoordinates[m] = UnitCellCoordinates[k][3] + UnitCellSizeZ * (i - 1);
      }   
    
    //setting coordinates of the free sites in the simulation box
    for(i = SubstrateSizeZ + 1; i <= BoxSizeZ; i++) 
      for(j = 1; j <= Nx; j++)
      {
        //setting atoms type in global array of free and occupied states
        l = j + (i - 1) * Nx;                   // current number in the SimulationBox[] array
        k = (l - 1) % UnitCellAtomsNumber + 1;  // current number of the atom in the unit cell
        m = l * 3;                              // current number in the SimulationBoxCoordinates[] array
        n = (int)((j - 1) / UnitCellAtomsNumber); // current number of the (column - 1)  in the simulatioon box
        //printf("i=%d  j = %d   k = %d   l = %d   m = %d  n = %d\n", i,j,k,l,m, n);
        SimulationBox[l] = 0;
        SimulationBoxCoordinates[m - 2] = UnitCellCoordinates[k][1] + UnitCellSizeX * n;
        SimulationBoxCoordinates[m - 1] = UnitCellCoordinates[k][2];
        SimulationBoxCoordinates[m] = UnitCellCoordinates[k][3] + UnitCellSizeZ * (i - 1);
      }   
  }  
  
   for(i = 1; i <= SimulationBoxAtomsNumber; i++) 
   {
     j = i * 3;
     SimulationBoxCoordinates[j - 2] = SimulationBoxCoordinates[j - 2] - Lx / 2;
     //SimulationBoxCoordinates[j] = SimulationBoxCoordinates[j] - Lz / 2;
   }
   
}

//generation of the flat deposition map for the simulation box with ZB structure
void SurfaceMapGenerationZB( int RoughnessIndex )
{
  int i, j, k;
  int NSurfX; // number of deposition/diffusion sites on the surface 
  
  SurfaceMapMemoryAllocation();   
  //printf("SurfaceMapMemoryAllocation\n");
  
  NSurfX = BoxSizeX;
  
  j = BoxSizeX * SubstrateSizeZ * UnitCellAtomsNumber;
  
  //generating flat interafce
  if(RoughnessIndex == 1)
  {
    //building of the DepositionMap[] array
    for(i = 1; i <= NSurfX; i++)
    {
      k = i * 2 - 1;
      SurfaceMap[i] = j + k;
      //printf("!!!!%d  %d\n", i, j + k);
      //getch(); 
    } 
    SurfaceMapLength = NSurfX;  
  } 
  
}

void SurfaceMapGeneration( int RoughnessIndex )
{
  if(SubstrateType == 5)
    SurfaceMapGenerationZB( RoughnessIndex );
}

void SurfaceMapEvolution2D( int SurfaceSiteNumber, int SurfaceMapNumber, int ProcessIndex )
{
  int i, j, k, l, m, n, p;
  int S; // number of the nearest neighbors of the potential site in SurfaceMap[]
  
                                    /* NB! #1 */
  /*
  evolution of the SurfaceMap in the case of the surface diffusion
  is equivalent to the 'evaporation + deposition' processes
  */
  
  									/* NB! #2 */
  /*
  when SurfaceMap updates in the case of the evaporation 									
  it is considered that evaporated atom doesn't leave his position!
  So setting corresponding element of the SimulationBox array to zero 
  must be done after SurfaceMap update!
  (to be changed in later versions)
  */
  
  
  //  
  //printf("Before %d  \n", SurfaceMapLength);
  //getch();

  /* working (not-optimized) version, 22.05.2014 */
  /*
  for(i = 1; i <= SurfaceMapLength; i++ )
  {
    
    if(SurfaceMap[i] == SurfaceSiteNumber && i != SurfaceMapLength)
    {
      SurfaceMap[i] = SurfaceMap[SurfaceMapLength];
      //printf("Main %d  %d  %d\n", SurfaceSiteNumber, i, SurfaceMap[SurfaceMapLength]);
      //getch();  
      SurfaceMapLength = SurfaceMapLength - 1;
    }
      
    if(SurfaceMap[i] == SurfaceSiteNumber && i == SurfaceMapLength)
    {
      SurfaceMap[i] = 0;
      //printf("Main %d  %d  %d\n", SurfaceSiteNumber, i, SurfaceMap[SurfaceMapLength]);
      //getch();  
      SurfaceMapLength = SurfaceMapLength - 1;
    }
  }*/ 
  
  //evilution of the surfaceMap in the case of the deposition
  if(ProcessIndex == 1)
  {
    //updating length of the SurfaceMap[] array
    if( SurfaceMapNumber != SurfaceMapLength)
    {
      SurfaceMap[SurfaceMapNumber] = SurfaceMap[SurfaceMapLength];
      SurfaceMapLength = SurfaceMapLength - 1;
    }
    if( SurfaceMapNumber == SurfaceMapLength)
    {
      SurfaceMap[SurfaceMapNumber] = 0;
      SurfaceMapLength = SurfaceMapLength - 1;
    }
    
    //adding SurfaceMap elements, which is not in the array yet
    for(i = 3; i >= 0; i--)
    {
    //printf("potential\n");
    //getch();
      k = SurfaceSiteNumber * 4 - i;  //number of the element in the NearestNeighborList[]
      l = NearestNeighborList[k];     // number of the site (value of the element in the NearestNeighborList[])which is the nearest neighbor to the SurfaceSiteNumber
      S = 0;
      for(m = 3; m >= 0; m--)
      {
        n = l * 4 - m;                   //number of the element in the NearestNeighborList[]
        p = NearestNeighborList[n];      // number of the site which is nearest to the nearest site to the SurfaceSiteNumber
            
        if(SimulationBox[p] >= 1)        //calculation of the nearest neighbors of the nearest neighbor of the potenrial site in SurfaceMap
        {
          S = S + 1;
        }
      }
      if(S == 1)                          // 
      {
        SurfaceMap[SurfaceMapLength + 1] = l;
        SurfaceMapLength = SurfaceMapLength + 1;
      }
    }
  } 
  
  //evolution of the surfaceMap in the case of the evaporation
  if(ProcessIndex == 2)
  {
     SurfaceMapLength = SurfaceMapLength + 1;
     SurfaceMap[SurfaceMapLength] = SurfaceSiteNumber; 
     //cout << "Map evolution1" << endl;
     //getch();
     
	 //removing SurfaceMap elements, which single NearestNeighbor was evaporated atom
     for(j = 3; j >= 0; j--)
     {
    //printf("potential\n");
    //getch();
      k = SurfaceSiteNumber * 4 - j;  //number of the element in the NearestNeighborList[]
      l = NearestNeighborList[k];     // number of the site (value of the element in the NearestNeighborList[])which is the nearest neighbor to the SurfaceSiteNumber
      S = 0;  
      for(m = 3; m >= 0; m--)
      {
        n = l * 4 - m;                   //number of the element in the NearestNeighborList[]
        p = NearestNeighborList[n];      // number of the site which is nearest to the nearest site to the SurfaceSiteNumber
            
        if(SimulationBox[p] >= 1)        //calculation of the nearestneighbors 
        {
          S = S + 1;
        }
      }
      
     //cout << "Map evolution2" << endl;
     //getch();
      if(S == 1)                          // one nearest neighbor (SurfaceSiteNumber)
      {
        for(i = 1; i <= SurfaceMapLength; i++ )
        {
          if(SurfaceMap[i] == l && i != SurfaceMapLength)
          {
            SurfaceMap[i] = SurfaceMap[SurfaceMapLength];
            SurfaceMapLength = SurfaceMapLength - 1;
          }
          if(SurfaceMap[i] == l && i == SurfaceMapLength)
          {
            SurfaceMap[i] = 0;
            SurfaceMapLength = SurfaceMapLength - 1;
          }  
        }  
      }
    }   
  } 
}

int SurfaceMapNumberCalculation( int SiteNumber )
{
  int i, j, k;
  
  i = 1;
  while(SiteNumber != SurfaceMap[i])
  {
  	i++;
  }	 
  
  return i;
}

int LastAtomCalculation( int AtomNumber, int ProcessIndex )  
{
  int i, k;
  // ProcessIndex = 1 - Deposition
  // ProcessIndex = 2 - Evaporation
  // ProcessIndex = 3 - Diffusion
  if(ProcessIndex == 1 || ProcessIndex == 3)
  {
    if(AtomNumber > LastAtomNumber)
      LastAtomNumber = AtomNumber;
  }
  
  //searching forthe last atom number in the case of the evaporation
  if(ProcessIndex == 2)
  {
    if(AtomNumber == LastAtomNumber)
    {
      i = LastAtomNumber;  
      k = 1;
      while(SimulationBox[i - k] != 1)
      {
        k++;
      }
      LastAtomNumber = i - k;
    }  
  }  
  
  return LastAtomNumber;
}

double AtomSiteEnergy( int AtomNumber )
{
  int i, j, k, l;
  int AtomNumberType, NearestNeighborType, NextNearestNeighborType;
  double SiteEnergy;
  int NearestNeighborsMax; //maximal number of the nearest neighbors for the given crystal structure
  int S;   //
  int NeiborsTemp[9];
  FILE *F;
  
  if(CrystalStructureType == 5)
    NearestNeighborsMax = 4;
  
  S = 0;
  SiteEnergy = 0;
  AtomNumberType = SimulationBox[AtomNumber];
  //cout << "AtomNumberType=" << AtomNumberType << endl;
  //getch();
  
  for(i = 3; i >= 0; i--)
  {
    j = 4 * AtomNumber - i;
    k = NearestNeighborList[j];
    l = NextNearestNeighborList[j];
    if(SimulationBox[k] != 0)
    {
      NearestNeighborType = SimulationBox[k]; 
      SiteEnergy = SiteEnergy + BondEnergy[AtomNumberType][NearestNeighborType];
      S = S + 1;
      //fprintf(FService1, "AtomNumber = %d  NearestNeighborType = %d  S = %d \n", AtomNumber, NearestNeighborType, S);
	  //cout << "SN=" << S << "  SE= "<< SiteEnergy << endl;
      //getch();
    } 
    if(SimulationBox[l] != 0)
    {
      NextNearestNeighborType = SimulationBox[l]; 
      SiteEnergy = SiteEnergy + BondEnergy[AtomNumberType][10 + NextNearestNeighborType];
      //fprintf(FService1, "AtomNumber = %d  NextNearestNeighborType = %d  S = %d \n", AtomNumber, NextNearestNeighborType, S);
      //cout << "SNN=" << S << "  SE= "<< SiteEnergy << endl;
      //getch();
    }
  }
  //printf("S = %d  SiteEnergy = %f\n", S, SiteEnergy);
  //getch();
  //cout << "ERROR" << endl;
  //getch();
  
  if(PhysicalProcessIndex[2] == 1)
    return SiteEnergy;
 
  if(PhysicalProcessIndex[2] == 0)
  {
    if(S == 4)
      return -1.0;
    if(S != 4)
    {
      return SiteEnergy;  
      //cout << SiteEnergy<< endl;
      //getch();
    }
     
  }
 
}

