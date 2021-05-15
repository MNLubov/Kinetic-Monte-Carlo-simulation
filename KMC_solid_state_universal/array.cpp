#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <conio.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>
//#include<vector.h> 

#include "array.h"
#include "parameters.h"

using namespace std;

void SimulationBoxMemoryAllocation( void )
{
  
  SimulationBox = new int[SimulationBoxAtomsNumber];
  SimulationBoxCoordinates = new double[SimulationBoxAtomsNumber * 3];
  AtomSitesConfiguration = new double[SimulationBoxAtomsNumber];
}

void SimulationBoxMemoryCleaning( void )
{
  //free(SimulationBox);
  // calculation simulation box parameters  
  delete [] SimulationBox;  
  delete [] SimulationBoxCoordinates;      
}


void NeighborListMemroyAllocation( void )
{
   //cout << "SimulationBoxAtomsNumber=" <<SimulationBoxAtomsNumber; 
   //getch();
   NearestNeighborList = new int[SimulationBoxAtomsNumber * NearestNeighbornumber];
   NextNearestNeighborList = new int[SimulationBoxAtomsNumber * NextNearestNeighbornumber] ;  
}


void SimulationBoxConfigurationMemoryAllocation( void )
{
  //SimulationBoxConfiguration = new int[BoxConfigurationNumber];
  //vector< vector<int> > SimulationBoxConfiguration;
  AtomConfigurationNumber = new int [2 * (SimulationBoxAtomsNumber - SubstrateBottomNumber) + 1];
  
}

void NeighborListMemroyCleaning( void )
{
   delete [] NearestNeighborList;
   delete [] NextNearestNeighborList;
}

void SurfaceMapMemoryAllocation ( void )
{
   SurfaceMap = new int [BoxSizeX * UnitCellAtomsNumber + 1];
}

void SurfaceMapMemoryCleaning ( void )
{
   delete [] SurfaceMap;
}

void AtomSitesConfigurationMemoryCleaning( void )
{
   delete [] AtomSitesConfiguration;
}

void SimulationBoxConfigurationMemoryCleaning( void )
{
   //delete [] SimulationBoxConfiguration;
   //vector<tempObject>().swap(tempVector);
   vector< vector<int> >().swap(SimulationBoxConfiguration);
   vector<double>().swap(EnergyList);
   vector<double>().swap(ProcessRate);
   vector<double>().swap(DiffusionRate);
   //vector<double>().swap(AtomConfigurationNumber);
   delete [] AtomConfigurationNumber;
   //vector<double>().swap(AtomConfigurationNumber);
}

void MemoryAllocation( void )
{
  /* allocating memory for the simulation data */
  SimulationBoxMemoryAllocation();
  printf("SimulationBoxMemoryAllocation\n");
  //NeighborListMemroyAllocation();
  //printf("NeighborListMemroyAllocation\n");
  //SurfaceMapMemoryAllocation();   
  //printf("allocation3\n");
}

void MemoryCleaning( void )
{
  SimulationBoxMemoryCleaning();
  NeighborListMemroyCleaning();
  SurfaceMapMemoryCleaning();
  AtomSitesConfigurationMemoryCleaning();
  SimulationBoxConfigurationMemoryCleaning();
}

void SetStaticArray( void )
{
  int i, j, k;
  
  for(i = 0; i <= 10; i++)
    for(j = 0; j <= 20; j++)
      BondEnergy[i][j] = 0;
  
  for(i = 0; i <= 10; i++)  
  {
  	Flux[i] = 0;
	HopEnergy[i] = 0;
  	EvaporationEnergy[i] = 0;
  }
    
    
  for(i = 0; i <= 99; i++)  
  {
  	NearestBondEnergy[i] = 0;
    NextNearestBondEnergy[i] = 0;
    AtomFrequency[i] = 0;
  }
  
}


