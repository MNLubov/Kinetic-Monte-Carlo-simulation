#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>

#include "array.h"
#include "deposition.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"

//function of the deposition only to the sites having nearest neighbor atoms
void DepositionNearest2D( int DepositionSiteNumber, int SurfaceMapNumber, int AtomType)
{
  int i, j, k;
  
  SimulationBox[DepositionSiteNumber] = AtomType;
  printf("%d   %d\n", DepositionSiteNumber, AtomType);
  //getch();
  SurfaceMapEvolution2D( DepositionSiteNumber,  SurfaceMapNumber, 1 );
  //printf("Error\n");   
  //getch();
}

void Evaporation( int EvaporationSiteNumber, int  SurfaceMapNumber )
{
  
  SurfaceMapEvolution2D( EvaporationSiteNumber,  SurfaceMapNumber, 2 );
  SimulationBox[EvaporationSiteNumber] = 0;
  //printf("Error\n");   
  //getch();
}
