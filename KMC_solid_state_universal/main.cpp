#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <vector>


#include "array.h"
#include "configurations.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"


using namespace std;

int main ( void )
{
  int i, j, k;
  int RandNumber1, //random numbers
      RandNumber2; 
      
  int FileErIndex, FileExIndex; //return values
  int FileClose1, FileClose2;   // indicates whether file is closed already or not 
  
  int CurrentTime,
      CurrentStep,
      CurrentAtom;
      
  int ErrorIndex1;
  
  vector <int> Vector1(10);
    
  FILE *FTemp;
  
  //FTemp = fopen("checkRNB","w+");
  /* open service files */
  FService1 = fopen("FS1","w+");
  /* setting initial parameters */
  //set static arrays values to zero before reading from the input files
  SetStaticArray();
  
  // open execution and error log files
  FError = fopen("ErrorLog","w+");
  FExec = fopen("ExecLog","w+");
  
  FileClose1 = 0;
  FileClose2 = 0;
  
  CurrentTime = 0;
  CurrentStep = 0;
  CurrentAtom = 0;
  /* read input data from the file */
  FileErIndex = ReadData();
  
  if(FileErIndex == 2)
  {
  	printf("Inappropriate input parameters. See Error log\n");
  	fclose(FError);
  	FileClose1 = 1;
  	return 0;
  }
  
  
                   /* calculation of the initial state of the system */
  /*                                                                         */
  PreDeposition();
  //getch();
  /* generation of the simulation box */
  printf("        Initialisation of the system");
  printf("\n generation of the simulation box......................Done \n");
  //printf("----------------------------------\n \n");
  SimulationBoxGeneration();
  
  /* initial neighborlist calculation */
  printf(" initial neighborlist calculation......................Done\n");
  //printf("----------------------------------\n \n");
  NeighborListCalc(1);  // DepositedStructureIndex = 1 - substrate and deposited material have the same crystqal structure
    
  /* output initial simulation data */
  printf(" output initial state of the simulation box............Done\n");
  //printf("----------------------------------\n \n");
  SimulationBoxOutput(0);
  //cout << "ERROR_MAIN";
  NeighborsListOutput();
    
  /* initialisation of random number generator */
  srand(time(NULL)); 
               
  /* generation of the surface map */
  printf(" generation of the surface map.........................Done\n");
  //printf("----------------------------------\n \n");
  SurfaceMapGeneration(1);  
  
  /* calculation of the initial configuration of the box */
  //cout << "ExistingConfigurationNumber=" << ExistingConfigurationNumber << endl;
  
  /* system evolution */
  /* 2D deposition of the first atom */
  printf(" deposition of the first atom..........................Done\n");
  //printf("----------------------------------\n \n");
  RandNumber1 = 1 + rand()% SurfaceMapLength;    //random atom position on the surface
  printf(" calculation of the initial configuration of the box...Done\n");
  //printf("----------------------------------\n \n");
  LastAtomCalculation(RandNumber1, 1); 
  BoxInitialConfiguration();
  //RandNumber2 = 1 + rand()% DepositedAtomNumber; // random atom type 
  DepositionNearest2D( SurfaceMap[RandNumber1], RandNumber1, 2);
  fclose(FService1);
  ConfigurationOutput("Config1");
  SystemStateOutput("System1");
  //AtomConfigurationNumberOutput("Dep_Check1");
  //return 1;
  //return 1;
  //update system information
  CurrentStep = 1;
  CurrentAtom = 1;
  printf(" initial values of process rates calculation...........Done\n");
  printf("----------------------------------\n \n");
  InitialRatesCalculation( 1 );
  RatesOutput("Rates1");
  //getch();
  //BoxConfigurationUpdate( -1, SurfaceMap[RandNumber1], 1 ); 
  printf("      System evolution calculation\n");
  /* calculate system evolution until the endof the time/time steps/deposited atoms */
  if(ProgrammCheck == 1)
  {
  	while(CurrentTime != TotalTime)
  	{
  	  SystemEvolution();
  	  CurrentTime = -1;
  	}
  }
  if(ProgrammCheck == 2)
  {
  	while(CurrentStep != TotalSteps)
	{
  	  
  	  //getch();
	  SystemEvolution();
  	  CurrentStep = CurrentStep + 1;
  	}
  }
  if(ProgrammCheck == 3)
  {
  	while(CurrentTime == TotalAtoms)
	{
  	  SystemEvolution();
  	  CurrentAtom = -1;
  	}
  }
  
  
  SimulationBoxOutput(1);
  /* cleaning memory from the simulation data */
  MemoryCleaning();
  
  if(FileClose1 == 0)
    fclose(FError);
  
  if(FileClose2 == 0)  
    fclose(FExec);
   
}
