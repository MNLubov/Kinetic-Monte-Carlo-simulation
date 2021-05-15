#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string.h>
#include <iostream>
#include <vector>

#include "parameters.h"
#include "array.h"
#include "init.h"
#include "unit_cell.h"

using namespace std;

int ReadData( void )
{
  int i, j, k;
  
  std::fstream InitFile, InitFile2, InitFile3;
  std::string str;
  int InitIndex;
  int DepositedAtomNumber;
  int Type;  // type of the atom
  double Freq; // Deby frequency, variable is not used in the programm
  int Temp1;

  /*
  InitFile.open("InputData",std::ios::in);
  getline(InitFile,str);
  InitFile >> BoxSizeX >> BoxSizeY >> BoxSizeZ >> SubstrateSizeZ; 
  getline(InitFile,str);
  InitFile >> Dimension >> SubstrateType >> UnitCellAtomsTypeNumber >>  SubstrateTermination;
  printf("%d  %d  %d  %d\n", Dimension, SubstrateType,  UnitCellAtomsTypeNumber, SubstrateTermination);
  //getch();
  
  InitFile.close();
  */
  //printf("D\n");
  //getch();
  fprintf(FError, "Initialisation error log\n");
  fprintf(FExec, "Initialisation log\n");
  InitIndex = 0;
  
  InitFile2.open("InputData", std::ios::in);
  //read parameters of the substrate 
  getline(InitFile2, str);
  InitFile2 >> str >> Dimension;
  if(Dimension < 2 || Dimension > 3 )
  {
    fprintf(FError, "Dimension of the simulation box could only be '2' or '3'\n");
    InitIndex = 1;
  }
  //printf("dimension = %d",Dimension);
  InitFile2 >> str >> SubstrateType;
  InitFile2 >> str >> UnitCellAtomsTypeNumber;
  InitFile2 >> str >> SubstrateTermination;
  InitFile2 >> str >> SubstrateLatticeParameter;
  //read size of the substrate and the box and number of diiferent atom types in the system
  getline(InitFile2, str);
  //cout << str << endl;
  //getch();
  getline(InitFile2, str);
  //cout << str << endl;
  //getch();
  InitFile2 >> str >> BoxSizeX;
  //cout << BoxSizeX << endl;
  //getch();
  InitFile2 >> str >> BoxSizeY;
  InitFile2 >> str >> BoxSizeZ;
  InitFile2 >> str >> SubstrateSizeZ;
  InitFile2 >> str >> AtomTypesNumber;
  if(SubstrateSizeZ >= BoxSizeZ)
  {
    fprintf(FError, "Substrate Size exceeds Simulation Box size\n");
    fprintf(FError, "SimulationBoxSize alon Z = %f   Substrate size along Z = %f\n",BoxSizeZ, SubstrateSizeZ);
    InitIndex = 1;
  }
  //read data of the deposited material parameters
  //getline(InitFile2, str);
  getline(InitFile2, str);
  getline(InitFile2, str);
  InitFile2 >> str >> LatticeParameter;
  //cout << LatticeParameter << endl;
  //getch();
  //read deposited atoms types
  InitFile2 >> str >> DepositedAtomNumber;
  //cout << "DepositedAtomNumber=" << DepositedAtomNumber;
  //getch();
  if(DepositedAtomNumber > AtomTypesNumber)
  {
    fprintf(FError, "Number of the types of the deposited atoms is more than types of the atoms in the system\n");
    InitIndex = 1;
  }
  for(i = 1; i <= AtomTypesNumber; i++)
    InitFile2 >> str >> Flux[i];
  InitFile2 >> str >> ProgrammCheck;
  if(ProgrammCheck == 1)
  {
  	InitFile2 >> str >> TotalTime;
	InitFile2 >> str >> Temp1;
  	InitFile2 >> str >> Temp1;
  }
  if(ProgrammCheck == 2)
  {
  	InitFile2 >> str >> Temp1;
  	InitFile2 >> str >> TotalSteps;
  	InitFile2 >> str >> Temp1;
  }
  if(ProgrammCheck == 3)
  {
  	InitFile2 >> str >> Temp1;
  	InitFile2 >> str >> Temp1;
  	InitFile2 >> str >> TotalAtoms;
  }
  getline(InitFile2, str);
  getline(InitFile2, str);
  InitFile2 >> str >> Temperature;
  //cout << str << "  !  " << Temperature;
  //getch();
  
  //cout << Flux[2] << endl;
  //getch();
    
  //printf("%d  %d  %d  %d\n", Dimension, SubstrateType,  UnitCellAtomsTypeNumber, SubstrateTermination);
  //getch();
  InitFile3.open("EnergyDebye", std::ios::in);
  getline(InitFile3, str);
  //cout << str << endl;
  //getch();
  getline(InitFile3, str);
  //cout << str << endl;
  //getch();
  for(i = 1; i <= AtomTypesNumber; i++)
    InitFile3 >> Type >> BondEnergy[i][1] >> BondEnergy[i][2] >> BondEnergy[i][3] >> BondEnergy[i][4] >>
                         BondEnergy[i][5] >> BondEnergy[i][6] >> BondEnergy[i][7] >> BondEnergy[i][8] >>
                         BondEnergy[i][9] >> BondEnergy[i][10] >> AtomFrequency[i];
  //cout << Type << " "<< BondEnergy[1][1] << endl;
  //getch();
  //cout << "Error" ;
  //getch();
  getline(InitFile3, str);
  getline(InitFile3, str);
  //cout << str << endl;
  //getch();
  getline(InitFile3, str);
  //cout << str << endl;
  //getch();
  for(i = 1; i <= AtomTypesNumber; i++)
    InitFile3 >> Type >> BondEnergy[i][11] >> BondEnergy[i][12] >> BondEnergy[i][13] >> BondEnergy[i][14] >>
                         BondEnergy[i][15] >> BondEnergy[i][16] >> BondEnergy[i][17] >> BondEnergy[i][18] >>
                         BondEnergy[i][19] >> BondEnergy[i][20] >> Freq;
  
  //cout << Type << " "<< BondEnergy[1][11] << endl;
  //getch();
  
  InitFile2.close();
  InitFile3.close(); 
  
  /*
  for(i = 1; i <=AtomTypesNumber; i++)
    for(j = 1; j <= 20; j++)
    {
      cout << "BondEnergy[" << i << "]" << "[" << j << "]=" <<BondEnergy[i][j] << endl; 
      //getch();	
    }  
  */
    
  if(InitIndex == 0)
  {
  	fprintf(FError, "Initialisation sucsesfully done!\n");
  	return 1;
  }
  if(InitIndex == 1)
    return 2;                    
}

void RandomNumberCheck( void )
{
  int i, j, diap;
  int rnm[1000];
  double j0;
  FILE *F;
  
  for(i = 0; i <= 1000; i++)
    rnm[i] = 0; 
  
  j = 100000000;
  diap = 1000;
  j0 = (double)(j / diap);
  while(j >= 0)
  {
    i = 1 + rand()%RAND_MAX;
    rnm[i] = rnm[i] + 1;
    j = j - 1;
    //printf("%d\n",j);
  }
  
  F = fopen("rnmcheck","w+");
  for(i = 0; i <= 1000; i++)
    fprintf(F,"%d  %f\n",i, rnm[i] / j0);
  
  fclose(F);  
    
  
}

int factorial( int n )
{
  if(n == 0)
    return 1;
    
  if(n > 1)
    return (n * factorial(n - 1)); 
}

void PreDeposition ( void )
{
  int i, j;
  
  //SimulationBox = (int *) malloc(SimulationBoxSize * sizeof(int));
  //set initial parameters of the system
  CrystalStructureType = SubstrateType;
  UnitCellTermination = SubstrateTermination;
  
  //set indexes
  PhysicalProcessIndex[4] = 0;
  
  if(Dimension == 3)
    ;
  if(Dimension == 2)
  {
    UnitCellGeneration2D(CrystalStructureType, UnitCellAtomsTypeNumber, UnitCellTermination); 	
  }
  if(CrystalStructureType == 5 && SubstrateType == 5)
  {
  	NearestNeighbornumber = 4;
    NextNearestNeighbornumber = 4;
  }
  
  if(PhysicalProcessIndex[4] == 0 && SubstrateSizeZ >= 3)
    SubstrateProcessSize = SubstrateSizeZ - 2;
  else
  {
  	printf("Error, substrate size Z must be greater than 3 UnitCells "); 
  	fprintf(FError, "Error, substrate size Z must be greater than 3 UnitCells\n");
  	fprintf(FError, "SubstrateProcessSize = %d  SubstrateSizeZ = %d\n",SubstrateProcessSize, SubstrateSizeZ);
  }  
  SubstrateAtomsNumber = UnitCellAtomsNumber * BoxSizeX * SubstrateSizeZ + 1;
  SimulationBoxAtomsNumber = UnitCellAtomsNumber * BoxSizeX * BoxSizeZ + 1;  
  /*
  cout << "SimulationBoxAtomsNumber=" <<SimulationBoxAtomsNumber << endl; 
  cout << "UnitCellAtomsNumber=" << UnitCellAtomsNumber << endl;
  cout << "BoxSizeX=" << BoxSizeX << endl;
  cout << "BoxSizeZ=" << BoxSizeZ << endl;
  //getch();
  */
  SubstrateBottomNumber = (SubstrateAtomsNumber - 1) - (BoxSizeX * UnitCellAtomsNumber) * (SubstrateSizeZ - SubstrateProcessSize);
  LastAtomNumber = SubstrateAtomsNumber - 1;
  fprintf(FExec,"SubstrateAtomsNumber = %d\n", SubstrateAtomsNumber);
  fprintf(FExec,"SimulationBoxAtomsNumber = %d\n", SimulationBoxAtomsNumber);
  fprintf(FExec,"SubstrateBottomNumber = %d\n", SubstrateBottomNumber);
  fprintf(FExec,"LastAtomNumber = %d\n", LastAtomNumber);
  //cout << "LastAtomNumber =" << LastAtomNumber << endl;
  //getch();
  
  
  BoxConfigurationNumber = factorial(AtomTypesNumber * 2 + (NearestNeighbornumber + NextNearestNeighbornumber) - 1) / 
                           (factorial(NearestNeighbornumber + NextNearestNeighbornumber) * factorial(AtomTypesNumber * 2 - 1));
                           
  
  SimulationBoxConfigurationMemoryAllocation();
  std::fstream FILE1;
  FILE1.open("AtomConfNumber",std::ios::out);
  //allocate memory for the AtomConfigurationNumber[] array
  //SimulationBoxConfigurationMemoryAllocation();
  for(i = 1; i <= 2 * (SimulationBoxAtomsNumber - SubstrateBottomNumber); i++)
  {
  	 AtomConfigurationNumber[i] = 0;
  	 FILE1 << i << "  " << AtomConfigurationNumber[i] <<  endl;
  }
  FILE1.close(); 
  
  //setting additional parameters
  deltaE = 0.00001; 
}


void SetFrequencies ( int FirstNewProcess, int LastNewProcess, int ConfigurationFrequency )
{
   int i;
   int FirstAtom, AtomType;
   
  
  DebayeFrequencies.push_back(int());
  
  //set Debye frequency for the configuration
  if(ConfigurationFrequency == 1)
  {
    for(i = FirstNewProcess; i <= LastNewProcess; i++)
    {
      FirstAtom = SimulationBoxConfiguration[i][1];
      AtomType = SimulationBox[FirstAtom];
      DebayeFrequencies.push_back(AtomFrequency[AtomType]);     
    }  
  }
    
}
