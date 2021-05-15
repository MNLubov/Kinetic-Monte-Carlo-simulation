#ifndef _H_PARAMETERS_H_
#define _H_PARAMETERS_H_

#include <istream>
#include <fstream> 
#include <ios>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

extern int temp;

extern int UnitCellAtomsNumber;
extern double UnitCellSizeX, UnitCellSizeY, UnitCellSizeZ;
extern int UnitCellAtomsType[101];
extern int UnitCellAtomsTypeNumber;
extern double UnitCellCoordinates[101][4];
extern int CrystalStructureType;
extern int UnitCellTermination;
extern double UnitCellLengthX, UnitCellLengthY, UnitCellLengthZ;
extern int NearestNeighbornumber;
extern int NextNearestNeighbornumber;
extern double LatticeParameter;
extern double SubstrateLatticeParameter;

extern int BoxSizeX, BoxSizeY, BoxSizeZ;
extern int SubstrateSizeZ;
extern int SubstrateAtomsNumber;
extern int SimulationBoxAtomsNumber;
extern int Dimension;
extern int SubstrateType;
extern int SubstrateLayers;
extern int SubstrateTermination;
extern int SurfaceMapLength;
extern int AtomTypesNumber;
extern int DepositedAtomNumber;
extern int DepositedAtomType[11];
extern int SiteStateBefore, SiteStateAfter;

extern int *SimulationBox;
extern double *SimulationBoxCoordinates; 
extern int *NearestNeighborList; 
extern int *NextNearestNeighborList;
extern int *SurfaceMap;
extern double *AtomSitesConfiguration;
extern int *AtomConfigurationNumber; 
extern vector< vector<int> > SimulationBoxConfiguration;
extern vector <double> EnergyList;
extern vector<int> NearestNeighborsMax;
extern int *AtomConfigurationNumber; 
extern vector<double> ProcessRate;
extern vector<double> DiffusionRate;
extern vector<double> DebayeFrequencies;

extern double Flux[11];
extern double Temperature;
extern double HopEnergy[11];
extern double EvaporationEnergy[11];
extern double NearestBondEnergy[100];
extern double NextNearestBondEnergy[100];
extern double AtomFrequency[100];

extern int ProgrammCheck;
extern int TotalTime;
extern int TotalSteps;
extern int TotalAtoms;

extern int LastAtomNumber;
extern double BondEnergy[11][21];
extern int BoxConfigurationNumber;
extern int ExistingConfigurationNumber;
extern int SubstrateBottomNumber;  
extern int SubstrateProcessSize;  
extern int TotalProcessNumber;
extern int TotalProcessNumberPrevious; 
extern double DepositionTotalRate;


extern int PhysicalProcessIndex[10];


extern FILE *FError;
extern FILE *FExec;
extern FILE *FService1;
extern FILE *FCheck, *FCheck2, *FCheck3;

extern double deltaE;

#endif
