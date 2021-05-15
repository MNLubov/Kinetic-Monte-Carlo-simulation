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

#include "parameters.h"


#define pi 3.1415926

using namespace std;

int temp;

         /* Variable names for unit cell */
// UnitCellAtomsNumber - number of the atoms in the unit cells
// UnitCellSizeX - translation distance along X direction      
// UnitCellSizeY - translation distance along Y direction
// UnitCellSizeZ - translation distance along Z direction 
// UnitCellCoordinates[300] - coordinates of unit cell atoms
// UnitCellAtomsType[100] - type of atoms in unit cell
// UnitCellAtomsTypeNumber - number of atoms of different type in the unit cell
// UnitCellTermination - index of the termination of the unit cell (simulation box), UnitCellTermination = 1,2 - first and second elements in chemical formula, respectively
// CrystalStructureType - crystal structure of the unit cell
// UnitCellLengthX, UnitCellLengthY, UnitCellLengthZ - length of the unit cell along X,Y,Z directions respectively
// NearestNeighbornumber - number of the nearest neighbors of the atom 
// NextNearestNeighbornumber - number of the next nearest neighbors of the atom 
// LatticeParameter - lattice parameter of the deposited material
// SubstrateLatticeParameter - lattice parameter of the substrate
int UnitCellAtomsNumber;
double UnitCellSizeX, UnitCellSizeY, UnitCellSizeZ;
double UnitCellCoordinates[101][4];
int UnitCellAtomsType[101];
int UnitCellAtomsTypeNumber;
int UnitCellTermination;
int CrystalStructureType;
double UnitCellLengthX, UnitCellLengthY, UnitCellLengthZ;
int NearestNeighbornumber;
int NextNearestNeighbornumber;
double LatticeParameter;
double SubstrateLatticeParameter;


       /* Variable names for simulation box */
// BoxSizeX - size of the simulation box along X direction in unit cells (number of unit cells)
// BoxSizeY - size of the simulation box along Y direction in unit cells
// BoxSizeZ - size of the simulation box along Z direction in unit cells
// SubstrateSizeZ - size of the substrate box along Z direction in unit cells
// SubstrateAtomsNumber - total number of atoms in the susbtrate box
// SimulationBoxAtomsNumber - total number of atoms in the simulation box (including free sites)
// SubstrateType - crystal type of the substrate
// SubstrateLayers - number of layers in the substrate
// SubstrateTermination - index of the termination of the substrate
// Dimension - dimension of the substrate
// DepostionMapLength - length of the SurfaceMap[] array 
// AtomTypesNumber - total amount of the atoms in the system
// DepositedAtomNumber - number of the atom types which is deposited on the surface
// DepositedAtomType - array of atoms types deposited on the substrate = 1 - deposited, = 0 - not deposited
// SiteStateBefore, SiteStateAfter - state of the site (empty, occupied) before and after the process occured
int BoxSizeX, BoxSizeY, BoxSizeZ;
int SubstrateSizeZ;
int SubstrateAtomsNumber;
int SimulationBoxAtomsNumber;
int Dimension;
int SubstrateType;
int SubstrateLayers;
int SubstrateTermination;
int SurfaceMapLength;
int AtomTypesNumber;
int DepositedAtomNumber;
int DepositedAtomType[11];
int SiteStateBefore, SiteStateAfter;

            /* Variables name for the arrays*/
// SimulationBox - array of free/occupied sites in the simulation box 
// SimulationBoxCoordinates - array of sites coordinates in the simulation box 
// NearestNeighbourList - list of the nearest neighbours in the simulation box
// NextNearestNeighbourList - list of the next nearest neighbours in the simulation box
// SurfaceMap - deposition positions for the atoms
// AtomSitesConfiguration - array of the atom sites energies
// SimulationBoxConfiguration - vector of vectors for storage of the different configurations existing in the system 
//                             (atom with 1 nearets neighbor, atom with 2 nearest neighbors etc)
// EnergyList - vector(array) of energies of all existing configurations in the system
// NearestNeighborsMax - vector(array) of numbers of the atoms with maximal value of the nearest neighbors 
// AtomConfigurationNumber - position of the atom in the vector of vectors
// ProcessRate - array of rates for all physical processes
// DiffusionRate - vector of rates for the diffusion processes. Array have the same size as 'EnergyList' array
int *SimulationBox;
double *SimulationBoxCoordinates; 
int *NearestNeighborList;
int *NextNearestNeighborList;
int *SurfaceMap; 
double *AtomSitesConfiguration;
vector< vector<int> > SimulationBoxConfiguration;
vector<double> EnergyList;
vector<int> NearestNeighborsMax; 
int *AtomConfigurationNumber; 
vector<double> ProcessRate;
vector<double> DiffusionRate;
vector<double> DebayeFrequencies;
//vector<int> GlobalProcessIndex;
//int *SimulationBoxConfiguration;
//int *DifferentAtomsConfigurations;


         /* Variables for the processes rate callculation */
// LastAtomNumber - number of the last atom in the simulation box
// BondEnergy[11][21]- bond energies between atoms, [1-10][1-10] - nearest neighbors; [1-10][11-20] - next nearest neighbors
// BoxConfigurationNumber - number of different atoms configurations (different neighborhoods) in the simulation box
// ExistingConfigurationNumber - number of configurations that exist in the simulation box
// SubstrateBottomNumber - The highest number of the atom in the substrate part where proccesses not allowed
// SubstrateProcessSize - Height (in unit cells, from the top of the substrate) of the substrate part where processes are allowed
// TotalProcessNumber - number of the processes in the system on current programm step
// TotalProcessNumberPrevios - number of the processes in the system on the previous programm step
// DepositionTotalRate - summ of all deposition rates
int LastAtomNumber;               
double BondEnergy[11][21];          
int BoxConfigurationNumber;       
int ExistingConfigurationNumber;  
int SubstrateBottomNumber;  
int SubstrateProcessSize;  
int TotalProcessNumber;
int TotalProcessNumberPrevious;
double DepositionTotalRate;


//
// BoxSizeLX - length of the simulation box along X direction in Angstrems
// BoxSizeLY - length of the simulation box along Y direction in Angstrems
// BoxSizeLZ - length of the simulation box along Z direction in Angstrems
//
//
// SitesNeighbours - array of neighhbours sites numbers of the simulation box sites
// OccupiedSites - array of indicators of empty and occupied sites in the simulation box
//
// 

       /* Variable names for physical parameters */
// Flux - deposition flux on the surface in 1/(cm2sec1) for different types of atoms
// Temperature - substrate temperature in K
// HopEnergy - energy of adatom hop
// EvaporationEnergy - additional energy/energy of the evaporation
// NearestBondEnergy - energy per 1 bond with the nearest neighbor
// NextNearestBondEnergy - energy per 1 bond with the next nearest neighbor  
// AtomFrequency - debye frequency of the atom
double Flux[11];
double Temperature;
double HopEnergy[11];
double EvaporationEnergy[11];
double NearestBondEnergy[100];
double NextNearestBondEnergy[100];
double AtomFrequency[100]; 

/* variable names for the time steps */
// ProgrammCheck - indication of the programm evaluation (in secs/programm steps/deposited atoms)
// TotalTime - total time in secs
// TotalSteps - total number of programm steps
// TotalAtoms - total number of the deposited atoms 
int ProgrammCheck;
int TotalTime;
int TotalSteps;
int TotalAtoms;


     /*  Variable names for physical processes   */
//  Probabilities - realization probabilities of the physical procceses
//  Rates - realization rates of the physical procceses

				/* System flags */
// PhysicalProcessesIndex - indicates whether specific physical process can occur in the system 
// PhysicalProcessesIndex[] = 1 can occur, PhysicalProcessesIndex[] = 0 - can't
// PhysicalProcessesIndex[1] - deposition
// PhysicalProcessesIndex[2] - evaporation
// PhysicalProcessesIndex[3] - surface diffusion
// PhysicalProcessesIndex[4] - volume diffusion
//
int PhysicalProcessIndex[10];

/* variables names for the log files */
// FError - file for the errors log
// FExec - file for the execution data log
// FService1 - ancilary file
// FCheck - file for checking differnet things
FILE *FError;
FILE *FExec;
FILE *FService1;
FILE *FCheck, *FCheck2, *FCheck3;

/* program additional parameters */
double deltaE;

