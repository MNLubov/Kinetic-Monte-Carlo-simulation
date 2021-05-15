#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <limits>

#include "array.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"

using namespace std;

//function outputs Simulation box in the file
void SimulationBoxOutput( int BeginEndIndex )
{
   int i, j, k;
   int XLength, ZLength;
   
   std::fstream file;
   std::fstream file2;
   std::string string2, string3, string4, string5;
   std::ostringstream string1;
   FILE *F1, *F2, *F3, *F4, *F5;
   std::fstream OutFile;
   
   
   //file.open("SimulBox",std::ios::out);
   //file2.open("SimulBoxCoordinates",std::ios::out);
   if(BeginEndIndex == 0)
   {
     F1 = fopen("SimulBoxInitial", "w+");
     F2 = fopen("SimulBoxCoordinatesInitial", "w+");
    
     for(i = 1; i <= SimulationBoxAtomsNumber - 1; i++)  
     {
       j = i * 3;
       fprintf(F1, "%d  %d\n", i, SimulationBox[i]);
       fprintf(F2, "%d  %f  %f  %f\n",SimulationBox[i], SimulationBoxCoordinates[j - 2], SimulationBoxCoordinates[j - 1], SimulationBoxCoordinates[j]);
     } 
   
     fclose(F1);
     fclose(F2);
     
     //cout << "ERRORRR!!!!!";
     //getch();
   //F5 = fopen("SimulationBoxNumbers","w+");
     OutFile.open("SimulationBoxNumbers", std::ios::out);
   
     for(i = 2 * BoxSizeZ; i >= 1; i--)
     {
       for(j = 1; j <= BoxSizeX; j++)
       {
         if((i % 2) == 1)
     	 {
     	   k = (i - 1) * BoxSizeX + 2 * j - 1;
     	   //cout << "ERROR k= " << k ;
     	   //getch();
     	   if(k < 10)
     	   {
     	    string1 << k;
     	  	string2 = "00";
     	  	string4 = string1.str();
			string3 = string2 + string4;		
			//cout << string3 << endl;
			//string3.clear();
			//cout << k ;
			//getch();
     	  	OutFile << string3 << "  ";	
     	  	string1.str("");
			string1.clear();
     	   }
     	   if(k < 100 && k >= 10)
     	   {
     	   	string1 << k;
     	  	string2 = "0";
     	  	string4 = string1.str();
     	  	string3 = string2 + string4;
     	  	OutFile << string3 << "  ";	
     	  	string1.str("");
			string1.clear();
     	   }
     	   if(k < 1000 && k >= 100)
     	   {
     	  	 string1 << k;
     	  	 string2 = string1.str();
     	  	 string3 = string2;
     	  	 OutFile << string3 << "  ";
		 	 string1.str("");
			 string1.clear();	
     	   }
		  
     	 }
		 if((i % 2) == 0)
     	 {
     	  k = (i - 2) * BoxSizeX + 2 * j;
		  if(k < 10)
     	  {
     	  	string1 << k;
     	  	string2 = "00";
     	  	string4 = string1.str();
     	  	string3 = string2 + string4;     	  	
     	  	//cout << string3 <<" ERROR";
     	  	//getch();
			OutFile << "  " << string3;	
			string1.str("");
			string1.clear();
     	  }
     	  if(k < 100 && k >= 10)
     	  {
     	  	string1 << k;
     	  	string2 = "0";
     	  	string4 = string1.str();
     	  	string3 = string2 + string4;     	  	
     	  	OutFile << "  " << string3;	
     	  	string1.str("");
			string1.clear();
     	  }
     	  if(k < 1000 && k >= 100)
     	  {     	  	
     	  	//string1.clear();
			string1 << k;
     	  	string2 = string1.str(); 
     	  	
			string3 = string2;
			//string3 = "1";
			//cout << string1 << endl;
     	  	//cout << "string3=" << string3 << endl;
			//cout << "string2=" << string2 << endl;   
			//getch();  	  	
     	  	OutFile << "  " << string3;	     	  	
			string1.str("");
			string1.clear();
     	  	//string2.clear();
     	  	//string3.clear();
     	  }
		  	
         }
       }
       OutFile << endl;
     }
     OutFile.close();
     
   }
   
   if(BeginEndIndex == 1)
   {
     F3 = fopen("SimulBoxFinal", "w+");
     F4 = fopen("SimulBoxCoordinatesFinal", "w+");
    
     for(i = 1; i <= SimulationBoxAtomsNumber - 1; i++)  
     {
       j = i * 3;
       fprintf(F3, "%d  %d\n", i, SimulationBox[i]);
       fprintf(F4, "%d  %f  %f  %f\n",SimulationBox[i], SimulationBoxCoordinates[j - 2], SimulationBoxCoordinates[j - 1], SimulationBoxCoordinates[j]);
     } 
   
     fclose(F3);
     fclose(F4);
   }
}

void NeighborsListOutput( void )
{
   int i, j, k;
   
   std::fstream file;
   std::fstream file2;
   FILE *F1, *F2;
   
   F1 = fopen("NearestNeighbors", "w+");
   F2 = fopen("NextNearestNeighbors", "w+");

   for(i = 1; i <= SimulationBoxAtomsNumber - 1; i++)  
   {
     j = i * 4;
     fprintf(F1, "%d  %d  %d  %d  %d\n",i, NearestNeighborList[j - 3], NearestNeighborList[j - 2], NearestNeighborList[j - 1], NearestNeighborList[j]);
     fprintf(F2, "%d  %d  %d  %d  %d\n",i, NextNearestNeighborList[j - 3], NextNearestNeighborList[j - 2], NextNearestNeighborList[j - 1], NextNearestNeighborList[j]);   
     //printf("%d  %d  %d  %d  %d\n",i, NearestNeighborList[j - 3], NearestNeighborList[j - 2], NearestNeighborList[j - 1], NearestNeighborList[j]);
     //printf("%d  %d  %d  %d  %d\n",i, NextNearestNeighborList[j - 3], NextNearestNeighborList[j - 2], NextNearestNeighborList[j - 1], NextNearestNeighborList[j]);
     //getch();
   }
   
   fclose(F1);
   fclose(F2);
     
}


void ConfigurationOutput( string ConfigurationName ) 
{
  int i, j, k;	
  double En;
  std::fstream ConfigOut;
  
  ConfigOut.open(ConfigurationName.c_str(), std::ios::out); 
  ConfigOut.precision(4);
  ConfigOut.setf(ios::fixed);
  ConfigOut.setf(ios::showpoint);
  for(i = 1; i <= ExistingConfigurationNumber; i++)
  {
  	ConfigOut << "Energy = " << EnergyList[i] << "   Atoms in Configuration = " << SimulationBoxConfiguration[i][0] << endl;
  	for(j = 1; j <= SimulationBoxConfiguration[i][0]; j++)
  	{
  	  //En = AtomSiteEnergy(SimulationBoxConfiguration[i][j]);
	  k = SimulationBoxConfiguration[i][j] * 4; 
	  ConfigOut << j <<  "  " << SimulationBoxConfiguration[i][j] << "  " <<
	  NearestNeighborList[k - 3] << "  " << NearestNeighborList[k - 2] << "  " << 
	  NearestNeighborList[k - 1] << "  " <<  NearestNeighborList[k] <<  endl;
  	//getch();
    }
  }
 
  ConfigOut.close(); 
}

void SystemStateOutput( string SystemOutName )
{
  int i, j, k;	
  std::fstream SystemOut;
  
  //cout << "SYSTEM STATE OUTPUT1" << endl;
  SystemOut.open(SystemOutName.c_str(), std::ios::out); 
  for(i = 1; i <= SimulationBoxAtomsNumber - 1; i++)  
  {
  	j = i * 3;
	SystemOut << SimulationBox[i] << "   " << SimulationBoxCoordinates[j - 2]
	          <<  "  " << SimulationBoxCoordinates[j - 1] << "   " <<  SimulationBoxCoordinates[j] << endl;
  	
  }
 
  //cout << "SYSTEM STATE OUTPUT2" << endl;
  SystemOut.close(); 
}

void RatesOutput( string RatesOutName )
{
	int i, j, k;
	int ProcNumb;
	std::fstream RatesOut;
	
	RatesOut.open(RatesOutName.c_str(), std::ios::out); 
	ProcNumb = AtomTypesNumber + ExistingConfigurationNumber;
	//cout << ProcNumb;
	//getch();
    for(i = 1; i <= ProcNumb; i++) 
    {
       RatesOut << i <<  "  " << ProcessRate[i] << endl;
    }
 
    RatesOut.close(); 
}

void SurfaceMapOutput( string SurfaceMapOutName )
{
	int i, j, k;
	int ProcNumb;
	std::fstream SurfaceMapOut;
    std::string str1;
	//cout << "Map output 1! " << SurfaceMapLength << endl;
	//getch();
	SurfaceMapOut.open(SurfaceMapOutName.c_str(), std::ios::out); 
	//SurfaceMapOut.open("Proverka", std::ios::out); 
	//ProcNumb = AtomTypesNumber + ExistingConfigurationNumber;
	//cout << "Map output 2! " << SurfaceMapLength << endl;
	//getch();
    for(i = 1; i <= SurfaceMapLength; i++) 
    {
       SurfaceMapOut << i <<  "  " << SurfaceMap[i] << endl;
       //cout << i <<  "  " << SurfaceMap[i] << endl;
       //getch();
    }
        
    SurfaceMapOut.close(); 
    
    //SurfaceMapOutName.str("");
	//SurfaceMapOutName.clear();
	
}


void AtomConfigurationNumberOutput( string AtomConfOutName )
{
	int i, j, k;
	int ProcNumb;
	std::fstream AtomConfOut;
	
	AtomConfOut.open(AtomConfOutName.c_str(), std::ios::out); 
	//ProcNumb = AtomTypesNumber + ExistingConfigurationNumber;
	//cout << SimulationBoxAtomsNumber << endl;
	//getch();
    for(i = 1; i <= 2 * (SimulationBoxAtomsNumber - SubstrateBottomNumber); i++)
    {
       if(i % 2 == 0)
	     AtomConfOut << i << "  " << (i / 2 + SubstrateBottomNumber) << "  " << AtomConfigurationNumber[i] <<  endl;
	   //if(i % 2 == 1)  
	     //AtomConfOut << i << "  " << AtomConfigurationNumber[i] <<  endl;
    }
 
    AtomConfOut.close(); 
    //cout << "config output 2" << endl;
    
}

