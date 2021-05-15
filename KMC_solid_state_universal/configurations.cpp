#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string.h>
#include <iostream>
#include <vector>

#include "array.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"

#define kT0 0.025 // im meV

using namespace std;

void BoxInitialConfiguration( void )
{
  int i, j, k;
  int FirstAtomIndex;
  int CurrentAtomEnergy;
  int ExistConfigurationIndex;     // indication TRUE or FALSE that current energy of atom site corresponds with existing configuration
  int ConfigurationSiteNumber;     // number of the sites in the configuration 
  FILE *F1;
    
  /* Substrate
  ---------------------
  |  Processes allowed |
  |                    |
  ---------------------* <-- SubstrateBottomNumber
  |     Processes      |
  |    not allowed     |
  ---------------------
  */
  
  //SubstrateProcessSize;    // Height (in unit cells, from the top of the substrate) of the substrate part where processes are allowed
  if(SubstrateProcessSize >= SubstrateSizeZ)
  {
    printf("Substrate Part where processes are allowed is equal to the Substrate Size\n");
    //getch();
    //return 0;
  }  
  
  k = 1;
  //cout << "SubstrateBottomNumber=" << SubstrateBottomNumber << endl;
  //cout << "LastAtomNumber=" << LastAtomNumber << endl;
  //getch();
  
  ExistingConfigurationNumber = 0; // set initial values of the configurations in the system
  //printf("SubstrateBottomNumber = %d\n",SubstrateBottomNumber);
  //printf("SubstrateProcessSize = %d\n",SubstrateProcessSize);
  //getch();
  FirstAtomIndex = 1;
  for(i = SubstrateBottomNumber + 1; i <= LastAtomNumber; i++)
  {
	 //calculation of the atom enrgy in the site of the simulation box
     AtomSitesConfiguration[i] = AtomSiteEnergy(i);
     //AtomConfigurationNumberOutput("Dep_Check001Pre");
	 //cout << "i = " << i << endl;
     //getch();
     //atom has 4 nearest neighbors
     /*
     if(AtomSitesConfiguration[i] == -1)
     {
       if(k == 1)
       {
         NearestNeighborsMax.push_back(int());
         NearestNeighborsMax.push_back(i);
         k++;
       }  
       else
       {
         NearestNeighborsMax.push_back(i);
         k++;
       }  
     }
     */
     
     //excluding atoms with full set of the nearest neighbors
     if(AtomSitesConfiguration[i] != -1)
     {     
       //increase number of rows in the SimulationBoxConfiguration vector of vectors by 2. SimulationBoxConfiguration[0] is not considered further
       if(ExistingConfigurationNumber == 0)
       {      
         //printf("ExistingConfigurationNumber = 0  i = %d\n",i);
         //printf("%d  %d  %d  %d\n", NearestNeighborList[i * 4 - 3], NearestNeighborList[i * 4 - 2], NearestNeighborList[i * 4 - 1], NearestNeighborList[i * 4]);
         //getch();
         FirstAtomIndex = i;
		 //create vectort number '0' in the configuration vector of vectors (not used)
  		 SimulationBoxConfiguration.push_back(vector<int>());
  		 //create first vector for the first configuration
         SimulationBoxConfiguration.push_back(vector<int>());  		
   		 //create element number '0' in the vector of energies (not used)
  		 EnergyList.push_back(int()); 
		 //put the energy of the first considered atom in the element of EnergyList[1]
         EnergyList.push_back(AtomSitesConfiguration[i]);
         //cout << EnergyList[1] << AtomSitesConfiguration[i];
         //getch();
         //put energy of the atom configuration by the number of the element in the EnergyList array 
         AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2 - 1] = 1;   //first existing configuration  
                  
         //setting number of the atoms in the first configuration to the zero element of vector in the vector of vectors
         SimulationBoxConfiguration[1].push_back(1);
         //setting number of the atom in the first configuration
         SimulationBoxConfiguration[1].push_back(i);

         //write number of the element in the SimulationBoxConfiguration vector with energy EnergyList[1]      
         AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2] = 1;
                  
         //set number of the configurations in the system equal to 1
         ExistingConfigurationNumber = 1;
         //cout << "ExistingConfigurationNumber=" << ExistingConfigurationNumber << endl;
         //getch();
         //FirstAtomIndex = 0;
         
       }
     
       //add all other element to the existing configuration or create new one
       if(ExistingConfigurationNumber > 0 && i != FirstAtomIndex)
       {
         ExistConfigurationIndex = 0; 
         //seek if appropriate configuration exists
         for(j = 1; j <= ExistingConfigurationNumber; j++)
         {
           if(AtomSitesConfiguration[i] > (EnergyList[j] - deltaE) && AtomSitesConfiguration[i] <= (EnergyList[j] + deltaE))
           {
             //adding atom with number i to the configuration j
             SimulationBoxConfiguration[j].push_back(i);   
             //indicating that corresponding configuration have been created before
             ExistConfigurationIndex = 1; 
             
             //increasing number of the atoms in the configuration j by 1
             SimulationBoxConfiguration[j][0] = SimulationBoxConfiguration[j][0] + 1;
             
             //setting the energy of the atom configuration by indicating the number of element for the corresponding energy in the EnergyList array 
             AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2 - 1] = j;
             //setting number of the element in the SimulationBoxConfiguration vector with energy equal to EnergyList[j]      
             AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2] = SimulationBoxConfiguration[j][0] ;
           }  
         }
       
         //create new configuration
         if(ExistConfigurationIndex == 0)
         {
           //create new vector in the vector of vectors 
           SimulationBoxConfiguration.push_back(vector<int>());
           //increase number of the configuratuions in the system
           ExistingConfigurationNumber = ExistingConfigurationNumber + 1;
           //put energy in the EnergyList[ExistingConfigurationNumber] and number of the site in the SimulationBoxConfiguration[j][1]
           EnergyList.push_back(AtomSitesConfiguration[i]);
           SimulationBoxConfiguration[ExistingConfigurationNumber].push_back(1);
           SimulationBoxConfiguration[ExistingConfigurationNumber].push_back(i);
           
           //setting the energy of the atom configuration by indicating the number of element for the corresponding energy in the EnergyList array 
           AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2 - 1] = ExistingConfigurationNumber;
           //setting number of the element in the SimulationBoxConfiguration vector with energy equal to EnergyList[j]  
           AtomConfigurationNumber[(i - SubstrateBottomNumber) * 2] = 1;
         }   
         
       }  
      }
    }
 
  F1 = fopen("InitialConfig","w+");
  for(i = 1; i <= ExistingConfigurationNumber; i++)
  {
  	fprintf(F1,"%f  %d\n", EnergyList[i], SimulationBoxConfiguration[i][0]);
  	for(j = 1; j <= SimulationBoxConfiguration[i][0]; j++)
	  fprintf(F1,"%d  %d \n",j, SimulationBoxConfiguration[i][j]); 
	          
  	//getch();
  }
  fclose(F1);
    
 
}

void DeleteConfigurationElement( int AtomNumber, int ProcessIndex )
{
  int i, j, k, l;
  
  //printf("DELETE_CONF_ERROR1  AtomNumber = %d  ProcessIndex = %d\n", AtomNumber, ProcessIndex);   
  //getch();
  AtomConfigurationNumberOutput("ACheck1");
  i = AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2 - 1];
  j = AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2];
  k = SimulationBoxConfiguration[i][0];
  //printf("DELETE_CONF_ERROR1  AtomNumber = %d  ProcessIndex = %d\n", AtomNumber, ProcessIndex);   
  cout << "AtomNumber = " << AtomNumber << " i = " << i << " k = " << k << endl;
  cout << "SimulationBoxConfiguration[" <<i<<"][0] = " << SimulationBoxConfiguration[i][0] << endl;
  //getch();
  
  //updating configuration 'i' by rewriting element of the configuration j by last element of the vector 
  SimulationBoxConfiguration[i][j] = SimulationBoxConfiguration[i][k];
  //deleting last element of the vector (configuration i)
  SimulationBoxConfiguration[i].pop_back();
  //decrease number of the atoms in the 'i' configuration
  SimulationBoxConfiguration[i][0]= SimulationBoxConfiguration[i][0] - 1;
  if(SimulationBoxConfiguration[i][0] < 0)
  {
    printf("Error! Number of the atoms in the configuration with energy %f can't be less than zero\n",EnergyList[i]);
    //getch();
  }  
  
  //updating information in the AtomConfigurationNumber array
  AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2 - 1] = 0;
  AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2] = 0;
  l = SimulationBoxConfiguration[i][j];
  AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] = j;
  AtomConfigurationNumberOutput("ACheck2");
  SystemStateOutput("SystemData2");
  cout << "SystemDataOut" << endl;
} 


void MoveConfigurationElement( int AtomNumber, double Energy, int ProcessIndex, int AddMoveIndex)
{
  int i, j, k, l;
  int EnergyListNumber;         // number of the configuration in the EnergyList array
  int ExistConfigurationIndex;  // indicate whether there is corresponding configuration exists;  
                                // = 2 configuration exists
                                // = 0 configuration doesn't exist
  
  if(AddMoveIndex == 0) // = 1 - adding new element without deleting old
    DeleteConfigurationElement(AtomNumber, ProcessIndex); 
  
  //printf("ERROR_MOVE1  AtomNumber = %d  AddMoveIndex = %d\n", AtomNumber, AddMoveIndex);   
  //getch();
  // search for the correposnding energy number in the energyList array
  i = 1;
  ExistConfigurationIndex = 1;
  while(ExistConfigurationIndex == 1)
  {
    //cout << EnergyList[i] << " AAAA " << Energy << endl;
    //cout << i << "  " << ExistingConfigurationNumber << endl;
    //getch();
    if(Energy > (EnergyList[i] - deltaE) && Energy <= (EnergyList[i] + deltaE))
    {
    	ExistConfigurationIndex = 2;
    	EnergyListNumber = i;
    }     
	if(i == ExistingConfigurationNumber && ExistConfigurationIndex == 1)
      ExistConfigurationIndex = 0;
     
    i++;
  }  
  
  //cout << ExistingConfigurationNumber << "  "  << ExistConfigurationIndex << endl;
  //getch();
  //configuration exists
  if(ExistConfigurationIndex == 2)
  {
    //cout << "!!!!Errorr2" << endl; 
	AtomConfigurationNumberOutput("ConfigCheck");
    //SurfaceMapOutput("MapCheck");
    ConfigurationOutput("ConfigCheck");
    //RatesOutput("RatesCheck");
	//printf("MOVE_CONFGURATION_ERROR1\n");
    //getch();
	//add element to the corresponding configuration
    SimulationBoxConfiguration[EnergyListNumber].push_back(AtomNumber);
    //increase number of the elements in the 'i' configuration
    SimulationBoxConfiguration[EnergyListNumber][0]= SimulationBoxConfiguration[EnergyListNumber][0] + 1;
  
    //update AtomConfigurationNumber array
    AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2 - 1] = EnergyListNumber;
    AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2] = SimulationBoxConfiguration[EnergyListNumber][0];
  }
  
  if(ExistConfigurationIndex == 0)
  {
    //cout << "!!!!Errorr0" << endl; 
	AtomConfigurationNumberOutput("ConfigCheck");
    //SurfaceMapOutput("MapCheck");
    ConfigurationOutput("ConfigCheck");
    //RatesOutput("RatesCheck");
	//printf("MOVE_CONFGURATION_ERROR2\n");
    //getch();
	//create new vector in the vector of vectors 
    SimulationBoxConfiguration.push_back(vector<int>());
    //printf("MOVE_CONFGURATION_ERROR3\n");       
    //increase number of the configuratuions in the system
    //cout << " ExistingConfigurationNumber = " << ExistingConfigurationNumber << endl;
	ExistingConfigurationNumber = ExistingConfigurationNumber + 1;
    //put energy in the EnergyList[ExistingConfigurationNumber] and number of the site in the SimulationBoxConfiguration[j][1]
    //cout << " ExistingConfigurationNumber = " << ExistingConfigurationNumber << endl;
    //cout << "Energy = " << Energy << endl;
	EnergyList.push_back(Energy);
    //printf("MOVE_CONFGURATION_ERROR4\n");
    SimulationBoxConfiguration[ExistingConfigurationNumber].push_back(1);
    SimulationBoxConfiguration[ExistingConfigurationNumber].push_back(AtomNumber);
    //printf("MOVE_CONFGURATION_ERROR5\n");       
    //setting the energy of the atom configuration by indicating the number of element for the corresponding energy in the EnergyList array 
    AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2 - 1] = ExistingConfigurationNumber;
    //setting number of the element in the SimulationBoxConfiguration vector with energy equal to EnergyList[j]  
    AtomConfigurationNumber[(AtomNumber - SubstrateBottomNumber) * 2] = 1;
    //printf("MOVE_CONFGURATION_ERROR6\n");
  }
  //printf("ERROR_MOVE2  AtomNumber = %d  AddMoveIndex = %d\n", AtomNumber, AddMoveIndex);   
  //getch();
}

// update configuration around specific site (PositionNumber) without Volume diffusion
void SiteConfigurationUpdate( int PositionNumber, int ProcessIndex )
{
  int i, j, k, l, m;
  int AtomNumberType, NearestNeighborType, NextNearestNeighborType;
  int NearestSumAfter, NearestSumBefore, NextNearestSumBefore;
  int ExistConfigurationIndex;
  double AtomEnergy, NearestNeighborEnergy, NextNearestNeighborEnergy;
  int NextNeighborNumber;
  
  /*  - обновление конфигурации вокруг сайта происходит после изменения массива SimulationBox */ 
  /*  изменение массива AtomConfigurationNumber происходит в ф-иях MoveConfigNumber и DeleteConfNumber */
  /*  !SimulationBox->!AtomConfigurationNumber */
  //updating atom configurations around PositionNumber
  for(i = 3; i >= 0; i--)
  {
    j = 4 * PositionNumber - i;
    //number of the nearest and next nearest neighbors of the Atom in NewPosition
    k = NearestNeighborList[j];
    l = NextNearestNeighborList[j];
    
    //printf("SITE_CONFIGURATION_UPDATE_ERROR1  k = %d  l = %d  %d  %d\n", k, l, 
	//                                                        SimulationBox[k], SimulationBox[l]);   
    //getch();
    //atom types of the nearest and next nearest neighbors of the Atom in NewPosition
    NearestNeighborType = SimulationBox[k];
    NextNearestNeighborType = SimulationBox[l];
    //printf("SITE_CONFIGURATION_UPDATE_ERROR2  k = %d  l = %d\n", k, l);    
    
    //delete atom from the confiugration array, since he was there and now number of his neighbors is equal to 4
	if(NearestNeighborType != 0)
	{
	  NearestSumAfter = NearestNeighborCalc(k);
	  NearestSumBefore = NearestNeighborBeforeProcess(k);
	  //cout << "NearestSumAfter = " << NearestSumAfter << "NearestSumBefore = " << NearestSumBefore << endl;
	  //only delete atom from the Configurqation array
	  //since atoms with 4 neighbors are considered fixed
	  if( NearestSumAfter == 4 && NearestSumBefore < 4)
        DeleteConfigurationElement(k, ProcessIndex);
      //move atom in Configuration array without deleting  
      if(NearestSumAfter < 4 && NearestSumBefore == 4)
      {
         NearestNeighborEnergy = AtomSiteEnergy(k);
    	 MoveConfigurationElement( k, NearestNeighborEnergy, ProcessIndex, 1 ); 
      }  
      //move atom in Configuration array with deleting since he was in array
	  if(NearestSumAfter < 4 && NearestSumBefore < 4)
      {
         //cout << "ERRORRRR!" << endl;
         //getch();
		 NearestNeighborEnergy = AtomSiteEnergy(k);
		 //cout << "AtomConfigurationNumber[(k - SubstrateBottomNumber) * 2] = " << 
		 //       AtomConfigurationNumber[(k - SubstrateBottomNumber) * 2]<< endl;
        // getch();
    	 if(AtomConfigurationNumber[(k - SubstrateBottomNumber) * 2] == 0)
    	 {
		   //cout << "ERRORRRR3!" << endl;
           //getch();
           MoveConfigurationElement( k, NearestNeighborEnergy, ProcessIndex, 1 );
		 }
		 if(AtomConfigurationNumber[(k - SubstrateBottomNumber) * 2] != 0)  
		 {
		   //AtomConfigurationNumberOutput("AtomConfNumber1");
		   //cout << "ERRORRRR4!" <<  (k - SubstrateBottomNumber) * 2 << "  "
		   //    << AtomConfigurationNumber[(k - SubstrateBottomNumber) * 2] << endl;
           //getch();
		   MoveConfigurationElement( k, NearestNeighborEnergy, ProcessIndex, 0 );
         }
      }
    }
	/*  old variant
	if(NearestNeighborType != 0)
    {
      //atom energies of the nearest and next nearest neighbors of the Atom in NewPosition
      //NearestSumm = NearestNeighborCalc(  )
      NearestNeighborEnergy = AtomSiteEnergy(k); 
      //printf("SITE_CONFIGURATION_UPDATE_ERROR3  k = %d  l = %d\n", k, l);    
      //getch();
	  if(NearestNeighborEnergy == -1.0)     //delete configuration if atom has maximum number of NearestNeighbors
	  {
	  	//printf("SITE_CONFIGURATION_UPDATE_ERROR4  k = %d  l = %d\n", k, l);    
		printf("SITE_CONFIGURATION_UPDATE_ERROR2  k = %d  l = %d\n", k, l);
		getch();
		DeleteConfigurationElement(k, ProcessIndex);
	  	//printf("SITE_CONFIGURATION_UPDATE_ERROR5  k = %d  l = %d\n", k, l);    
		  //printf("SITE_ERROR2\n");   
        //getch();
	  }
	  if(NearestNeighborEnergy != -1.0) 
	  {
	  	printf("SITE_CONFIGURATION_UPDATE_ERROR3  k = %d  l = %d\n", k, l); 
		NearestSum = NearestNeighborBeforeProcess( k, PositionNumber);
		//atom had had less than 4 neighbors before process occured, so it's in the Configuration array
		// it is necessary to remove atom from the configuration array, AddMoveIndex = 0
		if(NearestSum != 4)
		  MoveConfigurationElement( k, NearestNeighborEnergy, ProcessIndex, 0 );
		//atom had had 4 neighbors before process occured , so it's not in the configuration array and 
		//it doesn't need to remove him form the Configuration array, AddMoveIndex != 0
		if(NearestSum == 4)
		  MoveConfigurationElement( k, NearestNeighborEnergy, ProcessIndex, 1 );
	  	 
		  //printf("SITE_ERROR2\n");   
        //getch();
	  }  
      
    }
    */
    //Пришел дальний сосед!, конфигурация изменилась
    //calculation of the nearest neighbors of the NextNeighbor atom
    //cout << "NNext" << endl;
    if(NextNearestNeighborType != 0)
    {
       NearestSumAfter = NearestNeighborCalc(l);
	   NearestSumBefore = NearestNeighborBeforeProcess(l);
	   //NextNearestSumBefore = NextNearestNeighborBeforeProcess(l);
	   //cout << "NearestSumAfter = " << NearestSumAfter << " NearestSumBefore = " << NearestSumBefore << endl;
	   if( NearestSumAfter < 4 && NearestSumBefore == 4)
	   {
	     NextNearestNeighborEnergy = AtomSiteEnergy(l); 	
	     MoveConfigurationElement( l, NextNearestNeighborEnergy, ProcessIndex, 1 );
	   }
	   if( NearestSumAfter < 4 && NearestSumBefore < 4)
	   {
	     //cout << "AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2]" <<
	     //        AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] << endl;
		 NextNearestNeighborEnergy = AtomSiteEnergy(l); 	
	     if(AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] == 0)
	       MoveConfigurationElement( l, NextNearestNeighborEnergy, ProcessIndex, 1 );
	       
		 if(AtomConfigurationNumber[(l - SubstrateBottomNumber) * 2] != 0)
		   MoveConfigurationElement( l, NextNearestNeighborEnergy, ProcessIndex, 0 );
	   }
	   
	   if(NearestSumAfter == 4 && NearestSumBefore == 4)
	   {
	   	 //NextNearestNeighborEnergy = AtomSiteEnergy(l); 
	   	 ;
	   }  
	   //printf("ENERGY_SITE_ERROR2  %d  %f\n", NextNeighborNumber, NextNearestNeighborEnergy);   
       //getch(); 
	}
	
  }  
}

//updating configuration of the whole box
void BoxConfigurationUpdate( int OldPositionNumber, int NewPositionNumber, int ProcessIndex )
{ 
  // OldPositionNumber - number of the site before the process occurs
  // NewPositionNumber - number of the site after the  process occurs
  // OldPosition = -1 in case of the deposition
  // NewPosition = -1 in case of the evaporation
  
  int i, j, k, l, m;
  int AtomNumberType, NearestNeighborType, NextNearestNeighborType;
  int ExistConfigurationIndex;
  double AtomEnergy, NearestNeighborEnergy, NextNearestNeighborEnergy;
  
    
  if(OldPositionNumber == -1 && NewPositionNumber == -1)
  {
    printf("Deposition/Evaporation process Error\n");
    //getch();
  }
  
  //in the case of the deposition
  if(OldPositionNumber == -1)
  {
  	AtomEnergy = AtomSiteEnergy(NewPositionNumber); 
  	//cout << AtomEnergy << endl;
  	//getch();
	MoveConfigurationElement( NewPositionNumber, AtomEnergy, ProcessIndex, 1);  // AddMoveIndex = 1, 
	                                                                            // therefore we add deposited atom configuration
	SiteConfigurationUpdate( NewPositionNumber, ProcessIndex ); //update configuration around the site
  }
   
  //in the case of the evaporation
  if(NewPositionNumber == -1)
  {
  	DeleteConfigurationElement( OldPositionNumber, 2 ); // delete configuration of the evaporated atom
	SiteConfigurationUpdate( OldPositionNumber, ProcessIndex );
  }
   
  //in the case of the diffusion  
  if(OldPositionNumber != -1 && NewPositionNumber != -1)  
  {
    cout << "MOVE_CONF_ERROR1" << endl;
	//getch();
	
	//update old position information
	DeleteConfigurationElement( OldPositionNumber, 2 ); // delete configuration for the old position of diffused atom
	cout << "MOVE_CONF_ERROR2" << endl;
	//getch();
	
	SiteConfigurationUpdate( OldPositionNumber, ProcessIndex ); // update configuration around old position
	cout << "MOVE_CONF_ERROR3" << endl;
	//getch();
	//ne dobavljaem element v configuration array poskol'ku pri obnowlenii staroj
	//pozicii on popadaet tuda awtomaticheski
	SiteConfigurationUpdate( NewPositionNumber, ProcessIndex ); 
	
	//calculate new position energy
	/*
	AtomEnergy = AtomSiteEnergy(NewPositionNumber); //calculate energy of the new position, 
	if(AtomEnergy != -1)
	{
	  MoveConfigurationElement( NewPositionNumber, AtomEnergy, ProcessIndex, 1);  // AddMoveIndex = 1, 
	                                                                              // therefore we add new configuration
      SiteConfigurationUpdate( NewPositionNumber, ProcessIndex );                 //update information around new position
	}
	if(AtomEnergy == -1)
	{
	  SiteConfigurationUpdate( NewPositionNumber, ProcessIndex ); 
	}
	*/
  }
   
}



