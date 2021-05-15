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
#include "configurations.h"
#include "init.h"
#include "neighbors.h"
#include "outdata.h"
#include "parameters.h"
#include "processes.h"
#include "simulation_box.h"
#include "unit_cell.h"

#define kT0 0.025 // im meV

using namespace std;
									
									/*  NB!  */
/*  - обновление конфигурации вокруг сайта происходит после изменения массива simulationBox */ 
//function of the deposition only to the sites having nearest neighbor atoms
void DepositionNearest2D( int DepositionSiteNumber, int SurfaceMapNumber, int AtomType)
{
  int i, j, k;
  
  SimulationBox[DepositionSiteNumber] = AtomType;
  //printf("%d   %d\n", DepositionSiteNumber, AtomType);
  //getch();
  SurfaceMapEvolution2D( DepositionSiteNumber,  SurfaceMapNumber, 1 );
  //printf("DEPOSITION_ERROR1\n");
  //getch();
  BoxConfigurationUpdate( -1, DepositionSiteNumber, 1 );
  //printf("DEPOSITION_ERROR_LAST\n");
  //getch();
}

void Evaporation( int EvaporationSiteNumber )
{
  
  SurfaceMapEvolution2D( EvaporationSiteNumber,  -1, 2 ); // update of the surface map in the case of the evaporation
	                                                      // doesn't need SurfaceMapNumber  
  SimulationBox[EvaporationSiteNumber] = 0; 
  /* NB!  SurfaceMapEvolution2D must be done before setting SimulationBox[EvaporationSiteNumber] = 0
  since all calculations in this function are done when SimulationBox[EvaporationSiteNumber] != 0  
  */	                                                      
  BoxConfigurationUpdate( EvaporationSiteNumber, -1, 2 );
  
  //printf("Error\n");   
  //getch();
}

// function for calculation new diffusion position for the atom in the SiteNumber position
int NewDiffusionDirection( int SiteNumber )
{
  int i, j, k, l, m, n;
  int FreeSite[10],         //array of free sites aroundSiteNumber
      NewDiffusionSite[10]; // array of the free sites with at least 1 nearest neighbors
  int NearestSumm;          // number of the nearest neighbor for the atoms from the FreeSite  array    
  int RandNumber;
    
  m = 1;
  if(NearestNeighbornumber == NextNearestNeighbornumber)
  {
	for(i = (NearestNeighbornumber - 1); i >=0; i--)
	{ 
	  j = SiteNumber * 4 - i;
      k = NearestNeighborList[j];
      l = NextNearestNeighborList[j];
      
      //put empty nearest and next nesrest neighbor positions in the NewDiffusion array
	  if(SimulationBox[k] == 0)
      {
        FreeSite[m] = k;
        m++;		
      }
      if(SimulationBox[l] == 0)
      {
        FreeSite[m] = l;
        m++;		
      }
	}	
  }
  
  //decrease m by 1 to get total amount of empty sites around SiteNumber, since m = 1 at the beginning
  m = m - 1;
  
  n = 1;
  //atoms can hop only in the site with at least 1 nearest neighbor
  for(i = 1; i <= m; i++)
  {
    NearestSumm = NearestNeighborCalc( FreeSite[i] );
	if(NearestSumm >= 1)
  	{
  	  NewDiffusionSite[n] = FreeSite[i];
	  n = n + 1;		
  	}  	  
  }
  //decrease n by 1 to get total amount of empty sites with at least 1 nearest neighbor around SiteNumber
  n = n - 1;
  
  //choose new diffusion direction
  RandNumber = 1 + rand()%n;
  
  return NewDiffusionSite[RandNumber];
  
}

void SurfaceDiffusion2D( int DiffusionSiteNumber )
{
	int i, j, k, l, m, n;
	int NewDiffusionSiteNumber;
	int AtomType;
	 
	AtomType = SimulationBox[DiffusionSiteNumber];
	cout << "AtomType=" << AtomType << endl;
	cout << "DiffusionSiteNumber=" << DiffusionSiteNumber << endl;
	SiteStateBefore = 0;
	SiteStateAfter = 1;
	
	//cout << "SurfDiff_1" << endl;
	//getch();
	//update surface map and simulation box arrays for the new position of the diffusing atom	
	NewDiffusionSiteNumber = NewDiffusionDirection( DiffusionSiteNumber );
	cout << "NewDiffusionSiteNumber=" << NewDiffusionSiteNumber << endl;
    SimulationBox[NewDiffusionSiteNumber] = AtomType;
  
	//cout << "SurfDiff_2" << endl;
	//getch();
	//cout << " NewDiffusionSiteNumber = " << NewDiffusionSiteNumber << endl;
	//getch();
	//SurfaceMapOutput("SurfaceMap1");
	
	i = 1;
    while(NewDiffusionSiteNumber != SurfaceMap[i])
    {
   	   	  
	  i++;
	  
    }	 
    cout << "i = " << i << " SurfaceMapLength = " << SurfaceMapLength << endl;
    //cout << "After cycle" << endl;
    //cout << "NewDiffusionSiteNumber = " << NewDiffusionSiteNumber << "  SurfaceMap[i] = " << SurfaceMap[i] << endl;
    //getch(); 
	//cout << "SurfDiff_4" << endl;
    //cout << "SurfaceMap[i] = " << SurfaceMap[i] << endl;
	
	
	//diffusion is the sum of two processes evaporation + diffusion
	//update surface map and simulation box arrays for the old position of the diffusing atom	
	SurfaceMapEvolution2D( DiffusionSiteNumber,  -1, 2 );  // update of the surface map in the case of the evaporation
	                                                       // doesn't need SurfaceMapNumber      
	SimulationBox[DiffusionSiteNumber] = 0;
	
	//cout << "Diff_Evap" << endl;
	//getch();
	
	//Приоритеты между SurfaceMapEvolution и SiteConfigurationUpdate!!!!!!!!!
	// в SME считается, что еще на месте(при испарении)!!!!!!!!!!
	SurfaceMapEvolution2D( NewDiffusionSiteNumber, i, 1 );
	//cout << "SME2D"<< endl;
	//cout << " DiffusionSiteNumber = " << DiffusionSiteNumber << "  NewDiffusionSiteNumber = " << NewDiffusionSiteNumber << endl;
	//getch();
	
	BoxConfigurationUpdate( DiffusionSiteNumber, NewDiffusionSiteNumber, 3);
	//cout << "Diff_Conf" << endl;
	//getch();
	//cout << "SurfDiff_4" << endl;
	//getch();
	
}

void InitialRatesCalculation( int ConfigurationFrequency )
{
  int i, j, k;
  double Rate;             //value of the process rate
  double SummaryRate;      // summ of all processes in the system rates
  double FrequencyValue;   // frequency of the atom of the iwth the  type AtomTypeFreq
  int AtomTypeFreq;        // type of the atom for the frequency calculation

  //zero elements don't use
  ProcessRate.push_back(0.0);    
  //Flux[0] = 0.0;
  SummaryRate = 0.0;
  DepositionTotalRate = 0.0;
  
  //printf("Initial rates calculation ERROR1\n");
  //getch();
     
  if(Dimension == 2)
  {
    //REWRITE universally for diffrent crystal structures
    if(CrystalStructureType == 5)
    {
      for(i = 1; i <= AtomTypesNumber; i++) 
      {
        Flux[i] = Flux[i] * (1E-8 * LatticeParameter) * (double)BoxSizeX / 2.0;  // deposition flux in the Sites per second
        //cout << Flux[i] << endl;
		//cout << UnitCellSizeX * UnitCellSizeY * LatticeParameter * LatticeParameter / 2.0 << endl;
		SummaryRate = SummaryRate + Flux[i];
		ProcessRate.push_back(SummaryRate);  
		//cout <<  "Flux[i] = " << Flux[i] << endl;      
		//getch();
      }  
    }   
  }
  
  //if(Dimension == 3)
    //;
  //printf("Initial rates calculation ERROR2\n");
  //getch();  
  
  DepositionTotalRate = SummaryRate;
  //No evaporation, Debye frequencies doesn't depend on atom type for specific configuration 'i' 
  if( PhysicalProcessIndex[4] == 0 && ConfigurationFrequency == 1 )
  {
    
    TotalProcessNumberPrevious = AtomTypesNumber + ExistingConfigurationNumber;
    //printf("Initial rates calculation ERROR3\n");
    //printf("ExistingConfigurationNumber = %d\n", ExistingConfigurationNumber);
    //getch(); 
    //setting values for the rates of different diffusion processes
    for(i = 1; i <= ExistingConfigurationNumber; i++)
    {
      AtomTypeFreq = SimulationBox[SimulationBoxConfiguration[i][1]];
	  FrequencyValue = AtomFrequency[AtomTypeFreq];
	  Rate = SimulationBoxConfiguration[i][0] * FrequencyValue * exp( (-1.0 * EnergyList[i]) / ((Temperature / 300) * kT0) );
      SummaryRate = SummaryRate + Rate;
	  /*
	  cout << " SimulationBoxConfiguration[i][0] = " << SimulationBoxConfiguration[i][0] << endl;
	  cout << " FrequencyValue = " << FrequencyValue << endl;
	  cout << " EnergyList[i] = " << EnergyList[i] << endl;
	  cout << " Exp = " << exp( (-1.0 * EnergyList[i]) / ((Temperature / 300) * kT0)) << endl; 
	  cout << " Rate =  " << Rate << endl;
	  cout << " SummaryRate = " << SummaryRate << endl;
	  getch();
	  */
	  //DiffusionRate.push_back(Rate);
	  ProcessRate.push_back(SummaryRate);
    }
    
    for(j = 1; j <= TotalProcessNumberPrevious; j++)
    {
       ProcessRate[j] = ProcessRate[j] / SummaryRate;
	   //cout << "ProcessRate = " << ProcessRate[j] << "   SummaryRate = " << SummaryRate << endl;	  
    }
    
  }
   //printf("Initial rates calculation ERROR_LAST\n");
   //getch(); 
   //cout <<  "TotalProcessNumberPrevious = "  <<  TotalProcessNumberPrevious << endl;
   //getch();
   
}

void RatesUpdate( int ConfigurationFrequency )
{
	int i, j, k;
	double Rate;
	double SummaryRate;
	double FrequencyValue;   // frequency of the atom of the iwth the  type AtomTypeFreq
    int AtomTypeFreq;        // type of the atom for the frequency calculation
	
	//cout << "Rates update" << endl;
	//getch();
	SummaryRate = DepositionTotalRate;	
	
	//cout << " DepositionTotalRate = " << DepositionTotalRate << endl;
	//cout << " TotalProcessNumberPrevious = " << TotalProcessNumberPrevious << endl;
	//cout << " AtomTypesNumber = " << AtomTypesNumber << endl;
	//cout << " TotalProcessNumber = " << TotalProcessNumber << endl;
	//getch();
	
	TotalProcessNumber = AtomTypesNumber + ExistingConfigurationNumber;
	for(i = 1; i <= (TotalProcessNumberPrevious - AtomTypesNumber); i++)
	{
		AtomTypeFreq = SimulationBox[SimulationBoxConfiguration[i][1]];
	    FrequencyValue = AtomFrequency[AtomTypeFreq];
		//cout << " AtomTypeFreq =" << AtomTypeFreq << endl;
		//getch();
		Rate = SimulationBoxConfiguration[i][0] * FrequencyValue * exp( (-1.0 * EnergyList[i]) / ((Temperature / 300) * kT0) );
		SummaryRate = SummaryRate + Rate;
		ProcessRate[i + AtomTypesNumber] = SummaryRate;
		/*
		cout << " SimulationBoxConfiguration[i][0] = " << SimulationBoxConfiguration[i][0] << endl;
	    cout << " FrequencyValue = " << FrequencyValue << endl;
	    cout << " EnergyList[i] = " << EnergyList[i] << endl;
	    cout << " Exp = " << exp( (-1.0 * EnergyList[i]) / ((Temperature / 300) * kT0)) << endl; 
	    cout << " Rate =  " << Rate << endl;
	    cout << " SummaryRate = " << SummaryRate << endl;
	    getch();                
		*/
	}
	
	//cout << "EVOLUTION_RATES2" << endl;
	//getch();
	ConfigurationOutput("ConfigRates");
	if(TotalProcessNumber > TotalProcessNumberPrevious)
	{
	  i = TotalProcessNumberPrevious;
	  while( i <= (TotalProcessNumber - TotalProcessNumberPrevious))
	  {
	    AtomTypeFreq = SimulationBox[SimulationBoxConfiguration[i][1]];
	    FrequencyValue = AtomFrequency[AtomTypeFreq];
		Rate = SimulationBoxConfiguration[i][0] * FrequencyValue * exp( (-1.0 * EnergyList[i]) / ((Temperature / 300) * kT0) );
	    SummaryRate = SummaryRate + Rate;
	    ProcessRate.push_back(SummaryRate);		
	    i++;
	    //cout << "i = " << i;
		//cout << " TotalProcessNumberPrevious = " << TotalProcessNumberPrevious << endl;
		//cout << " TotalProcessNumber = " << TotalProcessNumber << endl;
		//getch();
	  }
	  TotalProcessNumberPrevious = TotalProcessNumber;
	}
	
	//cout << "EVOLUTION_RATES3" << endl;
	//getch();
	for(i = 1; i <= AtomTypesNumber; i++)
      ProcessRate[i] = Flux[i] / SummaryRate;	
	
	//cout << "EVOLUTION_RATES4" << endl;
	//getch();
	for(i = AtomTypesNumber + 1; i <= TotalProcessNumber; i++)
	  ProcessRate[i] = ProcessRate[i] / SummaryRate;
    //cout << "EVOLUTION_RATES5" << endl;
	//getch();
   
}

void SystemEvolution( void )
{
   int i, j, k, l;
   int ProcessFlag;
   int ProcessRealisationNumber; //number of the process to realise, corresponding to the process rate in ProcessRate[] array
   //int RandNumber1, RandNumber2;
   double RandNumber1;	
   int AtomNumber;
   int Type;
	
   //cout << "SYSTEM EVOLUTION " << endl;
   AtomConfigurationNumberOutput("EVOL_Check1");
   RatesUpdate( 1 );
   //cout << "System evolution 2" << endl;
   //getch();
   
   //no evaporation, Debye frequencies doesn't depend on atom type for specific configuration 'i' 
   RandNumber1 = ((double)rand()) / ((double)RAND_MAX);
   //cout << "SYSTEM EVOLUTION2" << endl;
   //getch();
   i = 0;
   ProcessFlag = 0;
   while(ProcessFlag != 1)
   {
   	 //cout  << "i = " << i << "  ProcessRate[i] = " << ProcessRate[i] <<"  RandNumber1 = " << RandNumber1 << 
	 //	  "  ProcessRate[i + 1] = " << ProcessRate[i + 1] << endl;
	 //getch();
	 if(RandNumber1 > ProcessRate[i] && RandNumber1 <= ProcessRate[i + 1])
	 {
	   ProcessFlag = 1;	
	   ProcessRealisationNumber = i + 1;
	   //cout << "ProcessRealisationNumber = " << ProcessRealisationNumber << endl;
       //getch();
   	 }
   	 
	 i = i + 1; 
   }
   //|cout << ProcessRealisationNumber << endl;
   //getch();
                  /* Evaporation and volume diffusion is excluded */
   if(ProcessRealisationNumber >= 1 && ProcessRealisationNumber <= AtomTypesNumber)
   {
   	 AtomNumber = 1 + rand()% SurfaceMapLength;    //random atom position on the surface
   	 Type = ProcessRealisationNumber; // random atom type 
     //cout << "Type = " << Type << endl;
     //cout << "AtomNumber = " << AtomNumber << endl;
     //getch();
	 DepositionNearest2D( SurfaceMap[AtomNumber], AtomNumber, Type);
     
   }  
   //cout << "SYSTEM EVOLUTION3" << endl;
   //getch();
   if(ProcessRealisationNumber >= (AtomTypesNumber + 1) && ProcessRealisationNumber <= TotalProcessNumber)
   {
   	 j = ProcessRealisationNumber - AtomTypesNumber;
	 AtomNumber = 1 + rand()%SimulationBoxConfiguration[j][0];
   	 //cout << "j = " << j << endl;
	 //cout << "AtomNumber = " << AtomNumber << endl;
	 //cout <<  "SimulationBoxConfiguration[" <<j <<"][0]=" << SimulationBoxConfiguration[j][0] << endl;
   	 //cout << "DIFF_ERROR0.5" << endl;
   	 /*
		FCheck = fopen("CHECK","w+");
   	 fprintf(FCheck,"0  %d\n",SimulationBoxAtomsNumber);	 
		for(l = 1; l <= SimulationBoxAtomsNumber; l++)
   	 {
   	 	
		fprintf(FCheck,"%d  %d\n",l, SimulationBox[l]);	 
   	 }
	  
	 fclose(FCheck);  
	 */
		//cout << SurfaceMapLength << endl;
		//getch();
		//cout << "j = " << j << "  AtomNumber = " << AtomNumber <<  "  " << SimulationBoxConfiguration[j][AtomNumber]<< endl;
   	 //cout << "SurfaceMapLength = " << SurfaceMapLength << endl;
	 //cout << "DIFF_ERROR2" << endl;
	 //SurfaceMapOutput("SurfaceMapEvolution1");
	 //cout << "DIFF_ERROR1.5" << endl;
	 /*
	 FCheck2 = fopen("CHECK2","w+");
   	 fprintf(FCheck2,"0  %d\n",SimulationBoxAtomsNumber);	 
		for(l = 1; l <= SimulationBoxAtomsNumber; l++)
   	 {
   	 	
		fprintf(FCheck2,"%d  %d\n",l, SimulationBox[l]);	 
   	 }
	  
	 fclose(FCheck2);  
	 */
	 //getch();
   	 //AtomConfigurationNumberOutput("Check1");
   	 cout << "DIFF_ERROR2.5" << endl;
	 //getch();
	 cout << "SimulationBoxConfiguration["<<j<<"]["<<AtomNumber<<"] =" << SimulationBoxConfiguration[j][AtomNumber]<< endl;
	 SurfaceDiffusion2D( SimulationBoxConfiguration[j][AtomNumber] );
   	 cout << "DIFF_ERROR3" << endl;
   }
   //cout << "SYSTEM EVOLUTION4" << endl;
   //getch();
 
}

