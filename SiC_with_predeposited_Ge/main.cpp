#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <graphics.h>
//#include <winbgim.h>
#include <iostream>
#include <fstream>

//includes of the KMC programm
#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"
#include "additional.h"
#include "init.h"
#include "impurity.h"
#include "cluster.h"
//#include "draw.h"
//#include "additional.h"



int main( void )
{
  //cheking!!!
  FILE *F1, *F2, *F3;
  
  // System announcments
  FILE *System;
  

  int i, j, k;
  //FILE *F;
  //covergae of the substrate surface
  double coverage = 0.5;
  
  
  double totd;
  double prob;
  
  //double type step counter for output
  double c_k = 0.0; 
  double r = 1.0;
  double r_count;
  
  //variable for the cheking output condition
  double check;
  //double output step
  double outstep;
  
  // counters of the time steps 
  //'count' - current time step
  // 'count_total' - total number of time steps
  int count_total;
  
  //current adatom number in main() function
  int ad_curr;
  //number of atom with same number of neighbors before and after the hop
  int ad_same;
  //total number of deopsited atoms
  int atoms_deposited;
  
  //
  int cluster_data_flaq;
  
  //testing
  int proc;
    
  FILE *Neigh;

  printf("begin\n");
  //opening of output files
  //F = fopen("boxandvancat","w+");
  //F1 = fopen("Deposition","w+");
  //F2 = fopen("Diffusion","w+");
  //Neigh = fopen("NeighborL","w+");
  //F3 = fopen("process","w+");  
  //RateData = fopen("RateData","w+");  
  //F5 = fopen("GroupList","w+");
  //System = fopen("SystemInf","w+");
  //Frates = fopen("Rates","w+");
  
  //inicialisation of the variables
  //calculation of total amount of deposited atoms from coverage
   xcell = 100;
   ycell = 100;
   zcell = 1;
   /*
   totd = 10; //(coverage * xcell * ycell * 2);  
   total = (int)totd;
   */
   
   /* 
   setting values of the counters and arrays
   */
   //setting values of the counters
   count = 0;
   //outstep = 100000.0;
   count_total = 5000;
   grouplist_count = 0;
   sum_check = 0;
   cc = 0;
   global_counter = 0;
   atom_counter = 1;
   atoms_deposited = 1;
   ngh = 0;
   
   //number of steps (atoms deposited) before output
   outstep = 50;
   
   //setting initial values of the nearest neighbors and next nearest neighbors arrays
   for(i = 0; i <= 4; i++)
   {
      neighbor[i] = Nextneigh [i] = 0;
   }

   //zadanie coordinat of simulation box i vacant layer 
   MakeTopLayerCoords(xcell,ycell,zcell);  
   VacantPosCoords(xcell, ycell, zcell);
   InitDone();
     
   //initialisation of random number generator
   srand(time(NULL));  
   
   cluster_data_flaq = 1;
   if(cluster_data_flaq == 1)
   {
     DataRead();
     ClusterData();
     ClusterDistribution();
     return 0;  
   }  
   
   /* ustanovlenie nachal'nogo sostojanija */
   /*
   ImpPreDeposition(100);
   while(ngh == 0)
     ngh = Deposition(1);
   */
  /* Calculation of an evolution of the system */
  //count = 1100;
  count = -1;
  while(count <= count_total)
  {
    //determination of the system configuration on step 'count'
    GroupList();
    //Calculation of the process, which will occur
    process = Rates(); 
    //if(process > 2)
      //fprintf(Frates,"%d\n",process);
    //fprintf(Frates,"r=%.10f   r0 = %.10f  r1=%.10f   rImp1=%.10f\n", r, r0, r1, rImp1);

    printf("deposited=%d  process=%d\n",count + 1, process);
    //printf("\n");
    //fprintf(Frates,"%d\n", process);
    
    /* deposition of the atom carbon on the substrate */
    if(process == 1)
    {
       ngh = 0;
       
       /*
       //deposition with reflection from occupied sites
       ngh = Deposition(1);
       if(ngh != 0) 
       {
         atoms_deposited++; 
         //increase of time step counetrs
         count++;    
         c_k = c_k + 1.0;
         
         if(c_k == 100.0)
         {
           printf("%d \n",adatom_k);
           //getch();
         }  
       }
       */
       // deposition without reflection from occupied sites  
       
       while(ngh == 0)
         ngh = Deposition(1);   
       
       count++;    
       c_k = c_k + 1.0;          
       //fprintf(System,"step=%d   process=%d,  atoms_deposited=%d\n",count,process,atoms_deposited);  
    }
    
    /* diffusion of the adatom or an atom on the substrate */ 
    if(process != 1)
      ngh = Diffusion(process,1);
    
    global_counter++;
    
    // proverka uslovie vyvoda simulation box'a
    check = r * outstep;    
    //output of simulation box after 'check' time steps
    if(c_k == check)
    {
       DegreeR(check,8.0);
       OutStep();
       r++;
    }

  }
 
  //output of final state of the system
   out();
   
  //fclose(F);  
  //fclose(F1);   
  //fclose(F2);  
  //fclose(output);
  //fclose(Diff); 
  //fclose(F3);   
  //fclose(RateData);
  //fclose(F5);
  //fclose(System);
  //fclose(Frates);
   
  printf("End of the calculation\n");
  return 0;
}
