#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <graphics.h>
//#include <winbgim.h>
//#include <iostream>

//includes of the KMC programm
#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"
#include "additional.h"
#include "init.h"
//#include "draw.h"
//#include "additional.h"



int main( void )
{
  //cheking!!!
  FILE *F1, *F2, *F3;
  
  // System announcments
  FILE *System;

  //printf("Hello");
  int i;
  //FILE *F;
  //covergae of the substrate surface
  double coverage = 0.5;
  
  double totd;
  double prob;
  
  char a;
  
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
  // 'count_draw' - redraw of picture each 'count_draw' steps
  int count, count_total, count_draw;
  
  
  //current adatom number in main()
  int ad_curr;
  //number of atom with same number of neighbors before and after the hop
  int ad_same;
  
  FILE *Neigh;
  
  //vyvod ASCII symvolov ot 1 do 100
  /*FILE *symbols;
  char c;
  
  symbols = fopen("symb","w+");
  for(i = 1; i < 100; i++)
  {
    c = i;
    fprintf(symbols,"%d  %c\n", i, c);
  }
  fclose(symbols);
  */
  
  printf("begin\n");
  //opening of output files
  //F = fopen("boxandvancat","w+");
  F1 = fopen("Deposition","w+");
  F2 = fopen("Diffusion","w+");
  Neigh = fopen("NeighborL","w+");
  F3 = fopen("process","w+");  
  RateData = fopen("RateData","w+");  
  F5 = fopen("GroupList","w+");
  System = fopen("SysAnounce","w+");
  ch = fopen("chlist","w+");
  
  //F2 = fopen();
  
  //inicialisation of the variables
  //calculation of total amount of deposited atoms from coverage
   xcell = 50;
   ycell = 50;
   zcell = 1;
   /*
   totd = 10; //(coverage * xcell * ycell * 2);  
   total = (int)totd;
   */
   
   grouplist_count = 0;
   sum_check = 0;
   cc = 0;
   global_counter = 0;
   
   for(i = 0; i <= 4; i++)
   {
      neighbor[i] = Nextneigh [i] = 0;
   }
   
   for(i = 1; i <= xcell * ycell * 2; i++)
   {   
      NeighborList(i);
      //fprintf(Neigh,"for %d nearest %d  %d  %d  %d\n",i,  neighbor[1], neighbor[2], neighbor[3],neighbor[4]);
      //fprintf(Neigh,"for %d next nearest neugh %d  %d  %d  %d\n",i,  Nextneigh[1], Nextneigh[2], Nextneigh[3],Nextneigh[4]);
      //printf("%d  %d  %d  %d  %d\n",i, neighbor[1], neighbor[2], neighbor[3],neighbor[4]);
   }
  
    
  
   for(i = 0; i < 8; i++)
     timedata[i] = '0';
   
   //diffsuion count number
   count = 0;
   outstep = 100.0;
   count_total = 10000;
   count_draw = 5;
   
   //zadanie coordinat of simulation box i vacant layer 
   MakeTopLayerCoords(xcell,ycell,zcell);  
   VacantPosCoords(xcell, ycell, zcell);
   InitDone();
   
   
   atom_counter = 1;
   
   //initialisation of random number generator
   srand(time(NULL));  
   //ustanovlenie schetchikov razlichnyh sortov v nol'
   //adatom_k = atom1_k = atom2_k = atom3_k = 0;
   //adatom_kN = atom1_kN = atom2_kN = atom3_kN = 0;
   
   
   //ustanovlenie nachal'nogo sostojanija
   GroupList();
   ngh = Deposition(xcell,ycell,zcell);
   //printf("ngh = %d\n",ngh);
   NeighborList(ngh);
   //printf("%d  %d  %d  %d\n", neighbor[1],neighbor[2],neighbor[3],neighbor[4]);
   
   out();
  
   //testing the rates  
   //process = Rates();
   //printf("process = %d  \n",process);
  
  
  
  //testing og hte programm
  //printf("continue\n");
  //printf("%d   %d\n",count,count_total);
  
  // drawing of the substrate
  //initwindow(810,810); 
  //DrawMesh();
  //while(!kbhit());  
  //testing of the programm before simulation of the growth
  //count = count_total + 1;
    
  //evolution of the system
  while(count <= count_total)
  {
    
    GroupList();
    process = Rates();
    printf("process = %d  \n",process);
    //updating the system
    
    //fprintf(F3,"%d   %d\n",count,process);
    
    
    
    if(process == 1)
       ngh = Deposition(xcell,ycell,zcell);
    if(process != 1)
      ngh = Diffusion(process);
    
    //increase of time step counetrs
    count++;
    //printf("count = %d  \n", count);
    c_k = c_k + 1.0;
    global_counter++;
    
    // proverka uslovie vyvoda simulation box'a
    check = r * outstep;
    
    //output of simulation box after 'check' time steps
    if(c_k == check)
    {
       DegreeR(check,8.0);
       OutStep( timedata, adcoords );
       r++;
       fprintf(System,"Outstep call   %f  %f  %f\n",c_k,r*outstep, check);
    }
    
   
    //visual testing of the programm
    count_draw = 1;
   
    /*
    if(count_draw == 1)
    {  
      DrawMesh(); 
      DrawEvolution(ngh);
      //delay(1);
    }
   
   */
    
    
  }
  
  
 
  //while(!kbhit());
  //closegraph();
  //out();
  
  //fclose(F);  
  fclose(F1);   
  fclose(F2);  
  fclose(output);
  //fclose(Diff); 
  fclose(F3);   
  fclose(RateData);
  fclose(F5);
  fclose(System);
  fclose(Neigh);
  fclose(ch);
  
  
  
  return 0;
}
