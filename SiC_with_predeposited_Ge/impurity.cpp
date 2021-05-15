#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>

#include "impurity.h"
#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"

//impurity atoms array
int Impur[90000];

//number of deposited impurities
int ImpNumb;

//number of the atoms with 1 impurity neighbor
int atom_imp1k;

//number of the atoms with 2 impurity neighbor
int atom_imp2k;

//number of the atoms with 3 impurity neighbor
int atom_imp3k;

//number of the atoms with 1 impurity neighbor and 1 C neighbor
int atom_imp1k_1;

//number of the atoms with 1 impurity neighbor and C neighbors
int atom_imp1k_2;

//number of the atoms with 2 impurity neighbors and 1 C neighbor
int atom_imp2k_1;

/* hop of atoms with impuritites on diagonal */
//number of the atoms with 1 impurity neighbor
int atom_imp1kN;

//number of the atoms with 2 impurity neighbor
int atom_imp2kN;

//number of the atoms with 3 impurity neighbor
int atom_imp3kN;

//number of the atoms with 1 impurity neighbor and 1 C neighbor
int atom_imp1k_1N;

//number of the atoms with 1 impurity neighbor and C neighbors
int atom_imp1k_2N;

//number of the atoms with 2 impurity neighbors and 1 C neighbor
int atom_imp2k_1N;

//array of atoms with 1 impurity neioghbor
int atom_imp1[20000];

//array of atoms with 1 impurity neioghbor and 1 Si neighbor
int atom_imp1_1[20000];

//array of atoms with 1 impurity neioghbor and 2 Si neighbors
int atom_imp1_2[20000];

//array of atoms with 2 impurity neioghbors
int atom_imp2[20000];

//array of atoms with 2 impurity neioghbors
int atom_imp2_1[20000];

//array of atoms with 3 impurity neioghbors
int atom_imp3[20000];

/* arrays of atoms with impuritites that can hop on diagonal */
//array of atoms with 1 impurity neioghbor
int atom_imp1N[20000];

//array of atoms with 1 impurity neioghbor and 1 C neighbor
int atom_imp1_1N[20000];

//array of atoms with 1 impurity neioghbor and 2 C neighbors
int atom_imp1_2N[20000];

//array of atoms with 2 impurity neioghbors
int atom_imp2N[20000];

//array of atoms with 2 impurity neioghbors
int atom_imp2_1N[20000];

//array of atoms with 3 impurity neioghbors
int atom_imp3N[20000];


void ImpPreDeposition ( int ImpNumb )
{
  //number of current deposited impurity
  int ImpDeposited;
  int Depos_flaq;

  ImpDeposited = 0;
  while(ImpDeposited < ImpNumb)
  {
     Depos_flaq = Deposition(2);
     if( Depos_flaq != 0)
       ImpDeposited++;
  }
  
}

//calculation of the bonds before the hop to the nearest position with impurity neigbors
int ImpBonds( int process_state )
{
   int random_atom;
   int atomNumber;
   int dircount;
   

   dircount = 1;
   
   if(process_state == 9) 
   {
     random_atom = 1 + rand()%(atom_imp1k);
     atomNumber = atom_imp1[random_atom];
     NeighborList(atom_imp1[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     //printf("%d",newpos);
     //getch();
     return newpos;
     
    } 
    
   
   if(process_state == 10) 
   {
     random_atom = 1 + rand()%(atom_imp1k_1);
     atomNumber = atom_imp1_1[random_atom];
     NeighborList(atom_imp1_1[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     return newpos;
    } 
    
   if(process_state == 11) 
   {
     random_atom = 1 + rand()%(atom_imp1k_2);
     atomNumber = atom_imp1_2[random_atom];
     NeighborList(atom_imp1_2[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     return newpos;
    } 
    
    
   if(process_state == 12) 
   {
     random_atom = 1 + rand()%(atom_imp2k);
     atomNumber = atom_imp2[random_atom];
     NeighborList(atom_imp2[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     return newpos;
    } 
    
   
   if(process_state == 13) 
   {
     random_atom = 1 + rand()%(atom_imp2k_1);
     atomNumber = atom_imp2_1[random_atom];
     NeighborList(atom_imp2_1[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     return newpos;
    } 
    
   if(process_state == 14) 
   {
     random_atom = 1 + rand()%(atom_imp3k);
     atomNumber = atom_imp3[random_atom];
     NeighborList(atom_imp3[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[neighbor[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[neighbor[dirImp]] = 1;
         newpos = neighbor[dirImp];
       }
     }
     return newpos;
    } 
    
    ///dgffgdfg
   if(process_state == 15) 
   {
     random_atom = 1 + rand()%(atom_imp1kN);
     atomNumber = atom_imp1N[random_atom];
     NeighborList(atom_imp1N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
    
   
   if(process_state == 16) 
   {
     random_atom = 1 + rand()%(atom_imp1k_1N);
     atomNumber = atom_imp1_1N[random_atom];
     NeighborList(atom_imp1_1N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
    
   if(process_state == 17) 
   {
     random_atom = 1 + rand()%(atom_imp1k_2N);
     atomNumber = atom_imp1_2N[random_atom];
     NeighborList(atom_imp1_2N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
    
    
   if(process_state == 18) 
   {
     random_atom = 1 + rand()%(atom_imp2kN);
     atomNumber = atom_imp2N[random_atom];
     NeighborList(atom_imp2N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
    
   
   if(process_state == 19) 
   {
     random_atom = 1 + rand()%(atom_imp2k_1N);
     atomNumber = atom_imp2_1N[random_atom];
     NeighborList(atom_imp2_1N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
    
   if(process_state == 20) 
   {
     random_atom = 1 + rand()%(atom_imp3kN);
     atomNumber = atom_imp3N[random_atom];
     NeighborList(atom_imp3N[random_atom]);
     
     while(dircount != 2)
     {
       dirImp = 1 + rand()%4;
       if(adsite[Nextneigh[dirImp]] == 0)
       {
         dircount = 2;
         adsite[atomNumber] = 0;
         adsite[Nextneigh[dirImp]] = 1;
         newpos = Nextneigh[dirImp];
       }
     }
     return newpos;
    } 
        
    
        
}
