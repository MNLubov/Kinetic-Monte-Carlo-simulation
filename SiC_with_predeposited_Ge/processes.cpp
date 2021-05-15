#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"
#include "outdata.h"
#include "impurity.h"

// index of occuring process (deposition, diffusion of adatom, diffusion of atom with 1 neighbor e.t.c)
int process;

//rates of the processes
//deposition rate
double DeposRate;
//adatom diffusuin rate
double AdatomDiff;
//atom with one neighbor rate
double Atom1Diff;
//atom with 2 neighbors rate
double Atom2Diff;
//atom with 3 neighbors rate
double Atom3Diff;
//total rate
double total_rate;


//diffusion data file pointer
FILE *Diff;

// opredelenie nomerov mest sosednih s dannym
void NeighborList( int atom_curr )
{
   int xcell_number, ycell_number, z;
   int i, j, k;
   
      
   cell_curr = 0;
   // calculation of cell number in global cell list, in x and y directions
   if(atom_curr % 2 != 0) 
     cell_curr = atom_curr / 2 + 1;
   else
     cell_curr = atom_curr / 2;
     
   if( cell_curr % ycell == 0)
   {
      ycell_number = ycell;
      
   }
   else
   {
      ycell_number = cell_curr % ycell;
      
   }
   if( cell_curr % xcell == 0)
   {
      xcell_number = cell_curr / xcell;
      
   }
   else
   {
      xcell_number = cell_curr / ycell + 1;
      
   }
   
   //postroenie of neighbor list
   //corner sites and edge sites have neighbors on the othe side of the simulation box 
   
   //left edge atoms neighbors, without corner atoms
   if( xcell_number == 1 && (atom_curr%2) != 0  && atom_curr != 1 && atom_curr != ycell * 2 - 1)
   {
     neighbor [1] = xcell * ycell * 2 - (ycell - ycell_number) * 2;
     neighbor [2] = atom_curr + 1;
     neighbor [3] = atom_curr - 1;
     neighbor [4] = xcell * ycell * 2 - (ycell - ycell_number) * 2 - 2;
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = xcell * ycell * 2 - (ycell - ycell_number) * 2 - 1;
     Nextneigh[2] = atom_curr + 2;
     Nextneigh[3] = atom_curr + ycell * 2;
     Nextneigh[4] = atom_curr - 2;
     //fprintf(F2, "%d   %d   %d  %d\n", Nextneigh [1], Nextneigh [2], Nextneigh [3], Nextneigh [4]);
   }
   
   
   //right edge atoms neighbors, without corners
   if( xcell_number == xcell && (atom_curr%2) == 0 && atom_curr != (2 * ycell * (xcell - 1) + 1) && atom_curr != ycell * xcell * 2)
   {
     neighbor [1] = atom_curr + 1; 
     neighbor [2] = ycell_number * 2 + 1;
     neighbor [3] = ycell_number * 2 - 1;
     neighbor [4] = atom_curr - 1; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = atom_curr - ycell * 2;
     if( atom_curr == xcell * ycell * 2 - 1)
        Nextneigh[2] = atom_curr - ycell * 2 + 2;        
     else 
        Nextneigh[2] = atom_curr + 2;
     Nextneigh[3] = ycell_number * 2;
     Nextneigh[4] = atom_curr - 2;
   }
      
   //top edge atoms neighbors, without corners
   if( ycell_number == ycell && (atom_curr%2) == 0 && atom_curr != (ycell * 2 - 1)  && atom_curr != ycell * xcell * 2)
   {
     neighbor [1] = ycell * (xcell_number - 1) * 2 + 1;
     neighbor [2] = ycell * xcell_number * 2 + 1;
     neighbor [3] = ycell * (xcell_number + 1) * 2 - 1;
     neighbor [4] = atom_curr - 1; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = atom_curr - ycell * 2;
     Nextneigh[2] = atom_curr - ycell * 2 + 2;
     Nextneigh[3] = atom_curr + ycell * 2;
     Nextneigh[4] = atom_curr - 2;
   }
   
   
   //buttom edge atoms neighbors, without corners
   if( ycell_number == 1 && (atom_curr%2) == 1 && atom_curr != 1  && atom_curr != (ycell * (xcell - 1) * 2 + 1))
   {
     neighbor [1] = ycell * (xcell_number - 2) * 2 + 2;
     neighbor [2] = atom_curr + 1;
     neighbor [3] = ycell * xcell_number * 2;
     neighbor [4] = ycell * (xcell_number - 1 ) * 2; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = atom_curr - ycell * 2;
     Nextneigh[2] = atom_curr + 2;
     Nextneigh[3] = atom_curr + ycell * 2;
     Nextneigh[4] = atom_curr +ycell * 2 - 2;
   }
   
   //left buttom corner atom neighbors
   if(atom_curr == 1)
   {
     neighbor [1] = ycell * (xcell-1) * 2 + 2; 
     neighbor [2] = 2;
     neighbor [3] = ycell * 2;
     neighbor [4] = xcell * ycell * 2; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = ycell * (xcell-1) * 2 + 1;
     Nextneigh[2] = atom_curr + 2;
     Nextneigh[3] = atom_curr + ycell * 2;
     Nextneigh[4] = atom_curr + ycell * 2 - 2;
   }
   
    //left top corner atom neighbors
    if(atom_curr == (ycell * 2 - 1))
    {
     neighbor [1] = ycell * (xcell-1) * 2 + 2; 
     neighbor [2] = ycell * 2;
     neighbor [3] = atom_curr - 1;
     neighbor [4] = xcell * ycell * 2 - 2; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
     //Next nearest neighbors
     Nextneigh[1] = ycell * xcell * 2 - 1;
     Nextneigh[2] = 1;
     Nextneigh[3] = atom_curr + ycell * 2;
     Nextneigh[4] = atom_curr - 2;
    }
    
    //right buttom corner atom neighbors
    if(atom_curr == (( xcell - 1) * ycell * 2 + 1) )
    {
     neighbor [1] = ycell * (xcell-2) * 2 + 2; 
     neighbor [2] = atom_curr + 1;
     neighbor [3] = xcell * ycell * 2;
     neighbor [4] = (xcell-1) * ycell * 2; 
     //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
          
     //Next nearest neighbors
     Nextneigh[1] = atom_curr - ycell * 2;
     Nextneigh[2] = atom_curr + 2;
     Nextneigh[3] = 2;
     Nextneigh[4] = xcell * ycell * 2;
    }
    
    
    //right top corner atom neighbors
    if(atom_curr == ycell * xcell * 2 )
    {
      neighbor [1] = ycell * (xcell-1) * 2 + 2; 
      neighbor [2] = 1;
      neighbor [3] = ycell * 2 - 1;
      neighbor [4] = ycell * xcell * 2 - 1; 
      //fprintf(F2, "%d   %d   %d  %d\n", neighbor [1], neighbor [2], neighbor [3], neighbor [4]);
     
      //Next nearest neighbors
      Nextneigh[1] = atom_curr - ycell * 2;
      Nextneigh[2] = xcell * ycell * 2 - ycell * 2 + 2;
      Nextneigh[3] = ycell * 2;
      Nextneigh[4] = xcell * ycell * 2 - 2;
    }
   
   //central part atoms, with chetnie numbers
   if(atom_curr % 2 == 0 && xcell_number != xcell && ycell_number != ycell)
   {
     neighbor [1] = atom_curr + 1; 
     neighbor [2] = ycell * xcell_number * 2 + ycell_number * 2 + 1;
     neighbor [3] = ycell * xcell_number * 2 + ycell_number * 2 - 1 ;
     neighbor [4] = atom_curr -1; 
     
     
     //Next nearest neighbors
     if(xcell_number == 1)
     {
       Nextneigh[1] = ycell * (xcell - 1) * 2 + ycell_number * 2;
       Nextneigh[2] = atom_curr + 2;
       Nextneigh[3] = atom_curr + ycell * 2;
       if(atom_curr == 2)
         Nextneigh[4] = ycell * 2;
       else  
         Nextneigh[4] = atom_curr - 2;
     } 
      
     if(ycell_number == 1)
     {
       if(atom_curr == 2)
         Nextneigh[1] = ycell * 2;
       else  
         Nextneigh[1] = atom_curr - ycell * 2;
       Nextneigh[2] = atom_curr + 2;
       Nextneigh[3] = atom_curr + ycell * 2;    
       Nextneigh[4] = atom_curr +ycell * 2 - 2;
     } 
      
      if( xcell_number != 1 && ycell_number != 1)
      {
         Nextneigh[1] = atom_curr - ycell * 2;
         Nextneigh[2] = atom_curr + 2;
         Nextneigh[3] = atom_curr + ycell * 2;
         Nextneigh[4] = atom_curr - 2;     
      }
      
   }
   
   //central part atoms, with nechetnie numbers
   if(atom_curr % 2 == 1 && ycell_number != 1 && xcell_number != 1)
   {
     neighbor [1] = ycell * (xcell_number - 2) * 2 + ycell_number * 2; 
     neighbor [2] = atom_curr + 1;
     neighbor [3] = atom_curr - 1 ;
     neighbor [4] = ycell * (xcell_number - 2) * 2 + (ycell_number - 1) * 2; 
     
      //Next nearest neighbors
     if(ycell_number == ycell)
     {
       Nextneigh[1] = atom_curr - ycell * 2;
       if( atom_curr == xcell * ycell * 2 - 1)
         Nextneigh[2] = atom_curr - ycell * 2 + 2;
       else
         Nextneigh[2] = atom_curr - ycell * 2 + 2;
       Nextneigh[3] = atom_curr + ycell * 2;
       Nextneigh[4] = atom_curr - 2;
       
      } 
      
     if(xcell_number == xcell)
     {
       Nextneigh[1] = atom_curr - ycell * 2;
       if( atom_curr == xcell * ycell * 2 - 1)
         Nextneigh[2] = atom_curr - ycell * 2 + 2;        
       if( atom_curr != xcell * ycell - 1) 
         Nextneigh[2] = atom_curr + 2;
       Nextneigh[3] = ycell_number * 2 - 1;    
       Nextneigh[4] = atom_curr - 2;
       //fprintf(ch,"%d  %d  %d  %d\n",Nextneigh[1], Nextneigh[2],Nextneigh[3],Nextneigh[4]);
      } 
      if( xcell_number != xcell && ycell_number != ycell)
      {
         Nextneigh[1] = atom_curr - ycell * 2;
         Nextneigh[2] = atom_curr + 2;
         Nextneigh[3] = atom_curr + ycell * 2;
         Nextneigh[4] = atom_curr - 2;
     
      }
   }    
   //right edge atoms neares neighbors
   if(atom_curr == xcell * ycell * 2 - 1)
   {
       Nextneigh[1] = atom_curr - ycell * 2;
       Nextneigh[2] = atom_curr - ycell * 2 + 2;
       Nextneigh[3] = ycell_number * 2 - 1;    
       Nextneigh[4] = atom_curr - 2;
   }
   
   if(atom_curr == ycell * 2)
   {
       Nextneigh[1] = xcell * ycell * 2;
       Nextneigh[2] = atom_curr - ycell * 2 + 2;
       Nextneigh[3] = atom_curr + ycell * 2;    
       Nextneigh[4] = atom_curr - 2;
   }
}


//deposition to first layer
int Deposition( int sort )
{
   int i, j, k;
   int atom_number;
   int numb;
   

   double totd;
   
   // total number of sites in the first vacant layer 
   numb = xcell*ycell*2;
   
  
   /* deposition */
   atom_number = 1 + rand()% numb; 
    
   if( adsite[atom_number] == 0 )// && vacantpos[atom_number][2] < ((zcell + 0.25)/2 + 0.0001) )
   {
      adcoords[atom_number][0] = vacantpos[atom_number][0];
      adcoords[atom_number][1] = vacantpos[atom_number][1];
      adcoords[atom_number][2] = vacantpos[atom_number][2];
      //fprintf(F, "2   %f  %f  %f\n", vacantpos[atom_number][0],vacantpos[atom_number][1],vacantpos[atom_number][2]);
      //fprintf(F1, "2   %f  %f  %f\n", deposcoords[atom_number][0],deposcoords[atom_number][1],deposcoords[atom_number][2]);
      //deposited[atom_number] = 1;
      adsite[atom_number] = 1;
      atom_counter++;
      if( sort == 2 )
      {
        Impur[atom_number] = 1;
        //printf("Impur=%d\n",atom_number);
      }  
      //if(sort == 1)
       // printf("%d  %d\n",adsite[atom_number],atom_number);
      //printf("atom_number=%d\n",atom_number);
      return atom_number;      
   } 

   if(adsite[atom_number] == 1)
     return 0;
  
} 


//calculation of the diffusion process
int Diffusion( int process_state, int diff_sort )
{
  int i, j, k;
  
  int dircount = 0;
  //number of the position of the atom before the hop
  int DiffAnumber;
  //number of the position of the atom after the hop
  int NewAnumber;
  //direction of diffusion for atom with impurity neighbor
  int ImpDiff;
  
    
  //number ofneighbors before and after the hop
  int sumBefore, sumAfter1, sumAfter2;
   
   //nahozhdenie svobodnyh napravlenij
   //updating adatom numbers of occupied adatom sites
   //adat_curr[ad_curr] = neighbor[direction];
   if(process_state >= 2 && process_state <= 8)
   {
       DiffAnumber = BondsNearest(process_state);   
         
       //choosing direction
       directionDiff = 1 + rand()%jDir; //random number from 1 to 4;
       //printf("directionDiff=%d\n",directionDiff);
       //updating the global atom list
       adsite[DiffAnumber] = 0; 
       k = Diffdir[directionDiff];
       //printf("k=%d\n",Diffdir[directionDiff]);
       
       if(process_state == 2 || process_state == 4 || process_state == 6 || process_state == 8) 
       {
         adsite[BhopG[k]] = 1; 
         adcoords[BhopG[k]][0] = vacantpos[BhopG[k]][0];
         adcoords[BhopG[k]][1] = vacantpos[BhopG[k]][1];
         adcoords[BhopG[k]][2] = vacantpos[BhopG[k]][2];
         return BhopG[k];
       }
       
       if(process_state == 3 || process_state == 5 || process_state == 7) 
       {
         adsite[BhopNG[k]] = 1; 
         adcoords[BhopNG[k]][0] = vacantpos[BhopNG[k]][0];
         adcoords[BhopNG[k]][1] = vacantpos[BhopNG[k]][1];
         adcoords[BhopNG[k]][2] = vacantpos[BhopNG[k]][2];
         return BhopNG[k];     
       }  
   }
  
   if(process_state > 8 && process_state <= 20)
   {
    
     ImpDiff = ImpBonds( process_state );
     //printf("ImpDiff=%d\n",ImpDiff);
     //getch();
     adcoords[ImpDiff][0] = vacantpos[ImpDiff][0];
     adcoords[ImpDiff][1] = vacantpos[ImpDiff][1];
     adcoords[ImpDiff][2] = vacantpos[ImpDiff][2];
     return ImpDiff;  
   } 
 
  
 
}


