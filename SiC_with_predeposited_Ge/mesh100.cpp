#include <stdio.h>
#include <stdlib.h>
#include "mesh100.h"
#include "global.h"

//total amount of atoms in substrate box
int totalatoms;

//total amount of vacant positions
int totalvacpos;

//total amount of deposited atoms
int total;

//substrate box cell number
int xcell;
int ycell;
int zcell;

//first indicator atom number, second - coordinate 0-x-coord, 1-y-coord
//coords of the substrate box
double toplayer[80000][3];

//coords of vacant positions
double vacantpos[20000][3];

//coords of elementar box
double basis [9][3];

//coords of elementar vacant box 
double vacbasis[9][3];

void BasisCoords ( void )
{
   FILE *F;  
   int i,j;
     
   basis[0][0]= 0;
   basis[0][1]= 0;
   
   // element box of topo layer
   basis[1][0]= 0;
   basis[1][1]= 0;
   basis[1][2]= 0.0;
   basis[2][0]= 0.5;
   basis[2][1]= 0;
   basis[2][2]= 0.5;
   basis[3][0]= 0.25;
   basis[3][1]= 0.25;
   basis[3][2] = 0.25;
   basis[4][0]= 0.75;
   basis[4][1]= 0.25;
   basis[4][2]=0.75;
   basis[5][0]= 0;
   basis[5][1]= 0.5;
   basis[5][2]=0.5;
   basis[6][0]= 0.5;
   basis[6][1]= 0.5;
   basis[6][2]=0.0;
   basis[7][0]= 0.25;
   basis[7][1]= 0.75;
   basis[7][2]=0.75;
   basis[8][0]= 0.75;
   basis[8][1]= 0.75;
   basis[8][2]=0.25;
   
   /*for(i=0; i<9; i++)
     basis[i][2]= 0.0;*/
     
   
   F = fopen("elementbox","w+");
   for(i=1; i<9; i++)
      fprintf(F, "1  %f   %f   %f\n", basis[i][0], basis[i][1], basis[i][2]);
      
   fclose(F);
}


void MakeTopLayerCoords( int xcell, int ycell, int zcell )
{
   int i,j;
   int cell;
   int k,l,m;
   FILE *F;
   
   
   F = fopen("toplayer","w+");
   
   totalatoms = 8 * xcell * ycell * zcell;
   //fprintf(F, "%d   \n", totalatoms);  
   
   
   BasisCoords();
   
   cell = 1;
   for(k = 1; k <=xcell; k++)
     for(l = 1; l <= ycell; l++)
       for(m = 1; m <= zcell; m++)
       {
         for(i = 1 + 8 * (cell - 1); i <= 8 * cell; i++)
         {
           toplayer[i][0] =  basis[i - 8 * (cell-1)][0] + 1.0 * (k - 1); 
           //fprintf(F, "log xcell= %d  x= %f , i=%d  basis =%f\n", k,basis[i - 8 * (cell-1)][0],i - 8 * (cell-1),basis[i][0]);    
           toplayer[i][1] =  basis[i- 8 * (cell-1)][1] + 1.0 * (l - 1);
           //fprintf(F, "log zcell= %d  i= %d   \n", zcell,i);  
           toplayer[i][2] =  basis[i- 8 * (cell-1)][2] + 1.0 * (m - 1);
         } 
         cell++;
       }
  
   
   
     
   for(i = 1; i <= totalatoms; i++)
   {
       toplayer[i][0] =  toplayer[i][0] - (0.75 + (xcell-1)) / 2;   
       toplayer[i][1] =  toplayer[i][1] - (0.75 + (ycell-1)) / 2; 
       toplayer[i][2] =  toplayer[i][2] - (0.75 + (zcell-1)) / 2; 
   }    
   for(i=1; i<=totalatoms; i++)
      fprintf(F, "1  %f   %f   %f\n", toplayer[i][0], toplayer[i][1], toplayer[i][2]);  
      
   fclose(F);

}


void VacBasisCoords( void )
{
     
  FILE *F;   
  int i;
  // first vacant layer
  vacbasis[1][0] = 0.0;
  vacbasis[1][1] = 0.0;
  vacbasis[1][2] = 1.0;
  vacbasis[2][0] = 0.5;
  vacbasis[2][1] = 0.5;
  vacbasis[2][2] = 1.0;
  //second vacant layer
  vacbasis[3][0]= 0.25;
  vacbasis[3][1]= 0.25;
  vacbasis[3][2] = 1.25;
  vacbasis[4][0]= 0.75;
  vacbasis[4][1]= 0.75;
  vacbasis[4][2]=1.25;
  //third vacant layer
  vacbasis[5][0]= 0.5;
  vacbasis[5][1]= 0;
  vacbasis[5][2]= 1.5;
  vacbasis[6][0]= 0;
  vacbasis[6][1]= 0.5;
  vacbasis[6][2]=1.5;
  //fourth layer
  vacbasis[7][0]= 0.75;
  vacbasis[7][1]= 0.25;
  vacbasis[7][2]=1.75;
  vacbasis[8][0]= 0.25;
  vacbasis[8][1]= 0.75;
  vacbasis[8][2]=1.75;


 F = fopen("elementvacbox","w+");
   for(i=1; i<9; i++)
      fprintf(F, "1  %f   %f   %f\n", vacbasis[i][0], vacbasis[i][1], vacbasis[i][2]);
      
   fclose(F);
}

void VacantPosCoords(int xcell, int ycell, int zcell)
{
  int i,j;
  int cell;
  int k,l,m;
  FILE *F;
  
  
  totalvacpos = xcell*ycell*8;
  F = fopen("vacpos","w+");
  // zvaclayer = 1 means half of the initial box height
   
   VacBasisCoords();
  //F = fopen("toplayer","w+");
  cell = 1;
  
  //first vac layer
  for(k = 1; k <=xcell; k++)
     for(l = 1; l <= ycell; l++)
       {
         for(i = 1 + 2 * (cell - 1); i <= 2 * cell; i++)
         {
           vacantpos[i][0] =  vacbasis[i - 2 * (cell-1)][0] + 1.0 * (k - 1);   
           vacantpos[i][1] =  vacbasis[i- 2 * (cell-1)][1] + 1.0 * (l - 1);
           vacantpos[i][2] =  vacbasis[i- 2 * (cell-1)][2]+ zcell-1;
         } 
         cell++;
       }
  //second vac layer
  cell = 1;
  for(k = 1; k <=xcell; k++)
     for(l = 1; l <= ycell; l++)
       {
         for(i = 2 * xcell * ycell + 1 + 2 * (cell - 1); i <= 2 * xcell * ycell + 2 * cell; i++)
         {
           vacantpos[i][0] =  vacbasis[i + 2 - 2 * xcell * ycell- 2 * (cell-1)][0] + 1.0 * (k - 1);   
           vacantpos[i][1] =  vacbasis[i+ 2- 2 * xcell * ycell - 2 * (cell-1)][1] + 1.0 * (l - 1);
           vacantpos[i][2] =  vacbasis[i + 2 - 2 * xcell * ycell - 2 * (cell-1)][2]+ zcell-1;
         } 
         cell++;
       }
  
  //third vac layer
  cell = 1;
  for(k = 1; k <=xcell; k++)
     for(l = 1; l <= ycell; l++)
       {
         for(i = 4 * xcell * ycell + 1 + 2 * (cell - 1); i <= 4 * xcell * ycell + 2 * cell; i++)
         {
           vacantpos[i][0] =  vacbasis[i + 4 - 4 * xcell * ycell- 2 * (cell-1)][0] + 1.0 * (k - 1);   
           vacantpos[i][1] =  vacbasis[i+ 4- 4 * xcell * ycell - 2 * (cell-1)][1] + 1.0 * (l - 1);
           vacantpos[i][2] =  vacbasis[i + 4 - 4 * xcell * ycell - 2 * (cell-1)][2]+ zcell-1;
         } 
         cell++;
       }
       
       
  //fourth vac layer
  cell = 1;
  for(k = 1; k <=xcell; k++)
     for(l = 1; l <= ycell; l++)
       {
         for(i = 6 * xcell * ycell + 1 + 2 * (cell - 1); i <= 6 * xcell * ycell + 2 * cell; i++)
         {
           vacantpos[i][0] =  vacbasis[i + 6 - 6 * xcell * ycell- 2 * (cell-1)][0] + 1.0 * (k - 1);   
           vacantpos[i][1] =  vacbasis[i+ 6 - 6 * xcell * ycell - 2 * (cell-1)][1] + 1.0 * (l - 1);
           vacantpos[i][2] =  vacbasis[i + 6 - 6 * xcell * ycell - 2 * (cell-1)][2]+ zcell-1;
         } 
         cell++;
       }   
  
  for(i = 1; i <= totalvacpos; i++)
   {
       vacantpos[i][0] =  vacantpos[i][0] - (0.75 + (xcell-1)) / 2;   
       vacantpos[i][1] =  vacantpos[i][1] - (0.75 + (ycell-1)) / 2; 
       vacantpos[i][2] =  vacantpos[i][2] - (0.75 + (zcell-1)) / 2; 
       
   }    
   
   for(i=1; i<=totalvacpos; i++)
      fprintf(F, "%d  %f   %f   %f\n",2, vacantpos[i][0], vacantpos[i][1], vacantpos[i][2]);  
   
    fclose(F);
}
