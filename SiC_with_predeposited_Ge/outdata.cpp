#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <graphics.h>
#include <math.h>
 

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"

FILE *output;

void out( void )
{
    int i, j, k;
   
   
   output = fopen("simout","w+");  
   totalatoms = 8 * xcell * ycell * zcell;
   //output of simulation box
   for(i=1; i<=totalatoms; i++)
      fprintf(output, "1  %f   %f   %f\n", toplayer[i][0], toplayer[i][1], toplayer[i][2]);  
   for(i=1; i<=totalvacpos; i++)
      fprintf(output, "5  %f   %f   %f\n", vacantpos[i][0], vacantpos[i][1], vacantpos[i][2]);  
   for(k = 1; k <= xcell * ycell * 2; k++)
   {
     if(adsite[k] == 1)
     {
        fprintf(output, "2  %f  %f  %f\n",adcoords[k][0],adcoords[k][1],adcoords[k][2]);
        printf("Hello %d  %f  %f  %f\n",k, adcoords[k][0],adcoords[k][1],adcoords[k][2]);
     }   
   }

  fclose(output);
}   


void out1( double x, double y, double z )
{
  
   
   fprintf(Diff, "2  %f  %f  %f\n",x,y,z);
   
  
}   


void DegreeR( double check, double degree )
{
  int i;
  
  double value;
  
   
  int r_int, r_int1, r_int2, r_int3;
   
   for(i = 0; i < 8; i++)
      timedata[i] = '0';
          
          
   value = pow(10,degree);
   r_int1 = (int)(check / value);
   timedata[0] = r_int1 + 48;
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 1.0);
   r_int1 = (int)(check / value);
   timedata[0] = r_int1 + 48;    
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 2.0);
   r_int1 = (int)(check / value);
   timedata[1] = r_int1 + 48;  
   check = check - ((double)r_int1) * value;
   
   ///
   value = pow(10,degree - 3.0);
   r_int1 = (int)(check / value);
   timedata[2] = r_int1 + 48;
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 4.0);
   r_int1 = (int)(check / value);
   timedata[3] = r_int1 + 48;    
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 5.0);
   r_int1 = (int)(check / value);
   timedata[4] = r_int1 + 48;  
   check = check - ((double)r_int1) * value;
   
   ///
   value = pow(10,degree - 6.0);
   r_int1 = (int)(check / value);
   timedata[5] = r_int1 + 48;
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 7.0);
   r_int1 = (int)(check / value);
   timedata[6] = r_int1 + 48;    
   check = check - ((double)r_int1) * value;
       
   value = pow(10,degree - 8.0);
   r_int1 = (int)(check / value);
   timedata[7] = r_int1 + 48;  
   check = check - ((double)r_int1) * value;
       
       
}

void OutStep( char timedata [8], double a[20000][3] )
{
   int i, j, k;
   char filename[11];
   //string filename;
   
   FILE *F;
   
   totalatoms = xcell * ycell * zcell * 8;
   for(i = 0; i < 3; i++)
     filename[i] = ' ';
     
   filename[0] = 'o';
   filename[1] = 'u';
   filename[2] = 't';
   
   for(i = 3; i < 11 ; i++)
     filename[i] = timedata[i-3]; 
     
   //printf("check");
   totalatoms = 8 * xcell * ycell * zcell;
   //F = *filename;
   
   F = fopen(filename,"w+");
   //printf("check");
   for(i=1; i<=totalatoms; i++)
      fprintf(F, "1  %f   %f   %f\n", toplayer[i][0], toplayer[i][1], toplayer[i][2]); 
       
   for(i = 0; i <= xcell*ycell*2; i++)
   { 
         
      if(adsite[i] == 1) 
      {   
        fprintf(F, "2  %f   %f   %f\n",adcoords[i][0],adcoords[i][1],adcoords[i][2] ); 
      }  
   }
     
   fclose(F);
}
