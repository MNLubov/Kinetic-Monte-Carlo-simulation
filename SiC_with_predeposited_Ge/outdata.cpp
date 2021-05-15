#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include "graphics.h"
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"
#include "impurity.h"

FILE *output;

void out( void )
{
    int i, j, k;
   
   
   output = fopen("simout","w+");  
   totalatoms = 8 * xcell * ycell * zcell;
   //output of simulation box
   for(i=1; i<=totalatoms; i++)
      fprintf(output, "1  %f   %f   %f\n", toplayer[i][0], toplayer[i][1], toplayer[i][2]);  
   //for(i=1; i<=totalvacpos; i++)
     //fprintf(output, "5  %f   %f   %f\n", vacantpos[i][0], vacantpos[i][1], vacantpos[i][2]);  
   for(k = 1; k <= xcell * ycell * 2; k++)
   {
     if(adsite[k] == 1 && Impur[k] == 0)
        fprintf(output, "2  %f  %f  %f\n",adcoords[k][0],adcoords[k][1],adcoords[k][2]);
     if(adsite[k] == 1 && Impur[k] == 1)   
        fprintf(output, "3  %f  %f  %f\n",adcoords[k][0],adcoords[k][1],adcoords[k][2]);
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

// Output of the system configuration each 
void OutStep( void )
{
   int i, j, k;
   char filename[11];
     
   //FILE *FOutStep;
   std::fstream file;
   FILE *f1;
   
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
   
   file.open(filename,std::ios::out);
   f1 = fopen("f1","w+");
   /* output of  the substrate */
   
   for(i=1; i<=totalatoms; i++)
   {
      file<< "1   " << toplayer[i][0] << "  "<< toplayer[i][1] << "  "<< toplayer[i][2] << "\n";
      //fprintf(f1, "1  %f   %f   %f\n", c);  
      //fprintf(FOutStep, "1  %f   %f   %f\n", toplayer[i][0], toplayer[i][1], toplayer[i][2]);  
   }   
   
   
   for(i = 1; i <= xcell*ycell*2; i++)
   { 
         
      if(adsite[i] == 1 && Impur[i] == 1)
        file<< "3   " << adcoords[i][0] << "  "<< adcoords[i][1] << "  "<< adcoords[i][2] << "\n";
   }
   
   for(i = 1; i <= xcell*ycell*2; i++)
   { 
         
      if(adsite[i] == 1 && Impur[i] == 0) 
      {   
        file<< "2   " << adcoords[i][0] << "  "<< adcoords[i][1] << "  "<< adcoords[i][2] << "\n";
        //fprintf(f1, "2  %f   %f   %f\n",adcoords[i][0],adcoords[i][1],adcoords[i][2] ); 
        //fprintf(FOutStep, "2  %f   %f   %f\n",adcoords[i][0],adcoords[i][1],adcoords[i][2] ); 
      }  
     
   }
   
   
     
   //fclose(FOutStep);
   fclose(f1);
}

//vyvod ASCII symvolov ot 1 do 100
void AsciiOut ( void )
{

  int i;
  FILE *symbols;
  char c;
  
  symbols = fopen("symb","w+");
  for(i = 1; i < 100; i++)
  {
    c = i;
    fprintf(symbols,"%d  %c\n", i, c);
  }
  fclose(symbols);
}
