#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"
#include "additional.h"
#include "init.h"
#include "impurity.h"
#include "cluster.h"


void InitDone( void )
{
   int i, j, k;
   int read_flaq;
   std::fstream file;
   std::fstream file1, file2;
   
   read_flaq = 2;
   // adsite[i] = 0 - means free adatom site, adsite[i] = 1 - means occupied adatom site
   for(i = 1; i < 90000; i++)
   {
      adsite[i] = 0;
      Impur[i] = 0;
   }
   
   for(i = 1; i < 50000; i++)
   {
     adcoords[i][0] = 0;
     adcoords[i][1] = 0;
     adcoords[i][2] = 0;
   }
   
   
   if(read_flaq == 2)
   {
     file1.open("final",std::ios::in);
     file.open("final_copy",std::ios::out);
    
     //VacantPosCoords(xcell,ycell,zcell);
     j = 1;
     while(!file1.eof())
     { 
      //printf("test\n");
      //getch();              
       file1 >> i >> sys_adcoords[j][0] >> sys_adcoords[j][1] >> sys_adcoords[j][2];
      
      //printf("%d %d\n",j,i);
       if(i == 2 || i == 3)
       {
         for(k = 1; k <= xcell * ycell * 2; k++)
         {
           if(sys_adcoords[j][0] < vacantpos[k][0] + 0.00001 && sys_adcoords[j][0] > vacantpos[k][0] - 0.00001 && sys_adcoords[j][1] <vacantpos[k][1] + 0.00001 && sys_adcoords[j][1] > vacantpos[k][1] - 0.00001 && sys_adcoords[j][2] <vacantpos[k][2] + 0.00001 && sys_adcoords[j][2] > vacantpos[k][2] - 0.00001) 
           {  
             adsite[k] = 1;
             if(i == 2)
             {
               adcoords[k][0] = sys_adcoords[j][0];
               adcoords[k][1] = sys_adcoords[j][1];
               adcoords[k][2] = sys_adcoords[j][2];
               file << "2" << "   " << adcoords[k][0] << "   "<< adcoords[k][1]<<"   "<< adcoords[k][2] << "\n";
             }
             if(i == 3)
             {
               Impur[k] = 1;             
               adcoords[k][0] = sys_adcoords[j][0];
               adcoords[k][1] = sys_adcoords[j][1];
               adcoords[k][2] = sys_adcoords[j][2];
               file << "3" << "   " << adcoords[k][0] << "   "<< adcoords[k][1]<<"   "<< adcoords[k][2] << "\n";
             }
           } 
         }  
       }  
       j++;
    }   
    file.close();   
  }
  
}

/* checking random number generator */
void Random_Number_Generator_Test( void )
{
  int i,j;
  
  //variables for checking random number generator
  double rand_max = (double)RAND_MAX;
  double var_rand[100000];
  double var;
  double counts_rand[100000];
  int flaq;
  double sum;
  FILE *F;
 
   F = fopen("GeneratorTest","w+");
   //setting massive of counts of "random value" in diapazon [0,1]  to zero 
   for(i = 0; i < 100000; i++)
   {
     counts_rand[i] = 0.0;
     var_rand[i] = i / 100000.0;
   }
   //calculating random_value distribution function  
   for(i = 0; i < 1000; i++)
   {
     //printf("%d\n",i);
     var = rand()/rand_max;
     j = 0;
     flaq = 0;
     while( flaq != 1)  
     {
       if(var > var_rand[j])
         j++;
       if(var <= var_rand[j])
       { 
         counts_rand[j] = counts_rand[j] + 1.0;
         flaq = 1;
         //printf("j=%d\n",j);
       }          
     }
   }
   sum = 0.0;
   // normirovka distribution function
   for(i = 0; i < 100000; i++)
   {
      sum = sum + counts_rand[i];
      //printf("%f  %f\n",sum, counts_rand[i]);
   }  
   
   sum = sum * 1 / 100000.0;
   //printf("sum=%f\n",sum);
   for(i = 0; i < 100000; i++)
   {
     //counts_rand[i] = counts_rand[i] / sum;
     fprintf(F,"%f   %f\n",var_rand[i], counts_rand[i]);
   }  
   
   fclose(F);
} 
