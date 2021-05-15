#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>



#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"

//defing physical constants
#define KbT0  0.025
#define T0  300.0
#define Debay  1.0E13

//define energies
#define Ea 1.65
#define EaDiag  2.8
#define EaN  0.3
#define EaNN  0.07

int AdsiteSum ( void )
{
   int i;
   
   int sum = 0;
   
   
   for(i = 0; i <= xcell * ycell * 2; i++)
     sum = sum + adsite[i];
     
   return sum;
}

//calculation of the bonds before and after the hop yo the nearest position
int BondsNearest ( int process_state )
{
  int i, j, k;
  int random_atom;
  int atomNumber;
  
  sumG = 0;
  sumNG = 0; 
  
    
  /*
  hops to the nearest positions
  */
  //hop with the same number of bonds before and after the hop
  for(i = 0; i <= 4; i++)
    direction[i] = 0;
  
  
  if(process_state == 2) 
  {
     random_atom = 1 + rand()%(adatom_k);
     atomNumber = Adatom[random_atom];
     NeighborList(Adatom[random_atom]);
     //calculation of number of nearest neighbors before the hop
     for(j = 1; j <= 4; j++)
     {
       sumG = sumG + adsite[neighbor[j]];
       BhopG[j] = neighbor[j];
     }
     // calculation of the nearest neighbors after the hop
      for(j = 1; j <= 4; j++)
      {  
        sumAG[j] = 0;
        NeighborList(BhopG[j]);
        for(k = 1; k <= 4; k++)
          sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is neccesary to include initial position of the hopping atom
        sumAG[j] = sumAG[j] - 1;
      }
      //hop of the adatom
      if(sumG == 0)
      {
        for(j = 1; j <= 4; j++)
          direction[j] = j; 
      }   
      //hop of the atom with the same number ofneighbors before and after the hop
      if(sumG != 0) 
      {
        for(j = 1; j <= 4; j++)
       {
         if((sumG - sumAG[j]) == 0 && adsite[BhopG[j]] == 0)
           direction[j] = j; 
         //else 
           //direction[j] = 0;
       }   
      }
     
  }
  
  //hop to the nearest position with breaking one bond
  if(process_state == 4)
  { 
    random_atom = 1 + rand()%(atom1_k);
    atomNumber = Atom1[random_atom];
    NeighborList(Atom1[random_atom]); 
    //testing the program
    printf("atomNumber = %d, atom1_k = %d\n",atomNumber, atom1_k); 
    //calculation of number of nearest neighbors before the hop
    for(j = 1; j <= 4; j++)
    {
      sumG = sumG + adsite[neighbor[j]];
      BhopG[j] = neighbor[j];
    }
    // calculation of the nearest neighbors after the hop
    for(j = 1; j <= 4; j++)
    {  
      sumAG[j] = 0;
      NeighborList(BhopG[j]);
      for(k = 1; k <= 4; k++)
        sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      sumAG[j] = sumAG[j] - 1;
      if((sumG - sumAG[j]) == 1 && adsite[BhopG[j]] == 0)
        direction[j] = j;
      //if((sumG-sumAG[j]) != 1 && (sumG-sumAG[j]) > 0)    
        //direction[j] = 0;  
      if(sumG == 1 && (sumG - sumAG[j] < 0) && adsite[BhopG[j]] == 0)
        direction[j] = j;
      
    }
  }
  
  //hop to the nearest position with breaking 2 bonds
  if(process_state == 6)
  { 
    random_atom = 1 + rand()%(atom2_k); 
    atomNumber = Atom2[random_atom];                  
    NeighborList(Atom2[random_atom]);  
    //calculation of number of nearest neighbors before the hop
    for(j = 1; j <= 4; j++)
    {
      sumG = sumG + adsite[neighbor[j]];
      BhopG[j] = neighbor[j];
    }
    // calculation of the nearest neighbors after the hop
    for(j = 1; j <= 4; j++)
    {  
      sumAG[j] = 0;
      NeighborList(BhopG[j]);
      for(k = 1; k <= 4; k++)
        sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      sumAG[j] = sumAG[j] - 1;
      if((sumG - sumAG[j]) == 2 && adsite[BhopG[j]] == 0)
        direction[j] = j;
      if(sumG == 2 && (sumG - sumAG[j] < 0) && adsite[BhopG[j]] == 0)
        direction[j] = j;
      //if((sumG - sumAG[j]) != 2 && (sumG - sumAG[j]) > 0)     
        //direction[j] = 0;
    }
  }
  
  //hop to the nearest position with breaking 3 bonds  
  if(process_state == 8)
  { 
    random_atom = 1 + rand()%(atom3_k); 
    atomNumber = Atom3[random_atom];                     
    NeighborList(Atom3[random_atom]);  
    //calculation of number of nearest neighbors before the hop
    for(j = 1; j <= 4; j++)
    {
      sumG = sumG + adsite[neighbor[j]];
      BhopG[j] = neighbor[j];
    }
    // calculation of the nearest neighbors after the hop
    for(j = 1; j <= 4; j++)
    {  
      sumAG[j] = 0;
      NeighborList(BhopG[j]);
      for(k = 1; k <= 4; k++)
        sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      sumAG[j] = sumAG[j] - 1;
      if((sumG - sumAG[j]) == 3 && adsite[BhopG[j]] == 0)
        direction[j] = j;
      //if(sumG != 3)     
       // direction[j] = 0;
    }
  }

  /*
  hops to the next nearest positions
  */
  //hop to the next nearest positions with the same number of the bonds before and after the hop
  if(process_state == 3) 
  {
     random_atom = 1 + rand()%(adatom_kN);    
     atomNumber = AdatomN[random_atom];                   
     NeighborList(AdatomN[random_atom]);
     //calculation of number of nearest neighbors before the hop
     for(j = 1; j <= 4; j++)
     {
       sumG = sumG + adsite[neighbor[j]];
       //sumNG = sumNG + adsite[Neighnext[j]];
       BhopG[j] = neighbor[j];
       BhopNG[j] = Nextneigh[j];
     }
     // calculation of the nearest neighbors after the hop
      for(j = 1; j <= 4; j++)
      {  
        sumAG[j] = 0;
        NeighborList(BhopNG[j]);
        for(k = 1; k <= 4; k++)
          sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is neccesary to include initial position of the hopping atom
      }
      //hop of the adatom
      if(sumG == 0)
      {
        for(j = 1; j <= 4; j++)
        {
          if(adsite[BhopNG[j]] == 0)
            direction[j] = j; 
          //else  
           // direction[j] = 0;           
        }  
      }   
      //hop of the atom with the same number ofneighbors before and after the hop
      if(sumG != 0) 
      {
        for(j = 1; j <= 4; j++)
        {
         if((sumG - sumAG[j]) == 0 && adsite[BhopNG[j]] == 0)
           direction[j] = j; 
         //else 
          // direction[j] = 0;
       }   
      }
  }
  
  //hop to the nearest position with breaking 1 bond
  if(process_state == 5)
  { 
    random_atom = 1 + rand()%(atom1_kN);     
    atomNumber = Atom1N[random_atom];               
    NeighborList(Atom1N[random_atom]);  
    //calculation of number of nearest neighbors before the hop
    for(j = 1; j <= 4; j++)
    {
      sumG = sumG + adsite[neighbor[j]];
      BhopG[j] = neighbor[j];
      BhopNG[j] = Nextneigh[j];
    }
    // calculation of the nearest neighbors after the hop
    for(j = 1; j <= 4; j++)
    {  
      sumAG[j] = 0;
      NeighborList(BhopNG[j]);
      for(k = 1; k <= 4; k++)
        sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      if((sumG - sumAG[j]) == 1 && adsite[BhopNG[j]] == 0)
        direction[j] = j;
      if(sumG == 1 && (sumG - sumAG[j]) < 0 && adsite[BhopNG[j]] == 0)
        direction[j] = j;
      //if(sumG != 1)    
        //direction[j] = 0;
    }
  }
  
  //hop to the nearest position with breaking 2 bonds
  if(process_state == 7)
  { 
    random_atom = 1 + rand()%(atom2_kN); 
    atomNumber = Atom2N[random_atom];                
    NeighborList(Atom2N[random_atom]);  
    //calculation of number of nearest neighbors before the hop
    for(j = 1; j <= 4; j++)
    {
      sumG = sumG + adsite[neighbor[j]];
      BhopG[j] = neighbor[j];
      BhopNG[j] = Nextneigh[j];
    }
    // calculation of the nearest neighbors after the hop
    for(j = 1; j <= 4; j++)
    {  
      sumAG[j] = 0;
      NeighborList(BhopNG[j]);
      for(k = 1; k <= 4; k++)
        sumAG[j] = sumAG[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      if((sumG - sumAG[j]) == 2 && adsite[BhopNG[j]] == 0)
        direction[j] = j;
      if(sumG == 2 && (sumG - sumAG[j]) < 0 && adsite[BhopNG[j]] == 0)
        direction[j] = j;
      //if(sumG != 2)     
      //  direction[j] = 0;
    }
  }
    
  i = 1;
  jDir = 0;
  while(i <= 4)
  {
    if(direction[i] == i)
    {
      jDir++;
      Diffdir[jDir] = i;
    }
    i++;
  }   
  printf("jDir=%d\n",jDir);
  return atomNumber;
}

//postroenie spiska jacheek, zanjatyh atomami
void GroupList( void )
{
  int i, j, k;
  //number of the nearest neighbors of current atom before the hop
  int sum;
  //number of the next nearest neighbors of current atom 
  int sumN;
  // number of the nearest neighbors of current after before the hop
  int sumA[5]; 
  // number of the next nearest neighbors of current atom after the hop
  int sumAN[5]; 
  //number of the nearest neighbors before the hop
  int Bhop[5];
  //number of the next nearest neighbors before the hop
  int BhopN[5];
  int sumDiff;
  
  //uniqflaq[1] - same number of bonds or adatom
  //uniqflaq[2] - break 1 bond, nearest neighbor
  //uniqflaq[3] - break 2 bonds, nearest neighbor
  //uniqflaq[4] - break 3 bonds, nearest neighbor
  //uniqflaq[5] - diagonal hops with same number of neighbors or adatom
  //uniqflaq[6] - diagonal hop, break 1 bond
  //uniqflaq[7] - diagonal hop, break 2 bonds
  //uniqflaq[8] - diagonal hop, break 3 bonds
  int uniqflaq[9];
   
  /*
  Bhop[0] = 0;
  BhopN[0] = 0;
  */
  
   
  for(i = 0; i <20000; i++)
   {
     

     Adatom[i] = 0;
     Atom1[i] = 0;
     Atom2[i] = 0;
     Atom3[i] = 0;
     AdatomN[i] = 0;
     Atom1N[i] = 0;
     Atom2N[i] = 0;
     Atom3N[i] = 0;
     //Aneigh[i] = 0;
     //ANneigh[i] = 0;
   }
  //inicializacija summy chisla sosedej atoma
  grouplist_count++;
  adatom_k = 0;
  atom1_k = 0;
  atom2_k = 0;
  atom3_k = 0;
  adatom_kN = 0;
  atom1_kN = 0;
  atom2_kN = 0;
  atom3_kN = 0;
  
  for(i = 1; i <= xcell * ycell * 2; i++)
  {
    for(j = 0; j <= 8; j++)
       uniqflaq[j] = 1;
    NeighborList( i );
    //printf("%d  %d  %d  %d",);
    sum = 0;
    sumN = 0;
    for(j = 1; j <= 4; j++)
    {
      sum = sum + adsite[neighbor[j]];
      sumN = sumN + adsite[Nextneigh[j]];  
      Bhop[j] = neighbor[j];
      BhopN[j] = Nextneigh[j];
      //printf("i = %d %d\n",i, neighbor[j]);
    }
        
    for(j = 1; j <= 4; j++)
    {  
      sumA[j] = 0;
      // calculating of the nearest neighbors after the hop
      NeighborList(Bhop[j]);
      for(k = 1; k <= 4; k++)
        sumA[j] = sumA[j] + adsite[neighbor[k]];   //"-1" because it is initial position of the hopping atom
      sumA[j] = sumA[j] - 1;
      sumAN[j] = 0;  
      NeighborList(BhopN[j]);  
      for(k = 1; k <= 4; k++)
        sumAN[j] = sumAN[j] + adsite[neighbor[k]];   
      //if(i == ngh)
       // printf("sumA[j]=%d\n",sumA[j]);  
    }
    
    
   
    /* claculating of the total number of the adatoms (or atoms) with the same number of tne bond before and after the hop */
    if(adsite[i] == 1 && sum == 0 && uniqflaq[1] != 2)
    {
      adatom_k++;
      Adatom[adatom_k] = i;
      uniqflaq[1] = 2;
    }
    if(adsite[i] == 1 && sum == 0 && sumN != 4 && uniqflaq[5] != 2)
    {
      adatom_kN++;
      AdatomN[adatom_kN] = i;
      uniqflaq[5] = 2;
      //printf("uniqflaq[5]=%d\n",uniqflaq[5]);
      //printf("ngh = %d\n",ngh);
    }
    
      
    for(k = 1; k <=4; k++)
    {
      sumDiff = sum - sumA[k];       
      if(adsite[i] == 1 && sumDiff == 0 && adsite[Bhop[k]] == 0 && uniqflaq[1] != 2)
      {
         adatom_k++;
         Adatom[adatom_k] = i;
         uniqflaq[1] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 1 && adsite[Bhop[k]] == 0 && uniqflaq[2] != 2)
      {
        atom1_k++;
        Atom1[atom1_k] = i;
        uniqflaq[2] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 2 && adsite[Bhop[k]] == 0 && uniqflaq[3] != 2)
      {
        atom2_k++;
        Atom2[atom2_k] = i;
        uniqflaq[3] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 3 && adsite[Bhop[k]] == 0 && uniqflaq[4] != 2)
      {
        atom3_k++;
        Atom3[atom3_k] = i;
        uniqflaq[4] = 2;
      } 
      
      if(adsite[i] == 1 && sumDiff < 0 && adsite[Bhop[k]] == 0)
      {
         if(sum == 1 && uniqflaq[2] != 2)
         {
           atom1_k++;
           Atom1[atom1_k] = i;
           uniqflaq[1] = 2;
         }
         if(sum == 2 && uniqflaq[3] != 2)
         {
           atom2_k++;
           Atom2[atom2_k] = i;
           uniqflaq[2] = 2;
         }
      }
      
      
    }
    
    /* claculating of the total number of the adatoms (or atoms) with the different number of the bonds before and after the hop */
    // if the number of the bonds after is higher the before the hop it ok
    //else calculate number of such atoms
    for(k = 1; k <=4; k++)
    {
      sumDiff = sum - sumAN[k];    
      if(adsite[i] == 1 && sumDiff == 0 && adsite[BhopN[k]] == 0 && uniqflaq[5] != 2)
      {
        adatom_kN++;
        AdatomN[adatom_kN] = i;
        uniqflaq[5] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 1 && adsite[BhopN[k]] == 0 && uniqflaq[6] != 2)
      {
        atom1_kN++;
        Atom1N[atom1_kN] = i;
        uniqflaq[6] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 2 && adsite[BhopN[k]] == 0 && uniqflaq[7] != 2)
      {
        atom2_kN++;
        Atom2N[atom2_kN] = i;
        uniqflaq[7] = 2;
      }
      if(adsite[i] == 1 && sumDiff == 3 && adsite[BhopN[k]] == 0 && uniqflaq[8] != 2)
      {
        atom3_kN++;
        Atom3N[atom3_kN] = i;
        uniqflaq[8] = 2;
      }
      if(adsite[i] == 1 && sumDiff < 0 && adsite[BhopN[k]] == 0)
      {
         if(sum == 1 && uniqflaq[6] != 2)
         {
           atom1_kN++;
           Atom1N[atom1_kN] = i;
           uniqflaq[6] = 2;
         }
         if(sum == 2 && uniqflaq[7] != 2)
         {
           atom2_kN++;
           Atom2N[atom2_kN] = i;
           uniqflaq[7] = 2;
         }
         
      }
    }
    //printf("i = %d\n",i);
  }
  //printf("adatom_k=%d\n",adatom_k);
  //delay(500);
}


//calculation of processes rates
int Rates( void )
{
   int i, j;
   double depos_flux;
   double T, alpha;
   double speed; //in ML/sec
   double prob;
   int probal;
   int adsites;
   FILE *Frates;
  
  
   double rand_max = (double)RAND_MAX;
   
   Frates = fopen("Rates","w+");
     
   T = 900;
   alpha = KbT0 * (T / T0);
   speed = 0.75;
   total_rate = 0;
   
   // rates of hops without keeping the number of neighbors including adatom hops rate 
   double r, r0, r1, r2, r3;
   // rates of hops with keeping the same number of nearest neighbors
   double r0N, r1N, r2N, r3N;


      
   
   //calculation of deposition flux in atoms per site in second
   depos_flux = deposition / sites_on_cm2;
   
   
   /* calculation of processes rates */
   
   /* calculation of processes rates of hops without keeping the number of neighbors including adatom hops rate */
   //rate of the deposition, i.e. number of atoms depositing on the simulation box in a second
   DeposRate = speed * xcell * ycell * 2;
   r = DeposRate;
   
   //rate of diffusion to the nearest positions with the same number of neighbors
   AdatomDiff = Debay * exp((-1 * Ea) / alpha );
   r0 = AdatomDiff;
   //rate of diffusion to the next nearest positions with the same number of neighbors
   r0N = Debay * exp((-1 * (EaDiag)) / alpha );
   
   //rate of diffusion to the nearest positions with breaking 1 bond
   Atom1Diff = Debay * exp((-1 * (Ea + EaN)) / alpha );
   r1 = Atom1Diff;
   //rate of diffusion to the next nearest positions with breaking 1 bond
   r1N = Debay * exp((-1 * (EaDiag + EaN)) / alpha );
   
   //rate of diffusion to the nearest positions with breaking 2 bonds
   Atom2Diff = Debay * exp((-1 * (Ea + 2 * EaN)) / alpha );
   r2 = Atom2Diff;
   //rate of diffusion to the next nearest positions with breaking 2 bonds
   r2N = Debay * exp((-1 * (EaDiag + 2 * EaN)) / alpha );
   
   //rate of diffusion to the nearest positions with breaking 3 bonds
   Atom3Diff = Debay * exp((-1 * (Ea + 3 * EaN)) / alpha );
   r3 = Atom3Diff;
   //rate of atom diffusion with 3 bonds to the next nearest positions is forbidden (to much energy to brake)
   /*r3N = Debay * exp((-1 * (EaDiag + 3 * EaN)) / alpha );*/
   r3N = 0;
  
     
   /*fprintf(F,"DeposRate=%f\n",DeposRate);
   fprintf(F,"AdatomDiff=%f\n",AdatomDiff);
   fprintf(F,"Atom1Diff=%f\n",Atom1Diff);
   fprintf(F,"Atom2Diff=%f\n",Atom2Diff);
   fprintf(F,"Atom4Diff=%f\n",Atom3Diff);
   fprintf(F,"deposition=%f\n",deposition);
   fprintf(F,"Debay=%f\n",Debay);
   fprintf(F,"xcell*ycell*2=%d\n",xcell*ycell*2);
   fprintf(F,"exp((-1 * Ea) / alpha )=%f\n",exp((-1 * Ea) / alpha ));
   */
   GroupList();
   
   //calculating of rates ataking into account number of the neighbors
   //r = DeposRate; // / total_rate;
   r0 = adatom_k * r0;
   r0N = adatom_kN * r0N;
   r1 = atom1_k * r1;
   r1N = atom1_kN * r1N;
   r2 = atom2_k * r2;
   r2N = atom2_kN * r2N;
   r3 = atom3_k * r3;
   r3N = atom3_k * r3N;
   
   total_rate = total_rate + r + r0 + r0N + r1 + r1N + r2 + r2N + r3;
   
   //calculating of the normalized rate
   r = r / total_rate;
   r0 = r0 / total_rate;
   r0N =  r0N / total_rate;
   r1 =  r1 / total_rate;
   r1N =  r1N / total_rate;
   r2 =  r2 / total_rate;
   r2N =  r2N / total_rate;
   r3 =  r3 / total_rate;
   r3N = 0.0;
   
   /*
   printf("r=%1.10f\n",r);
   printf("r0=%1.15f, adatom_k = %d\n",r0, adatom_k);
   printf("r0N=%1.10f, adatom_kN = %d\n",r0N, adatom_kN);
   */
   //printf("r1=%f\n",r1);
    
   /*
   fprintf(F,"AdatomDiff=%f\n",AdatomDiff);
   fprintf(F,"Atom1Diff=%f\n",Atom1Diff);
   fprintf(F,"Atom2Diff=%.10f\n",Atom2Diff);
   fprintf(F,"Atom3Diff=%.10f\n",Atom3Diff);
   fprintf(F,"deposition=%f\n",deposition);
   fprintf(F,"Debay=%f\n",Debay);
   fprintf(F,"xcell*ycell*2=%d\n",xcell*ycell*2);
   fprintf(F,"exp((-1 * Ea) / alpha )=%f\n",exp((-1 * Ea) / alpha ));
   */  
   
   // generating of the random number in the interval from [0,1]
   probal = rand() % RAND_MAX;
   prob = probal / rand_max;
   
   
   cc = cc + 1;
   fclose(Frates);
   if( prob >= 0 && prob <= r)
   {
     sum_check = sum_check + 1;
     adsites = AdsiteSum();
     //fprintf(RateData,"%d   %d  %d %.10f\n",cc,sum_check, adsites, DeposRate/(AdatomDiff+1.0));
     //printf("%d   %d  %d %.10f\n",cc,sum_check, adsites, DeposRate/(AdatomDiff+1.0));
     return 1;
   }
   
   //diffusion of atom with thesame number of the neghbors
   if( prob > r && prob <= (r + r0))
     return 2;
   if( prob > (r + r0) && prob <= (r + r0 + r0N))
     return 3;
   if( prob > (r + r0 + r0N) && prob <= (r + r0 + r0N + r1))
     return 4;
   if( prob > (r + r0 + r0N + r1) && prob <= (r + r0 + r0N + r1 + r1N))
     return 5;
   if( prob > (r + r0 + r0N + r1 + r1N) && prob <= (r + r0 + r0N + r1 + r1N + r2))
     return 6;  
   if( prob > (r + r0 + r0N + r1 + r2) && prob <= (r + r0 + r0N + r1 + r1N + r2 + r2N))
     return 7;
   if( prob > (r + r0 + r0N + r1 + r1N + r2 + r2N) && prob <= (r + r0 + r0N + r1 + r1N + r2 + r2 + r2N + r3))
     return 8;  
}





