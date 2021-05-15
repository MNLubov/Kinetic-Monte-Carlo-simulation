 #include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>



#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"
#include "impurity.h"

//defing physical constants
#define KbT0  0.025
#define T0  300.0
#define Debay  1.0E13

//define energies
#define Ea 0.7
#define EaDiag  1.0
#define EaN  0.3
#define EaNN  0.07

#define dEa 0.3

//intercation of impurity and atom
#define dEaN -0.3  //dEaN>0 attractive impurity, dEaN <0 repulsive impureity

int AdsiteSum ( void )
{
   int i;
   
   int sum = 0;
   
   
   for(i = 0; i <= xcell * ycell * 2; i++)
     sum = sum + adsite[i];
     
   return sum;
}

//calculation of the bonds before and after the hop to the nearest position without impurity neigbors
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
      //hop of the atom with the same number of neighbors before and after the hop
      if(sumG != 0) 
      {
        for(j = 1; j <= 4; j++)
       {
         if((sumG - sumAG[j]) == 0 && adsite[BhopG[j]] == 0)
           direction[j] = j; 
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
    //printf("atomNumber = %d, atom1_k = %d\n",atomNumber, atom1_k); 
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
        }  
      }   
      //hop of the atom with the same number ofneighbors before and after the hop
      if(sumG != 0) 
      {
        for(j = 1; j <= 4; j++)
        {
         if((sumG - sumAG[j]) == 0 && adsite[BhopNG[j]] == 0)
           direction[j] = j; 
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
  //printf("jDir=%d\n",jDir);
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
  //number of impurity neighbors
  int sumImp;
  //number of the nearest neighbors before the hop
  int Bhop[5];
  //number of the next nearest neighbors before the hop
  int BhopN[5];
  int sumDiff;
  
  FILE *Fimp;
  
  //uniqflaq[1] - same number of bonds or adatom
  //uniqflaq[2] - break 1 bond, nearest neighbor
  //uniqflaq[3] - break 2 bonds, nearest neighbor
  //uniqflaq[4] - break 3 bonds, nearest neighbor
  //uniqflaq[5] - diagonal hops with same number of neighbors or adatom
  //uniqflaq[6] - diagonal hop, break 1 bond
  //uniqflaq[7] - diagonal hop, break 2 bonds
  //uniqflaq[8] - diagonal hop, break 3 bonds
  int uniqflaq[9];

  //impurity flaq
  //impflaq[1] - 1 impurity neighbor
  //impflaq[2] - 2 impurity neighbors
  //impflaq[3] - 3 impurity neighbors
  //impflaq[4] - 4 impurity neighbors
  //uniqflaq[5] - 
  //uniqflaq[6] - 
  //uniqflaq[7] - 
  //uniqflaq[8] - 
  int impflaq[5];
  
  //
  int uniqImpflaq[12]; 
   
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
  //for hops of atoms with impurities neighbors
  atom_imp1k = 0;
  atom_imp2k = 0;
  atom_imp3k = 0;
  atom_imp1k_1 = 0;
  atom_imp1k_2 = 0;
  atom_imp2k_1 = 0;
  atom_imp1kN = 0;
  atom_imp2kN = 0;
  atom_imp3kN = 0;
  atom_imp1k_1N = 0;
  atom_imp1k_2N = 0;
  atom_imp2k_1N = 0;
  
  Fimp = fopen("imp","w+");
  for(i = 1; i <= xcell * ycell * 2; i++)
  {
    for(j = 0; j <= 8; j++)
       uniqflaq[j] = 1;
    for(j = 0; j <= 11; j++)
       uniqImpflaq[j] = 1;   
       
          
    if(adsite[i] == 1 && Impur[i] == 0)
    { 
        NeighborList(i);
        sum = 0;
        sumN = 0;
        sumImp = 0;
        for(j = 1; j <= 4; j++)
        {
          sum = sum + adsite[neighbor[j]];
          sumN = sumN + adsite[Nextneigh[j]];  
          Bhop[j] = neighbor[j];
          BhopN[j] = Nextneigh[j];
        }
            
    
        for(j = 1; j <= 4; j++)
        {
          if(Impur[Bhop[j]] == 1)
            sumImp = sumImp + 1; 
        }  
        //printf("sumImp = %d\n",sumImp);
        
        //determination of characteristics for the atoms with impurity neighbors     
        if(sumImp != 0)
        {
          //hops to the nearest positions of atoms with impurities neighbors
          for(k = 1; k <= 4; k++)
          {                  
              if(sumImp == 1 && sum == 1 && adsite[Bhop[k]] == 0 && uniqImpflaq[0] != 2)
              {
                atom_imp1k++;
                atom_imp1[atom_imp1k] = i;
                uniqImpflaq[0] = 2; 
              }
              if(sumImp == 1 && sum == 2 && adsite[Bhop[k]] == 0 && uniqImpflaq[1] != 2)
              {
                atom_imp1k_1++;
                atom_imp1_1[atom_imp1k_1] = i;
                uniqImpflaq[1] = 2;
              }
              if(sumImp == 1 && sum == 3 && adsite[Bhop[k]] == 0 && uniqImpflaq[2] != 2)
              {
                atom_imp1k_2++;
                atom_imp1_2[atom_imp1k_2] = i;
                uniqImpflaq[2] = 2;
              }          
              if(sumImp == 2 && sum == 2 && adsite[Bhop[k]] == 0 && uniqImpflaq[3] != 2)
              {
                atom_imp2k++;  
                atom_imp2[atom_imp2k] = i; 
                uniqImpflaq[3] = 2; 
              }
              if(sumImp == 2 && sum == 3 && adsite[Bhop[k]] == 0 && uniqImpflaq[4] != 2)
              {
                atom_imp2k_1++;  
                atom_imp2_1[atom_imp2k_1] = i;
                uniqImpflaq[4] = 2;  
              }
              if(sumImp == 3 && adsite[Bhop[k]] == 0 && uniqImpflaq[5] != 2)
              {
                atom_imp3k++;   
                atom_imp3[atom_imp3k] = i; 
                uniqImpflaq[5] = 2;           
              } 
           }//end of determination of hops to the nearest positions of atoms with impurities
           
           //hops to the next nearest positions of atoms with impurities neighbors
           for(k = 1; k <= 4; k++)
           {                  
              if(sumImp == 1 && sum == 1 && adsite[BhopN[k]] == 0 && uniqImpflaq[6] != 2)
              {
                atom_imp1kN++;
                atom_imp1N[atom_imp1kN] = i;
                uniqImpflaq[6] = 2; 
              }
              if(sumImp == 1 && sum == 2 && adsite[BhopN[k]] == 0 && uniqImpflaq[7] != 2)
              {
                atom_imp1k_1N++;
                atom_imp1_1N[atom_imp1k_1N] = i;
                uniqImpflaq[7] = 2;
              }
              if(sumImp == 1 && sum == 3 && adsite[BhopN[k]] == 0 && uniqImpflaq[8] != 2)
              {
                atom_imp1k_2N++;
                atom_imp1_2N[atom_imp1k_2N] = i;
                uniqImpflaq[8] = 2;
              }          
              if(sumImp == 2 && sum == 2 && adsite[BhopN[k]] == 0 && uniqImpflaq[9] != 2)
              {
                atom_imp2kN++;  
                atom_imp2N[atom_imp2kN] = i; 
                uniqImpflaq[9] = 2; 
              }
              if(sumImp == 2 && sum == 3 && adsite[BhopN[k]] == 0 && uniqImpflaq[10] != 2)
              {
                atom_imp2k_1N++;  
                atom_imp2_1N[atom_imp2k_1N] = i;
                uniqImpflaq[10] = 2;  
              }
              if(sumImp == 3 && adsite[BhopN[k]] == 0 && uniqImpflaq[11] != 2)
              {
                atom_imp3kN++;   
                atom_imp3N[atom_imp3kN] = i; 
                uniqImpflaq[11] = 2;           
              } 
           }//end of determination of hops to the next nearest positions of atoms with impurities
              
        }
          
        //determination of characteristics for the atoms without impurity neighbors 
        if(sumImp == 0)
        {    
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
            }
            
            
           
            /* claculating of the total number of the adatoms (or atoms) with the same number of tne bond before and after the hop */
            if(sum == 0 && uniqflaq[1] != 2)
            {
              adatom_k++;
              Adatom[adatom_k] = i;
              uniqflaq[1] = 2;
            }
            if(sum == 0 && sumN != 4 && uniqflaq[5] != 2)
            {
              adatom_kN++;
              AdatomN[adatom_kN] = i;
              uniqflaq[5] = 2;
            }
            
              
            for(k = 1; k <=4; k++)
            {
              sumDiff = sum - sumA[k];       
              if(sumDiff == 0 && adsite[Bhop[k]] == 0 && uniqflaq[1] != 2)
              {
                 adatom_k++;
                 Adatom[adatom_k] = i;
                 uniqflaq[1] = 2;
              }
              if(sumDiff == 1 && adsite[Bhop[k]] == 0 && uniqflaq[2] != 2)
              {
                atom1_k++;
                Atom1[atom1_k] = i;
                uniqflaq[2] = 2;
              }
              if(sumDiff == 2 && adsite[Bhop[k]] == 0 && uniqflaq[3] != 2)
              {
                atom2_k++;
                Atom2[atom2_k] = i;
                uniqflaq[3] = 2;
              }
              if(sumDiff == 3 && adsite[Bhop[k]] == 0 && uniqflaq[4] != 2)
              {
                atom3_k++;
                Atom3[atom3_k] = i;
                uniqflaq[4] = 2;
              } 
              
              if(sumDiff < 0 && adsite[Bhop[k]] == 0)
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
              if(sumDiff == 0 && adsite[BhopN[k]] == 0 && uniqflaq[5] != 2)
              {
                adatom_kN++;
                AdatomN[adatom_kN] = i;
                uniqflaq[5] = 2;
              }
              if(sumDiff == 1 && adsite[BhopN[k]] == 0 && uniqflaq[6] != 2)
              {
                atom1_kN++;
                Atom1N[atom1_kN] = i;
                uniqflaq[6] = 2;
              }
              if(sumDiff == 2 && adsite[BhopN[k]] == 0 && uniqflaq[7] != 2)
              {
                atom2_kN++;
                Atom2N[atom2_kN] = i;
                uniqflaq[7] = 2;
              }
              if(sumDiff == 3 && adsite[BhopN[k]] == 0 && uniqflaq[8] != 2)
              {
                atom3_kN++;
                Atom3N[atom3_kN] = i;
                uniqflaq[8] = 2;
              }
              if(sumDiff < 0 && adsite[BhopN[k]] == 0)
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
        } //end of grouplist count for atom without impurity neighbor
     }//end of atom occupation check     
  }// end of grouplist procedure for the current atom
  fclose(Fimp);
}//end of the function 


//calculation of processes rates
int Rates( void )
{
   int i, j;
   double depos_flux;
   double T, alpha;
   double speed; //in ML/sec
   int probal;
   int adsites;
   
   double total_rate1,                  // total rate of processes for atoms without impurity neighbors 
          total_rate2, total_rate3;     // total rates of processes for atoms with impurity neighbors
  
  
   double rand_max = (double)RAND_MAX;
   
   
     
   T = 350;
   alpha = KbT0 * (T / T0);
   //speed = 0.75;
   total_rate1 = 0;
   total_rate2 = 0;
   total_rate3 = 0;
   total_rate = 0;
   
   //testing sum 
   double sum1t, sum2t, sum3t, sum4t, sum5t, sum6t, sum7t, sum8t, sum9t, sum10t, sum11t, sum12t;
   //current sum
   double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13;

      
   
   //calculation of deposition flux in  atoms per site in second
   depos_flux = deposition / sites_on_cm2;
   //calculation of deposition speed
   speed = depos_flux / xcell * ycell * 2;
   //printf("%f\n",speed);
  // getch();
      
   /* calculation of processes rates  */
   
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
   
   /* hops with for atoms with impurities neighbors */
   //rate of diffusion to the nearest positions for atoms with 1 impurity neighbor
   rImp1 = Debay * exp((-1 * (Ea + dEaN)) / alpha );
   rImp1_1 = Debay * exp((-1 * (Ea + dEaN + EaN)) / alpha );
   rImp1_2 = Debay * exp((-1 * (Ea + dEaN + 2 * EaN)) / alpha );   
   //rate of diffusion to the next nearest positions for atoms with 2 impurity neighbors
   rImp2 = Debay * exp((-1 * (Ea + 2 * dEaN)) / alpha );
   rImp2_1 = Debay * exp((-1 * (Ea + 2 * dEaN + EaN)) / alpha );
   //rate of diffusion to the next nearest positions for atoms with 3 impurity neighbors
   rImp3 = Debay * exp((-1 * (Ea + 3 * dEaN)) / alpha );

   //rate of diffusion to the next nearest positions for atoms with 1 impurity neighbor
   rImp1N = Debay * exp((-1 * (EaDiag + dEaN)) / alpha );
   rImp1_1N = Debay * exp((-1 * (EaDiag + dEaN + EaN)) / alpha );
   rImp1_2N = Debay * exp((-1 * (EaDiag + dEaN + 2 * EaN)) / alpha );   
   //rate of diffusion to the next nearest positions for atoms with 2 impurity neighbors
   rImp2N = Debay * exp((-1 * (EaDiag + 2 * dEaN)) / alpha );
   rImp2_1N = Debay * exp((-1 * (EaDiag + 2 * dEaN + EaN)) / alpha );
   //rate of diffusion to the next nearest positions for atoms with 3 impurity neighbors
   rImp3N = Debay * exp((-1 * (EaDiag + 3 * dEaN)) / alpha );
  
   GroupList();
   
   //calculating of rates taking into account number of the neighbors for atoms without impurity neighbors
   //r = DeposRate; // / total_rate;
   r0 = adatom_k * r0;
   r0N = adatom_kN * r0N;
   r1 = atom1_k * r1;
   r1N = atom1_kN * r1N;
   r2 = atom2_k * r2;
   r2N = atom2_kN * r2N;
   r3 = atom3_k * r3;
   r3N = atom3_k * r3N;

   /* calculating of rates taking into account number of the neighbors for atoms with impurity neighbors */
   //hops to the nearst positions
   rImp1 = atom_imp1k  * rImp1;
   rImp1_1 = atom_imp1k_1 * rImp1_1;
   rImp1_2 = atom_imp1k_2 * rImp1_2;  
   rImp2 = atom_imp2k * rImp2;
   rImp2_1 = atom_imp2k_1 * rImp2_1;
   rImp3 = atom_imp3k * rImp3; 
   //hops to the next nearst positions    
   rImp1N = atom_imp1kN  * rImp1N;
   rImp1_1N = atom_imp1k_1N * rImp1_1N;
   rImp1_2N = atom_imp1k_2N * rImp1_2N;  
   rImp2N = atom_imp2kN * rImp2N;
   rImp2_1N = atom_imp2k_1N * rImp2_1N;
   rImp3N = atom_imp3kN * rImp3N; 
      
   total_rate1 = r0 + r0N + r1 + r1N + r2 + r2N + r3;
   total_rate2 = rImp1+ rImp1_1 + rImp1_2 + rImp2 + rImp2_1 + rImp3;
   total_rate3 = rImp1N+ rImp1_1N + rImp1_2N + rImp2N + rImp2_1N + rImp3N;
   total_rate = r + total_rate1 + total_rate2 + total_rate3;
   
   /* calculating of the normalized rate */
   r = r / total_rate;
   //printf("r=%f\n",r);
   //getch();
   r0 = r0 / total_rate;
   r0N =  r0N / total_rate;
   r1 =  r1 / total_rate;
   r1N =  r1N / total_rate;
   r2 =  r2 / total_rate;
   r2N =  r2N / total_rate;
   r3 =  r3 / total_rate;
   r3N = 0.0;
   
   rImp1 = rImp1 / total_rate;
   rImp1_1 = rImp1_1 / total_rate;
   rImp1_2 = rImp1_2 / total_rate;  
   rImp2 = rImp2 / total_rate;  
   rImp2_1 = rImp2_1 / total_rate;  
   rImp3 = rImp3 / total_rate;  
    
   rImp1N = rImp1N / total_rate;  
   rImp1_1N = rImp1_1N / total_rate;  
   rImp1_2N = rImp1_2N / total_rate; 
   rImp2N = rImp2N / total_rate;  
   rImp2_1N = rImp2_1N / total_rate;  
   rImp3N = rImp3N / total_rate;  
   
   sum1t = r;
   sum2t = sum1t + r0;
   sum3t = sum2t + r0N;
   sum4t = sum3t + r1;
   sum5t = sum4t + r1N;
   sum6t = sum5t + r2;
   sum7t = sum6t + r2N;
   sum8t = sum7t + r3;
   
   sum1 = r + r0 + r0N + r1 + r1N + r2 + r2N + r3;
   //rates for hops to the nearests position for atoms with impuritties
   sum2 = sum1 + rImp1;
   sum3 = sum2 + rImp1_1;
   sum4 = sum3 + rImp1_2;
   sum5 = sum4 + rImp2;   
   sum6 = sum5 + rImp2_1;
   sum7 = sum6 + rImp3;
   //rates for hops to the next nearest positions for atoms with impuritties   
   sum8 = sum7 + rImp1N;
   sum9 = sum8 + rImp1_1N;
   sum10 = sum9 + rImp1_2N;
   sum11 = sum10 + rImp2N;
   sum12 = sum11 + rImp2_1N;
   sum13 = sum12 + rImp3N;
   
   
   //printf("r=%1.10f\n",r);
   //printf("r0=%1.15f, adatom_k = %d\n",r0, adatom_k);
   //printf("rImp1=%1.10f, atom_imp1k = %d\n",rImp1, atom_imp1k);
   
   
   // generating of the random number in the interval from [0,1]
   probal = rand() % RAND_MAX;
   prob = probal / rand_max;
   //printf("prob=%f\n",prob);
   //getch();
   
   //fprintf(Frates,"%d  %f\n",count, prob);
   //prob = (probal*probal) / (rand_max * rand_max);
   //printf("%f \n",prob);
   cc = cc + 1;

   if( prob >= 0 && prob <= sum1t)
   {
     sum_check = sum_check + 1;
      //printf("prob=%f  r=%f\n",prob,r);
      //getch();
     adsites = AdsiteSum();
     return 1;
   }
   
   //diffusion of atom with thesame number of the neghbors without impurity neighbors
   if( prob > sum1t && prob <= sum2t)
   {
     //printf("prob=%f  r0=%f\n",prob,r0);
     return 2;
   }
   if( prob > sum2t && prob <= sum3t)
     return 3;
   if( prob > sum3t && prob <=  sum4t)
     return 4;
   if( prob > sum4t && prob <= sum5t)
     return 5;
   if( prob > sum5t && prob <= sum6t)
     return 6;  
   if( prob > sum6t && prob <= sum7t)
     return 7;
   if( prob > sum7t && prob <= sum8t)
     return 8;  
   //hops of atoms with impurity neighbors  
   if( prob > sum1 && prob <= sum2)
     return 9;
   if( prob > sum2 && prob <= sum3)     
     return 10;
   if( prob > sum3 && prob <= sum4)
     return 11;
   if( prob > sum4 && prob <= sum5)     
     return 12;
   if( prob > sum5 && prob <= sum6)
     return 13;
   if( prob > sum6 && prob <= sum7)     
     return 14;
   if( prob > sum7 && prob <= sum8)
     return 15;
   if( prob > sum8 && prob <= sum9)     
     return 16;
   if( prob > sum9 && prob <= sum10)
     return 17;
   if( prob > sum10 && prob <= sum11)     
     return 18;
   if( prob > sum11 && prob <= sum12)
     return 19;
   if( prob > sum12 && prob <= sum13)     
     return 20;  
   if( prob > sum13)    
     printf("Bounds exceed!\n");
}





