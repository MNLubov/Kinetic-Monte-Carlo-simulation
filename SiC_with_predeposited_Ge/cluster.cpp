 #include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "additional.h"
#include "impurity.h"
#include "cluster.h"

int cluster[50000];
int cluster_atoms[50000];
int total_cluster_atoms;
int noncluster_atoms[50000];
int total_noncluster_atoms;
//coordinates of the atoms in the system
double sys_adcoords[90000][3];
double atoms_adcoords[90000][3];

//cluster matrix
int cluster_direct[5000][5000];
int cluster_sum[5000][5000];

int cluster_column[20010];

int total_atom_number;

// coincedence flaq of the current column
int coin_flaq;

int row_index;
int row_coin_flaq;

FILE *F_cluster;


//neighbors of the current atom
int curr_neigh[5];

int curr_size [10000];
int cluster_function[10000];
//int atom_cluster[90000];

void DataRead( void )
{
  int read_flaq;
  FILE *F, *F2;
  int i,j,k,l;
  int Nb = 10000;
  std::fstream file;
  
  std::fstream file1, file2;
  
  l = 1;
  for(i = 1; i <= xcell * ycell * 2; i++) 
  {
    cluster[i] = 0;
    adsite[i] = 0;
  }
  read_flaq = 1;  
  if(read_flaq == 1)
  {
    file1.open("final",std::ios::in);
    file.open("final_copy",std::ios::out);
    
    VacantPosCoords(xcell,ycell,zcell);
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
            cluster[l] = k;
            atoms_adcoords[k][0] = sys_adcoords[j][0];
            atoms_adcoords[k][1] = sys_adcoords[j][1];
            atoms_adcoords[k][2] = sys_adcoords[j][2];
            file << "1" << "   " << atoms_adcoords[k][0] << "   "<< atoms_adcoords[k][1]<<"   "<< atoms_adcoords[k][2] << "\n";
            l++;
          } 
        }  
      }  
      j++;
    }   

    
    total_atom_number = l - 1;
   
    file.close();
  }

    
}

void MatrixMade( void )
{
  int i,j,k,l;
  std::fstream file;
  
  k = 1;
  l = 1;
  total_cluster_atoms = 0;
  total_noncluster_atoms = 0;
  
  for(i = 0; i <= total_cluster_atoms; i++)
    for(j = 0; j <= total_cluster_atoms; j++)
      cluster_direct[i][j] = 0;
  
  for(i = 1; i <= total_atom_number; i++)
  {
    NeighborList(cluster[i]);
    //calculating total number of atoms in cluster (cluster is group of atom with more than 2)
    if(adsite[neighbor[1]] == 1 || adsite[neighbor[2]] == 1 || adsite[neighbor[3]] == 1 || adsite[neighbor[4]] == 1)
    {
      cluster_atoms[k] = cluster[i];
      total_cluster_atoms = total_cluster_atoms + 1;
      k++;
    }  
    //calculating total number of adatoms
    if(adsite[neighbor[1]] == 0 && adsite[neighbor[2]] == 0 && adsite[neighbor[3]] == 0 && adsite[neighbor[4]] == 0)
    {
      noncluster_atoms[l] = cluster[i];
      total_noncluster_atoms = total_noncluster_atoms + 1;
      l++;
    } 
  }  
  
  //setting values of direct matrix
  i = 1;  
  while(i <= total_cluster_atoms)
  {
     j = 2;
     cluster_direct[i][1] = cluster_atoms[i];
     NeighborList(cluster_atoms[i]);
     for(k = 1; k <= 4; k++)
     {
        //printf("%d  %d\n",k, adsite[neighbor[k]]);
        //getch();
        if(adsite[neighbor[k]] == 1)
        {
          cluster_direct[i][j] = neighbor[k];
          //printf("%d\n",neighbor[k]);
          j++;
        }  
        cluster_direct[i][total_cluster_atoms] = j - 1;
     }   
     i++;
  }
 
  //getch();
  //output surface atoms data
  F_cluster = fopen("cluster_data","w+");
  fprintf(F_cluster,"total atoms in clusters = %d\n",total_cluster_atoms);  
  fprintf(F_cluster,"total number of adatoms = %d\n",total_noncluster_atoms);
  fprintf(F_cluster,"total number of atoms = %d\n",total_atom_number);
  file.open("direct_matrix",std::ios::out);

  for(i = 1; i <= total_cluster_atoms; i++)
  {
    for(j = 1; j <= total_cluster_atoms; j++)
    {
      file << cluster_direct[i][j] << "  ";
    }
    file <<"\n";
    //printf("output row %d\n",i);
    //getch();
  }
  
  //fclose(F_cluster);
  file.close();
}


int RowCompare( int a, int b )
{
  int i,j,k;
  int row_data[5000];
  int row_ind_a, row_ind_b;
  int items_coin;
  
  coin_flaq = 0;
  items_coin = 1;
  
  //opredelenie sovpadenija rjadov
  row_coin_flaq = 0;
  row_ind_a = cluster_direct[a][total_cluster_atoms]; 
  row_ind_b = cluster_direct[b][total_cluster_atoms]; 
  i = 1;
  j = 2;
  while(coin_flaq != 1)
  {
    if(cluster_direct[a][i] == cluster_direct[b][j])   
    {
      coin_flaq = 1;
      row_coin_flaq = 1;
    }  
    
    if(cluster_direct[a][i] != cluster_direct[b][j])
      j++;
    
    if(j == (row_ind_b+1))
    {
      i++;
      j = 2;
    }
    
    if(i == (row_ind_a + 1) && j == 2)
    {
      coin_flaq = 1;
      row_coin_flaq = 0;
      //printf("row_coin_flaq=%d\n",row_coin_flaq);
      //getch();
    }       
  }
  

  //setting new values for row "a"
  if(row_coin_flaq == 1)
  {
     k = 1;     
     cluster_direct[b][total_cluster_atoms] = -1;
     for(j = 1; j <= row_ind_b; j++)
     {
       items_coin = 0;
       //look for same elements in rows "a" and "b"
       for(i = 1; i <= row_ind_a; i++)
       {  
        if(cluster_direct[b][j] == cluster_direct[a][i])
           items_coin = items_coin + 1;
       }
       
       //updating "a" row
       if(items_coin == 0)
       {
           cluster_direct[a][row_ind_a + k] = cluster_direct[b][j];
           k++;
       }  

     }    
     cluster_direct[a][total_cluster_atoms] = cluster_direct[a][total_cluster_atoms] + k - 1;
  }
  
  //printf("row_coin_flaq=%d\n",row_coin_flaq);
  return row_coin_flaq; 
 

}


void ClusterData( void )
{
  int i,j,k,l,m,n,p,t;
  int j_index;
  int sum = 0;
  int rel_index;
  int global_coin_flaq = 0;
  
  int coin_index;
  
  std::fstream file;
  int temp;
  int rows[5000];

  
  std::fstream file1, file2;
  /*
  for(i = 1; i <= total_cluster_atoms; i++)
    rows[i] = 0;  //row[i] = 0 means stroka ne sravnivalas'
  */
  //temp = RowCompare(1,2);
  //printf("%d",temp);
  //getch();
  
  MatrixMade();
  //calculating of the direct matrix
  i = 1;  
  j = 2;
  p = 1;
  row_index = 0;
  while(global_coin_flaq != 1)
  {
    sum = -1;
    p = i+1;
    rel_index = 0;

    while(sum < 0)
    {
      //printf("i=%d\n",i);
      if(cluster_direct[p][total_cluster_atoms] < 0)   
      {
        p++;
        //printf("p=%d\n",p);
      }
      if(cluster_direct[p][total_cluster_atoms] > 0)
      {
        sum = 0;    
        rel_index = 0;
        //printf("p=%d\n",p);
      }
      if(p == total_cluster_atoms + 1)
      {
        rel_index = 1;
        sum = 1;
        //printf("p=%d\n",p);
      }  
          
    }
    
    if(rel_index == 1)
      global_coin_flaq = 1;
    
    //printf("%d  %d\n",p,i);  
    if(rel_index == 0 )
    {
    if(cluster_direct[i][total_cluster_atoms] == -1 || cluster_direct[i][total_cluster_atoms-1] == -5)
    {
      i++;
      p = i;
      //printf("i=%d\n",i);
    } 
    
    if(cluster_direct[i][total_cluster_atoms] != -1 && i <= total_cluster_atoms )
    {
      //printf("i=%d  cluster_direct[i][total_cluster_atoms]=%d \n",i,cluster_direct[i][total_cluster_atoms]);
      //nahodim blizhaishuju sovpadajuschuju stroku
      j_index = i+1;
      //printf("j=%d\n",j);
      while(coin_index != 1)
      {
        //printf("coin_index=%d\n",coin_index);
        //printf("i=%d  j_index=%d  cluster_direct=%d\n",i,j_index,cluster_direct[j_index][total_cluster_atoms]);
        if(cluster_direct[j_index][total_cluster_atoms] == -1 && j_index <= total_cluster_atoms) 
        {
          j++;
          //printf("third line j=%d\n",j);
          //getch();
        }
        
        if(cluster_direct[j_index][total_cluster_atoms] != -1 && j_index <= total_cluster_atoms ) 
        {
          row_index = RowCompare(i,j_index);
          //printf("row_index=%d  %d  %d\n",row_index, i, j);
          
        }  
      
        if(row_index == 1)
        {
          coin_index = 1;  
          //printf("coin_index=%d\n",coin_index);
        }
        if(row_index == 0)
        {
          j_index = j_index + 1;  
          //printf("");
        }
        if(j_index == total_cluster_atoms + 1)
        {
          //printf("j=%d cluster_direct=%d  total=%d\n",j,cluster_direct[i][total_cluster_atoms],total_cluster_atoms);
          coin_index = 1; 
          cluster_direct[i][total_cluster_atoms-1] = -5;
          //printf("j=%d coin_index=%d\n",j,coin_index);
        }
        //printf("j=%d coin_index=%d\n",j,coin_index);  
      }
    }
    //printf("i=%d\n",i);
    //printf("p=%d  i =%d\n",p,i);
    coin_index = 0;
    row_index = 0;
    
    if(i == total_cluster_atoms+1)
      global_coin_flaq = 1; 
    }  
  
  }//end of total cluster matrix setting 
  
  //printf("calc...");
  //getch();
  /*
  for(i = 1; i <= ; i++)
  {
    file2<< i << "  "<<   curr_size[i] <<"\n";
    file1<< i << "  "<<   cluster_function[i] <<"\n";
  }  
  
  file1.close();
  file2.close();
  */
  file.open("direct_matrix_rewrite",std::ios::out);
  for(i = 1; i <= total_cluster_atoms; i++)
  {
    for(j = 1; j <= total_cluster_atoms; j++)
    {
      file << cluster_direct[i][j] << "  ";
    }
    file <<"\n";
    printf("output row %d\n",i);
  }
  
  
  file.close();
  
}

void ClusterDistribution( void )
{
  int i,j,k;
  int cluster_size;
  int cluster_number;
  double distribution_function[500];
  double mean_cluster_size;
  double sum1, sum2;
  int max_cluster_size;
  
  FILE *F;
  
  for(i = 0; i < 500; i++)
    distribution_function[i] = 0.0;

  
  cluster_number = 0;
  mean_cluster_size = 0.0;
  sum1= 0;
  sum2 = 0;
  max_cluster_size = 0;
  for(i = 1; i <= total_cluster_atoms; i++)
  {
    if(cluster_direct[i][total_cluster_atoms] > 0)
    {
      
      cluster_size = cluster_direct[i][total_cluster_atoms];
      //printf("%d\n",cluster_size);
      //getch();     
      cluster_number = cluster_number + 1;
      distribution_function[cluster_size] = distribution_function[cluster_size] + 1.0;
      if(cluster_size > max_cluster_size)
        max_cluster_size = cluster_size;
    }      
  }  
  
  //printf("%d\n",max_cluster_size);
  //getch();
  for(i = 1; i <= max_cluster_size; i++)
  { 
    sum1 = sum1 + distribution_function[i] * i;
    sum2 = sum2 + distribution_function[i];
    printf("%d  %d\n",i, distribution_function[i]);
  }  
  
  //printf("%d  %d",sum1, sum2);
  //getch();
  
  mean_cluster_size = sum1 / sum2;
  
  F = fopen("distribution","w+");
  for(i = 0; i <= max_cluster_size; i++)
    fprintf(F,"%d  %f\n",i,distribution_function[i]);
  
  fprintf(F_cluster,"mean clusetr size  %f\n",mean_cluster_size);
  fprintf(F_cluster,"number of clusters  %d\n",cluster_number);
  fprintf(F_cluster,"maximum cluster size  %d\n",max_cluster_size);
  fclose(F_cluster);
  fclose(F);
}
