#include <stdio.h>

#include "global.h"

double adcoords[20000][3];

//adatom coords
double adatom_coords[20000][3];
/*
//atom with 1 neighbor coords
double atom1_coords[20000][3];
*/
double attom1_coords[20000][3];

//atom with 2 neighbor coords
double atom2_coords[20000][3];

//atom with 3 neighbor coords
double atom3_coords[20000][3];

//current amount of deposited atoms
int atom_counter;

//new x, y, z coordinates of atom in simulation box 
double atom_coords_new[3];

//current x, y, z coordinates of atom in simulation box 
double atom_coords_curr[3];

//old x, y, z coordinates of atom in simulation box 
double atom_coords_old[3];

//nearest neighbors list of adatom
int neighbor[5];

//next nearest neighbors list of adatom
int Nextneigh[5];

//adatom sites
int adsite[20000];

//number of the cell 
int cell_curr;

//current numbor of the adatom
int atom_curr;

int global_counter;


//groups of different atoms: adatoms, 1-bonded atomes, 2-bonded atoms, 3-bonded atoms
//adatom list
int Adatom[20000];
//1-bonded atoms list
int Atom1[20000];
//2-bonded atoms list
int Atom2[20000];
//3-bonded atoms list
int Atom3[20000];
// group of adatoms and atoms having next nearest neighbors
//adatom list
int AdatomN[20000];
//1-bonded atoms list
int Atom1N[20000];
//2-bonded atoms list
int Atom2N[20000];
//3-bonded atoms list
int Atom3N[20000];

//schetchiki razlichnyh tipov atomov
int adatom_k, atom1_k, atom2_k, atom3_k;
//schetchiki atomov s ucheyom 'next nearest neighbors'
int adatom_kN, atom1_kN, atom2_kN, atom3_kN;

//number of nearest neighbors of each atom 0 - no neghbors, 1 - one neighbor and so on
//int Aneigh[20000];
//number of next nearest neighbors of each atom
//int ANneigh[20000];

//number of the nearest neighbors after the hop
//[i][][] -number of the site
//[i][], 0- not used, 1 - left up diag (nearest), 2 - right up diagonal (nearest), 3 - right down diagonal (nearest), 4 - left down diagonal (nearest),
//5 - left horizontal (next nearest), 6 - vertical up (next nearest), 7 - right horizontal (next nearest), 8 - vertical down (next nearest),
//[][i] - number of the bonds to break
int NeighHopdiff[20000][9][1];

//number of atom with same number of neighbors before and after the hop
extern int Nsame[20000];

//direction of diffusion
int directionDiff;

//massiv razlozhenija chisla time steps po stepenjam
int timedegree[20];

// filename time values
char timedata[8];

int atomcell[1000][2];

int atomN[20000][5];

//
FILE *RateData, *F5;

//group list function counter
int grouplist_count;


//sum cheking, in any functions
int sum_check;
int cc;

FILE *ch;

//Global variables
//number of the nearest neighbors of current atom before the hop
int sumG;
//number of the next nearest neighbors of current atom 
int sumNG;
// number of the nearest neighbors of current after before the hop
int sumAG[5]; 
// number of the next nearest neighbors of current atom after the hop
int sumANG[5]; 
//number of the nearest neighbors before the hop
int BhopG[5];
//number of the next nearest neighbors before the hop
int BhopNG[5];
int sumDiffG;

//number of the directions to hop
//0- not used, 1 - left up diag (nearest), 2 - right up diagonal (nearest), 3 - right down diagonal (nearest), 4 - left down diagonal (nearest),
//5 - left horizontal (next nearest), 6 - vertical up (next nearest), 7 - right horizontal (next nearest), 8 - vertical down (next nearest),
//int Hopdir[20000][9];

//allowed directions of hops
int direction[5];

//counter of the allowed diffusion direction
int jDir;

//array of the allowed diffusion directions
int Diffdir[5];

int ngh;

