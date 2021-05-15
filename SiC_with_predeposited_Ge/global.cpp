#include <stdio.h>

#include "global.h"

double adcoords[90000][3];

//adatom coords
double adatom_coords[90000][3];

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
int adsite[90000];

//number of the cell 
int cell_curr;

//current numbor of the adatom
int atom_curr;

int global_counter;

//groups of different atoms: adatoms, 1-bonded atomes, 2-bonded atoms, 3-bonded atoms
//adatom list
int Adatom[50000];
//1-bonded atoms list
int Atom1[50000];
//2-bonded atoms list
int Atom2[50000];
//3-bonded atoms list
int Atom3[50000];
// group of adatoms and atoms having next nearest neighbors
//adatom list
int AdatomN[50000];
//1-bonded atoms list
int Atom1N[50000];
//2-bonded atoms list
int Atom2N[50000];
//3-bonded atoms list
int Atom3N[50000];

//schetchiki razlichnyh tipov atomov
int adatom_k, atom1_k, atom2_k, atom3_k;
//schetchiki atomov s ucheyom 'next nearest neighbors'
int adatom_kN, atom1_kN, atom2_kN, atom3_kN;

//direction of diffusion
int directionDiff;

//massiv razlozhenija chisla time steps po stepenjam
int timedegree[20];

// filename time values
char timedata[8];

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

//name of the system configuration output file
extern char timedata[8];

//direction of the hopto the next nearest position of atom with impurity neighbor
int dirImp;

//direction of the hop to the next nearest position of atom with impurity neighbor
int dirImpN;

//new atom position for diffusion of atom with impurity neighbor
int newpos;

//probability in the range [0,1]
double prob;

//number of the simulation step
int count;

//output file for 'int Rates(void)' function
FILE *Frates;

// rates of hops without keeping the number of neighbors including adatom hops rate 
double r, r0, r1, r2, r3;
// rates of hops with keeping the same number of nearest neighbors
double r0N, r1N, r2N, r3N;
// rates of hops to the nearest positions for atoms with impurity neighbors
double rImp1, rImp2, rImp3, rImp1_1, rImp1_2, rImp2_1; 
//rates of hops to the next nearest positions for atoms with impurity neighbors
double rImp1N, rImp2N, rImp3N, rImp1_1N, rImp1_2N, rImp2_1N; 
