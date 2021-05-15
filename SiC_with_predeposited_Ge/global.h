#ifndef _H_GLOBAL_H_
#define _H_GLOBAL_H_


#define lattice_constant 5.43



//occupied site coordinate
extern double adcoords[20000][3];

//adatom coords
extern double adatom_coords[20000][3];
/*
//atom with 1 neighbor coords
double atom1_coords[20000][3];
*/
extern double attom1_coords[20000][3];

//atom with 2 neighbor coords
extern double atom2_coords[20000][3];

//atom with 3 neighbor coords
extern double atom3_coords[20000][3];

//current amount of deposited atoms
extern int atom_counter;

//new x, y, z coordinates of atom in simulation box 
extern double atom_coords_new[3];

//current x, y, z coordinates of atom in simulation box 
extern double atom_coords_curr[3];

//old x, y, z coordinates of atom in simulation box 
extern double atom_coords_old[3];

//nearest neighbors list of adatom
extern int neighbor[5];

//next nearest neighbors list of adatom
extern int Nextneigh[5];

//adatom sites
extern int adsite[20000];

//number of the cell 
extern int cell_curr;

//current numbor of the adatom
extern int atom_curr;

extern int global_counter;


//groups of different atoms having nearest neighbors: adatoms, 1-bonded atomes, 2-bonded atoms, 3-bonded atoms 
//adatom list
extern int Adatom[20000];
//1-bonded atoms list
extern int Atom1[20000];
//2-bonded atoms list
extern int Atom2[20000];
//3-bonded atoms list
extern int Atom3[20000];
// group of adatoms and atoms having next nearest neighbors
//adatom list
extern int AdatomN[20000];
//1-bonded atoms list
extern int Atom1N[20000];
//2-bonded atoms list
extern int Atom2N[20000];
//3-bonded atoms list
extern int Atom3N[20000];

//schetchiki razlichnyh tipov atomov
extern int adatom_k, atom1_k, atom2_k, atom3_k;
//schetchiki atomov s ucheyom 'next nearest neighbors'
extern int adatom_kN, atom1_kN, atom2_kN, atom3_kN;

//number of nearest neighbors of each atom 0 - no neghbors, 1 - one neighbor and so on
//extern int Aneigh[5];
//number of next nearest neighbors of each atom
extern int ANneigh[5];

//number of the nearest neighbors after the hop
//[i][][] -number of the site
//[][i][], 0- not used, 1 - left up diag (nearest), 2 - right up diagonal (nearest), 3 - right down diagonal (nearest), 4 - left down diagonal (nearest),
//5 - left horizontal (next nearest), 6 - vertical up (next nearest), 7 - right horizontal (next nearest), 8 - vertical down (next nearest),
//[][][i] - number of the neighbors
//extern int NeighHopDiff[20000][9][1];

//number of atom with same number of neighbors before and after the hop
//extern int Nsame[20000];

//direction of diffusion
extern int directionDiff;

//massiv razlozhenija chisla time steps po stepenjam
extern int timedegree[20];

// filename time values
extern char timedata[8];


//
extern FILE *RateData, *F5;

//group list function counter
extern int grouplist_count;

//massiv nomerov jacheek, v kotoroj raspolozhen atom
extern int atomcell[1000][2];

//atom-neighbors of the adatom 
extern int atomN[20000][5];


//sum cheking, in any functions
extern int sum_check;
extern int cc;

extern FILE *ch;

//Global variables
//number of the nearest neighbors of current atom before the hop
extern int sumG;
//number of the next nearest neighbors of current atom 
extern int sumNG;
// number of the nearest neighbors of current after before the hop
extern int sumAG[5]; 
// number of the next nearest neighbors of current atom after the hop
extern int sumANG[5]; 
//number of the nearest neighbors before the hop
extern int BhopG[5];
//number of the next nearest neighbors before the hop
extern int BhopNG[5];
extern int sumDiffG;

//number of the directions to hop
//0- not used, 1 - left up diag (nearest), 2 - right up diagonal (nearest), 3 - right down diagonal (nearest), 4 - left down diagonal (nearest),
//5 - left horizontal (next nearest), 6 - vertical up (next nearest), 7 - right horizontal (next nearest), 8 - vertical down (next nearest),
//extern int Hopdir[20000][9];

//directions of hops
extern int direction[5];

//counter of the allowed diffusion direction
extern int jDir;

//array of the allowed diffusion directions
extern int Diffdir[5];

extern int ngh;


#endif
