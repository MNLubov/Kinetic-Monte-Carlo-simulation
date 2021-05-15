#ifndef _H_GLOBAL_H_
#define _H_GLOBAL_H_


#define lattice_constant 5.43



//occupied site coordinate
extern double adcoords[90000][3];

//adatom coords
extern double adatom_coords[90000][3];

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
extern int adsite[90000];

//number of the cell 
extern int cell_curr;

//current numbor of the adatom
extern int atom_curr;

extern int global_counter;


//groups of different atoms having nearest neighbors: adatoms, 1-bonded atomes, 2-bonded atoms, 3-bonded atoms 
//adatom list
extern int Adatom[50000];
//1-bonded atoms list
extern int Atom1[50000];
//2-bonded atoms list
extern int Atom2[50000];
//3-bonded atoms list
extern int Atom3[50000];
// group of adatoms and atoms having next nearest neighbors
//adatom list
extern int AdatomN[50000];
//1-bonded atoms list
extern int Atom1N[50000];
//2-bonded atoms list
extern int Atom2N[50000];
//3-bonded atoms list
extern int Atom3N[50000];

//schetchiki razlichnyh tipov atomov
extern int adatom_k, atom1_k, atom2_k, atom3_k;
//schetchiki atomov s ucheytom 'next nearest neighbors'
extern int adatom_kN, atom1_kN, atom2_kN, atom3_kN;

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

/* sum cheking, in any functions */
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

//name of the system configuration output file
extern char timedata[8];

//direction of the hopto the next nearest position of atom with impurity neighbor
extern int dirImp;

//direction of the hop to the next nearest position of atom with impurity neighbor
extern int dirImpN;

//new atom position for diffusion of atom with impurity neighbor
extern int newpos;

//probability in the range [0,1]
extern double prob;

//number of the simulation step
extern int count;


//output file for 'int Rates(void)' function
extern FILE *Frates;

// rates of hops without keeping the number of neighbors including adatom hops rate 
extern double r, r0, r1, r2, r3;
// rates of hops with keeping the same number of nearest neighbors
extern double r0N, r1N, r2N, r3N;
// rates of hops to the nearest positions for atoms with impurity neighbors
extern double rImp1, rImp2, rImp3, rImp1_1, rImp1_2, rImp2_1; 
//rates of hops to the next nearest positions for atoms with impurity neighbors
extern double rImp1N, rImp2N, rImp3N, rImp1_1N, rImp1_2N, rImp2_1N; 

#endif
