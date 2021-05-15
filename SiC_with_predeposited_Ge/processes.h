#ifndef _H_PROCESSES_H_
#define _H_PROCESSES_H_



// index of occuring process (deposition, diffusion of adatom, diffusion of atom with 1 neighbor e.t.c)
extern int process;

//rates of the processes
//deposition rate
extern double DeposRate;
//adatom diffusuin rate
extern double AdatomDiff;
//atom with one neighbor rate
extern double Atom1Diff;
//atom with 2 neighbors rate
extern double Atom2Diff;
//atom with 3 neighbors rate
extern double Atom3Diff;
//total rate
extern double total_rate;

//
//extern double Atom1Diff1


//diffusion data file pointer
extern FILE *Diff;

int Deposition( int sort );

int Diffusion(  int process_state, int diff_sort );

void NeighborList( int atom_curr );

#endif
