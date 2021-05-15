#ifndef _H_IMPURITY_H_
#define _H_IMPURITY_H_

//impurity atoms array
extern int Impur[90000];

//number of deposited impurities
extern int ImpNumb;

//number of the atoms with 1 impurity neighbor
extern int atom_imp1k;

//number of the atoms with 2 impurity neighbor
extern int atom_imp2k;

//number of the atoms with 3 impurity neighbor
extern int atom_imp3k;

//number of the atoms with 1 impurity neighbor and 1 C neighbor
extern int atom_imp1k_1;

//number of the atoms with 1 impurity neighbor and C neighbors
extern int atom_imp1k_2;

//number of the atoms with 2 impurity neighbors and 1 C neighbor
extern int atom_imp2k_1;

/* hop of atoms with impuritites on diagonal */
//number of the atoms with 1 impurity neighbor
extern int atom_imp1kN;

//number of the atoms with 2 impurity neighbor
extern int atom_imp2kN;

//number of the atoms with 3 impurity neighbor
extern int atom_imp3kN;

//number of the atoms with 1 impurity neighbor and 1 C neighbor
extern int atom_imp1k_1N;

//number of the atoms with 1 impurity neighbor and C neighbors
extern int atom_imp1k_2N;

//number of the atoms with 2 impurity neighbors and 1 C neighbor
extern int atom_imp2k_1N;

//array of atoms with 1 impurity neioghbor
extern int atom_imp1[20000];

//array of atoms with 1 impurity neioghbor and 1 C neighbor
extern int atom_imp1_1[20000];

//array of atoms with 1 impurity neioghbor and 2 C neighbors
extern int atom_imp1_2[20000];

//array of atoms with 2 impurity neioghbors
extern int atom_imp2[20000];

//array of atoms with 2 impurity neioghbors
extern int atom_imp2_1[20000];

//array of atoms with 3 impurity neioghbors
extern int atom_imp3[20000];

/* arrays of atoms with impuritites that can hop on diagonal */
//array of atoms with 1 impurity neioghbor
extern int atom_imp1N[20000];

//array of atoms with 1 impurity neioghbor and 1 C neighbor
extern int atom_imp1_1N[20000];

//array of atoms with 1 impurity neioghbor and 2 C neighbors
extern int atom_imp1_2N[20000];

//array of atoms with 2 impurity neioghbors
extern int atom_imp2N[20000];

//array of atoms with 2 impurity neioghbors
extern int atom_imp2_1N[20000];

//array of atoms with 3 impurity neioghbors
extern int atom_imp3N[20000];


//number of the next nearest impurity neighbors before the hop
extern int BHopImp[5];



void ImpPreDeposition ( int ImpNumb );

int ImpBonds( int process_state );

#endif
