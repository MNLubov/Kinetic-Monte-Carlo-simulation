#ifndef _H_MESH100_H_
#define _H_MESH100_H_

//total amount of atoms in substrate box
extern int totalatoms;

//total amount of vacant positions
extern int totalvacpos;

//total amount of deposited atoms
extern int total;

//substrate box cell number
extern int xcell;
extern int ycell;
extern int zcell;

//first indicator atom number, second - coordinate 0-x-coord, 1-y-coord
//coords of the substrate box
extern double toplayer[80000][3];

//coords of elementar box
extern double basis [9][3];

//coords of vacant positions
extern double vacantpos[20000][3];

//coords of elementar vacant box 
extern double vacbasis[9][3];

void BasisCoords ( void );

void MakeTopLayerCoords( int xcell, int ycell, int zcell );

void VacBasisCoords( void );

void VacantPosCoords(int xcell, int ycell, int zcell);


#endif
