#ifndef _H_OUTPUT_H_
#define _H_OUTPUT_H_

//output state  file of the simulation box
extern FILE *output;


void out( void );

void out1( double x, double y, double z );

void DegreeR( double check, double degree );

void OutStep( void );

//vyvod ASCII symvolov ot 1 do 100
void AsciiOut ( void );

#endif
