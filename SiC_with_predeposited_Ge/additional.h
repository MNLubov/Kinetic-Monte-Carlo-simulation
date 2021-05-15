#ifndef _H_ADDITIONAL_H_
#define _H_ADDITIONAL_H_


//number of sites per 1cm2
#define sites_on_cm2  0.66E15

//defing physical constants
#define KbT0  0.025
#define T0  300.0
#define Debay  1.0E13

//define deposition parameters
#define deposition 6.0E13



int AdsiteSum( void );

void GroupList( void);

int Rates( void );

//calculation of the bonds before and after the hop yo the nearest position
int BondsNearest ( int process_state );

#endif
