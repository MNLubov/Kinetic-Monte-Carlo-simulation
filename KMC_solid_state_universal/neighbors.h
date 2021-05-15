#ifndef _H_NEIGHBORS_H_
#define _H_NEIGHBORS_H_

void NeighbourListZB2D ( void );

int NearestNeighborCalc( int SiteNumber );

int NearestNeighborBeforeProcess( int SiteNumber );

void NeighborListCalc ( int DepositedStructureIndex );

int NextNearestNeighborBeforeProcess( int SiteNumber );

#endif
