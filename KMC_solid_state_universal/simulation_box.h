#ifndef _H_SIMULATION_BOX_H_
#define _H_SIMULATION_BOX_H_

void SimulationBoxGeneration( void );

void SurfaceMapGenerationZB( int RoughnessIndex );

void SurfaceMapGeneration( int RoughnessIndex );

void SurfaceMapEvolution2D( int SurfaceSiteNumber, int SurfaceMapNumber, int ProcessIndex );

int LastAtomCalculation( int AtomNumber, int ProcessIndex );

double AtomSiteEnergy( int AtomNumber );

#endif
