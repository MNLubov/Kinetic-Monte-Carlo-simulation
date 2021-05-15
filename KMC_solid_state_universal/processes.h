#ifndef _H_PROCESSES_H_
#define _H_PROCESSES_H_

void DepositionNearest2D( int DepositionSiteNumber, int SurfaceMapNumber, int AtomType );

void Evaporation( int EvaporationSiteNumber );

int NewDiffusionDirection( int SiteNumber );

void SurfaceDiffusion2D( int DiffusionSiteNumber );

void InitialRatesCalculation( int ConfigurationFrequency );

void RatesUpdate( int ConfigurationFrequency );

void SystemEvolution( void );

#endif
