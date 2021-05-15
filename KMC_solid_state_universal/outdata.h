#ifndef _H_OUTDATA_H_
#define _H_OUTDATA_H_

#include <string>

using namespace std;

void SimulationBoxOutput( int BeginEndIndex );

void NeighborsListOutput( void );

void ConfigurationOutput( string ConfigurationName );

void SystemStateOutput( string SystemOutName );

void RatesOutput( string RatesOutName );

void SurfaceMapOutput( string SurfaceMapOutName );

void AtomConfigurationNumberOutput( string AtomConfOutName );

#endif
