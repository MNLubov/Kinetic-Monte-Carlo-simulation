#ifndef _H_CONFIGURATIONS_H_
#define _H_CONFIGURATIONS_H_

void BoxInitialConfiguration( void );

void DeleteConfigurationElement( int AtomNumber, int ProcessIndex );

void MoveConfigurationElement( int AtomNumber, double Energy, int ProcessIndex, int AddMoveInde );

void BoxConfigurationUpdate( int OldPositionNumber, int NewPositionNumber, int ProcessIndex);

#endif
