#include <stdio.h>
#include <conio.h>

#include "parameters.h"
#include "unit_cell.h"
#include "init.h"

// CrystalStructureType = 1 - simple cubic structure
// CrystalStructureType = 2 - bcc (body-centered cubic) structure
// CrystalStructureType = 3 - fcc (face-centered cubic) structure
// CrystalStructureType = 4 - hexagonal structure
// CrystalStructureType = 5 - zb (zinc-blend) structure
// CrystalStructureType = 
// 
void UnitCellGeneration2D (int CrystalStructureType, int UnitCellAtomsTypeNumber, int UnitCellTermination) 
{ 
  int i;
  
  // 2D diamond/zb structure
  if(CrystalStructureType == 5)  
  {
    NearestNeighbornumber = 4;
    NextNearestNeighbornumber = 4;
    
    //diamond structure
    if(UnitCellAtomsTypeNumber == 1)
    {
      //setting unit cell parameters
      UnitCellAtomsNumber = 2; 
      UnitCellSizeX = 1; 
      UnitCellSizeY = 1; 
      UnitCellSizeZ = 1;
      UnitCellLengthX = 0.5;
      UnitCellLengthZ = 0.5;
      
      //setting type of the atoms in the unit cell
      UnitCellAtomsType[1] = 1; 
      UnitCellAtomsType[2] = 1; 
      
      //setting coordinates of atoms in the unit cell
      UnitCellCoordinates[1][1] = 0;                          
      UnitCellCoordinates[1][2] = 0;
      UnitCellCoordinates[1][3] = 0;    
      UnitCellCoordinates[2][1] = 0.5;                          
      UnitCellCoordinates[2][2] = 0.0;
      UnitCellCoordinates[2][3] = 0.5;      
    }
    
    //zb structure
    if(UnitCellAtomsTypeNumber > 1)
    {
      //setting unit cell parameters
      UnitCellAtomsNumber = 2; 
      UnitCellSizeX = 1; 
      UnitCellSizeZ = 1;
      UnitCellLengthX = 0.5;
      UnitCellLengthZ = 0.5;
      
      //setting type of the atoms in the unit cell
      if(UnitCellTermination == 1)
      {
        UnitCellAtomsType[1] = 2;   
        UnitCellAtomsType[2] = 1;   
      }  
      //setting type of the atoms in the unit cell
      if(UnitCellTermination == 2)
      {
        UnitCellAtomsType[1] = 1;   
        UnitCellAtomsType[2] = 2;   
      }  
      //setting coordinates of atoms in the unit cell
      UnitCellCoordinates[1][1] = 0;                          
      UnitCellCoordinates[1][2] = 0;
      UnitCellCoordinates[1][3] = 0;    
      UnitCellCoordinates[2][1] = 0.5;                          
      UnitCellCoordinates[2][2] = 0.0;
      UnitCellCoordinates[2][3] = 0.5;      
    }          
    
  }
}

// 
void UnitCellGeneration3D(int CrystalStructureType) 
{ 
  int i;
  
  if(CrystalStructureType == 1)  // 3D fcc structure
  {
    UnitCellAtomsNumber = 4; 
    UnitCellSizeX = 1; 
    UnitCellSizeY = 1; 
    UnitCellSizeZ = 1;
    for(i = 1; i <= UnitCellAtomsNumber + 1; i++) 
      UnitCellAtomsType[i] = 1;   

    UnitCellCoordinates[1][1] = 0;                          
    UnitCellCoordinates[1][2] = 0;
    UnitCellCoordinates[1][3] = 0;    
    UnitCellCoordinates[2][1] = 0.5;                          
    UnitCellCoordinates[2][2] = 0.5;               
    UnitCellCoordinates[2][3] = 0;    
    UnitCellCoordinates[3][1] = 0.5;                          
    UnitCellCoordinates[3][2] = 0;
    UnitCellCoordinates[3][3] = 0.5;    
    UnitCellCoordinates[4][1] = 0;                          
    UnitCellCoordinates[4][2] = 0.5;
    UnitCellCoordinates[4][3] = 0.5;            
    
  }
}
