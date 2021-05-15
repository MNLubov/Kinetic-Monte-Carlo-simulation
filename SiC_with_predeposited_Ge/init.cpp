#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "mesh100.h"
#include "global.h"
#include "processes.h"
#include "outdata.h"
#include "additional.h"
#include "init.h"


void InitDone( void )
{
   int i, j, k;
   
   
   // adsite[i] = 0 - means free adatom site, adsite[i] = 1 - means occupied adatom site
   for(i = 1; i < 20000; i++)
      adsite[i] = 0;
 
   
   for(i = 1; i < 20000; i++)
   {
     adcoords[i][0] = 0;
     adcoords[i][1] = 0;
     adcoords[i][2] = 0;
   }
}
