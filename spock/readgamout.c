#include <stdio.h>
#include "util.h"
#include "gamout.h"

int readgamout(char *strbuffer, struct GAMOUT *gamoutptr)
{
  char *gamoutflagptr;

  gamoutflagptr = (char *)strstr(strbuffer, "LAYGAM=");
  if (gamoutflagptr != NULL)
    {
      sscanf(gamoutflagptr + 7, "%d", &gamoutptr->laygam);
      printf("LAYGAM: print confinement factor of layer: %d \n",
	     gamoutptr->laygam);
    }

  gamoutflagptr = (char *)strstr(strbuffer, "COMPGAM=");
  if (gamoutflagptr != NULL)
    {
      sscanf(gamoutflagptr + 8, "%d", &gamoutptr->compgam);
      printf("COMPGAM: calculate complex gamma: %d \n", gamoutptr->compgam);
    }

  gamoutflagptr = (char *)strstr(strbuffer, "GAMALL=");
  if (gamoutflagptr != NULL)
    {
      sscanf(gamoutflagptr + 7, "%d", &gamoutptr->gamall);
      printf("GAMALL: print gamma for all layers: %d \n", gamoutptr->gamall);
    }

  return 0;
}
