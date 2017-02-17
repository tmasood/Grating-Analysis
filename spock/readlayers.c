#include <stdio.h>
#include "util.h"
#include "layer.h"

int readlayers(char *strbuffer, struct LAYER *layerptr,
	       struct LAYER **headlayerptr)
{
  char *layerflagptr;
  char *matsysflagptr;

  struct LAYER *tmplayerptr;

  layerptr->nextptr = NULL;
  layerptr->prevptr = NULL;

  if (*headlayerptr == NULL)
    {
      *headlayerptr = layerptr;
    } 
  else
    {
      tmplayerptr = *headlayerptr;
      while(tmplayerptr->nextptr != NULL)
	{
	  tmplayerptr = tmplayerptr->nextptr;
	}
      tmplayerptr->nextptr = layerptr;
      layerptr->prevptr = tmplayerptr;
    }

  layerflagptr = (char *)strstr(strbuffer, "MATSYS=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 7, "%d", &layerptr->matsys);
      printf("MATSYS: material system: %d \t", layerptr->matsys);
      matsysflagptr = (char *)strstr(strbuffer, "XPERC=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 6, "%f", &layerptr->xperc);
	  printf("XPERC- X mole fraction: %g \t", layerptr->xperc);	  
	}
      matsysflagptr = (char *)strstr(strbuffer, "YPERC=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 6, "%f", &layerptr->yperc);
	  printf("YPERC- Y mole fraction: %g", layerptr->yperc);	  
	}
    }

  layerflagptr = (char *)strstr(strbuffer, "NREAL=");
  if (layerflagptr != NULL)
    {
      layerptr->matsys = 0;
      sscanf(layerflagptr + 6, "%lf", &layerptr->nreal);
      printf("NREAL: real index: %lf  \t", layerptr->nreal);
    }

  layerflagptr = (char *)strstr(strbuffer, "NLOSS=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 6, "%lf", &layerptr->nloss);
      printf("NLOSS: imag index: %lf  \t", layerptr->nloss);
    }

  layerflagptr = (char *)strstr(strbuffer, "TL=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 3, "%lf", &layerptr->tl);
      printf(" TL: layer thickness: %lf  \n", layerptr->tl);
    }

  return 0;
}
