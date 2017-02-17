#include <stdio.h>
#include <string.h>
#include "util.h"
#include "layer.h"

int createlayerfile(char *filename, struct LAYER *headlayerptr)
{
  struct LAYER *tmplayerptr;
  FILE *fplayerindexfile;
  char layerfile[FILESIZE];

  strcpy(layerfile, "");
  strcat(layerfile, filename);
  strcat(layerfile, ".ly");

  fplayerindexfile = fopen(layerfile,"w");
  if (fplayerindexfile == NULL)
    {
      printf("Cannot open layer file for write \n");
      exit(1);
    }

  tmplayerptr = headlayerptr;
  while (tmplayerptr != NULL)
    {
      fprintf(fplayerindexfile,"%g \t %g \t %g \n",
	      tmplayerptr->nreal, tmplayerptr->nloss, tmplayerptr->tl);

      tmplayerptr = tmplayerptr->nextptr;
    }

  fclose(fplayerindexfile);
  return 0;
}
