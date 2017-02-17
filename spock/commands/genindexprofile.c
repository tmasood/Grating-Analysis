#include <stdio.h>
#include <string.h>
#include "util.h"
#include "layer.h"

int genindexprofile(char *filename, struct LAYER *headlayerptr)
{
  struct LAYER *tmplayerptr;
  FILE *fplayerindexfile;
  int layernumber;
  double thickness;
  char layerplotfile[FILESIZE];

  strcpy(layerplotfile, filename);
  strcat(layerplotfile, ".dat");

  layernumber = 0;
  thickness = 0.0;

  fplayerindexfile = fopen(layerplotfile,"w");
  if (fplayerindexfile == NULL)
    {
      printf("Cannot open layer file for write \n");
      exit(1);
    }

  tmplayerptr = headlayerptr;
  while (tmplayerptr != NULL)
    {
      if (layernumber == 0)
	{
	  fprintf(fplayerindexfile,"%lf \t %lf \t %lf \n",
		  thickness, tmplayerptr->nreal, tmplayerptr->nloss);
	}
      else
	{
	  fprintf(fplayerindexfile,"%lf \t %lf \t %lf \n",
		  thickness, tmplayerptr->nreal, tmplayerptr->nloss);
	  thickness+= tmplayerptr->tl;
	  fprintf(fplayerindexfile,"%lf \t %lf \t %lf \n",
		  thickness, tmplayerptr->nreal, tmplayerptr->nloss);
	}
      layernumber++;
      tmplayerptr = tmplayerptr->nextptr;
    }

  fclose(fplayerindexfile);
  return 0;
}
