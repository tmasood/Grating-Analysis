#include <stdio.h>
#include <string.h>
#include "util.h"
#include "struct.h"
#include "layer.h"

int genplotindexprofile(char *filename, struct LAYER *headlayerptr,
			struct STRUCT *structptr)
{
  struct LAYER *tmplayerptr;
  FILE *fplayerindexfile;
  FILE *fpindexgpfile;
  int layernumber;
  double thickness;
  char layerplotfile[FILESIZE];
  char gpfile[FILESIZE];
  char psfile[FILESIZE];
  char sysgnuplot[2*MAXCHARS];

  layernumber = 0;
  thickness = 0.0;

  strcpy(layerplotfile, "");
  strcat(layerplotfile, filename);
  strcat(layerplotfile, ".dat");

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

  strcpy(gpfile, "");
  strcat(gpfile, filename);
  strcat(gpfile, ".gp");

  fpindexgpfile = fopen(gpfile,"w");
  if (fpindexgpfile == NULL)
    {
      printf("Cannot open layer file for write \n");
      exit(1);
    }

  strcpy(psfile, "");
  strcat(psfile, filename);
  strcat(psfile, ".ps");

  fprintf(fpindexgpfile,"set term postscript landscape color ");
  fprintf(fpindexgpfile,"\"Times-Roman-bold\" \n");
  fprintf(fpindexgpfile,"set time \n");
  fprintf(fpindexgpfile,"set output \"");
  fprintf(fpindexgpfile,"%s\" \n",psfile);
  fprintf(fpindexgpfile,"set xlabel \"X (um)\" \n");
  fprintf(fpindexgpfile,"set ylabel \"Refractive index\" \n");
  fprintf(fpindexgpfile,"set title \"Plot of Refractive index ");
  fprintf(fpindexgpfile,".vs. X (um)\" \n");
  fprintf(fpindexgpfile,"plot \"");
  fprintf(fpindexgpfile,"%s\" using 1:2 title \"Refractive Index ",
	  layerplotfile);
  fprintf(fpindexgpfile,"\" with lines lw 4 \n");
  fclose(fpindexgpfile);

  strcpy(sysgnuplot, "gnuplot ");
  strcat(sysgnuplot, gpfile);
  system(sysgnuplot);

  system("display -rotate 90 test.ps");
  return 0;
}
