#include <stdio.h>
#include <string.h>
#include "util.h"
#include "struct.h"
#include "modcon.h"
#include "layer.h"
#include "output.h"
#include "case.h"

extern int readinputfile(char *, struct CASE **, struct STRUCT **,
			 struct LAYER **, struct MODCON **,
			 struct OUTPUT **);
extern int commandoptions(char *, struct LAYER *, struct STRUCT *);
extern int guideparams(struct CASE *, struct STRUCT *, struct LAYER *,
		       struct MODCON *, struct OUTPUT *, char *);
extern int numzlp(double, int, double, int *, int *, struct CASE *);
extern int zincr(int, int, struct CASE *, struct STRUCT *);
extern int search(FILE *, struct CASE *, struct STRUCT *, struct LAYER *,
		  struct MODCON *, struct OUTPUT *);
extern int fields(struct CASE *, struct STRUCT *, struct LAYER *,
		  struct MODCON *, struct OUTPUT *, char *);

struct CHARMATRX Smatrix, Tmatrix; /* transfer matrix */

int main()
{
  /* **************************************************************
     --MAIN PROGRAM IS PRIMARILY A CALLING PROGRAM
     ************************************************************ */

  char* PROGRAM="WAVEGUIDE";
  char* VERSION="version 0.04";
  char* REVDATE="April 4, 2005";
  char filename[FILESIZE], curfile[FILESIZE];
  char command;
  int j, n, m, nerr;

  struct LAYER *headlayerptr;
  struct CASE *caseptr;
  struct MODCON *modconptr;
  struct STRUCT *structptr;
  struct OUTPUT *outputptr;

  FILE *fpdb;
  char dbfile[FILESIZE];

  caseptr = NULL;
  structptr = NULL;
  modconptr = NULL;
  headlayerptr = NULL;

  printf("                                      SS           \n");
  printf("                                SSSSSSSSSS         \n");
  printf("                              SSSSSSSSSSSSS        \n");
  printf("                            SSSSSSSSS  SSSSS       \n");
  printf("                         SSSSSSSSSSS     SSS       \n");
  printf("S   SSSS    SSSSSSSSSSSSSSSSSSSSSSS                \n");
  printf("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS                 \n");
  printf(" SSS      SSSSSSSSS SMU SSSSSSSSSS                 \n");
  printf("           SSSSSSSSSSSSSSSSSSSSSS                  \n");
  printf("           SSSSS    SSSSS SSS  SSS                 \n");
  printf("            SSSSS        SSS    SSS                \n");
  printf("              SSSSS    SSSSSSSSSS                  \n");
  printf("                  SSS  SS                          \n");
  printf("                     SS                            \n");
  printf("                     SS                            \n");
  printf("\n");
  printf("             AND PHOTODIGM INC.                    \n");
  printf("\n");
  printf("%s \n", PROGRAM);
  printf("%s \n", VERSION);
  printf("last updated: %s \n",REVDATE);
  printf("\n");
  printf("Copyright 2003, 2004 Photodigm Inc. and Southern Methodist");
  printf(" University, \n                      Dallas, Texas \n");

  printf("\n \n \n ");
  printf("\244 Enter input filename \244 \n \t OR \n");
  printf("\244 Enter <q> to quit OR <c> for other commands \244 \n");
  printf("> ");
  scanf("%s",filename);

  while (!((!strcmp(filename,"Q")) || (!strcmp(filename,"q"))))
    {
      if ((!strcmp(filename,"H")) || (!strcmp(filename,"h")))
	{
	  printf("Help is available in the User Manual \n");
	  printf("WAVEGUIDE program exiting \n");
	  exit(1);
	}
      else if ((!strcmp(filename,"C")) || (!strcmp(filename,"c")))
	{
	  printf("\n");
	  printf("C - Commands \n");
	  printf("E - Edit \n");
	  printf("H - Help (manual) \n");
	  printf("N - New \n");
	  printf("Q - Quit \n");
	  printf("V - View \n");
	  printf("X - eXit \n");
	  printf("\n");
	  printf("> ");
	  scanf("%c",&command);
	  scanf("%c",&command);

	  putchar(command);
	  switch(command)
	    {
	    case 'C':
	      commandoptions(curfile, headlayerptr, structptr);
	      break;
	    case 'c':
	      commandoptions(curfile, headlayerptr, structptr);
	      break;
	    default:
	      break;
	    }
	}
      else
	{
	  strcpy(curfile, filename);
	  readinputfile(curfile, &caseptr, &structptr, &headlayerptr,
			&modconptr, &outputptr);
	  calcindex(headlayerptr, structptr);

	  /* count the number of z loops */
	  n = 0;
	  m = 1;

	  strcpy(dbfile, filename);
	  strcat(dbfile, ".db");

	  fpdb = fopen(dbfile,"w");

	  if (fpdb == NULL)
	    {
	      printf("Cannot open db file to write \n");
	      exit(1);
	    }
	  fprintf(fpdb,"   QZMR \t QZMI \t\t phm \t\t neff \t\t\t\t km \t it \n");
	  fprintf(fpdb,"   ---- \t ---- \t\t --- \t\t ---- \t\t\t\t -- \t -- \n");

	  if (structptr->loopzv)
	    {
	      numzlp(structptr->zfinv[m][1], structptr->ilz[m][1],
		     structptr->zinc[m][1], &n, &nerr, caseptr);
	    }
	  for (j = 0; j <= n; j++)
	    {
	      zincr(j, m, caseptr, structptr);
	      guideparams(caseptr, structptr, headlayerptr, modconptr,
			  outputptr, curfile);
	      search(fpdb, caseptr, structptr, headlayerptr, modconptr, outputptr);
	      fields(caseptr, structptr, headlayerptr, modconptr,
		     outputptr, curfile);
	    }
	}

      printf("\244 Enter <q> to quit OR <c> for other commands \244 \n");
      printf("> ");
      scanf("%s",filename);
    }

  printf("WAVEGUIDE program exiting \n");
  close(fpdb);
  free(caseptr);
  free(structptr); 
  free(headlayerptr);

  return 0;
}

