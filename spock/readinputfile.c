/*[
 * Copyright 2004, 2005  Taha Masood
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
* with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

#include <stdio.h>
#include "util.h"
#include "layer.h"
#include "case.h"
#include "modcon.h"
#include "struct.h"
#include "output.h"
#include "gamout.h"

int extern readlayers(char *, struct LAYER *, struct LAYER **);
int extern readgradedlayer(char *, struct LAYER **, struct STRUCT *);
int extern readcase(char *, struct CASE **);
int extern readoutput(char *, struct OUTPUT *);
int extern readgamout(char *, struct GAMOUT *);
int extern readloopz(char *, struct STRUCT *);
int extern initialize(struct CASE *, struct STRUCT *, struct MODCON *);

int readinputfile(char *filename, struct CASE **caseptr,
		  struct STRUCT **structptr, struct LAYER **headlayerptr,
		  struct MODCON **modconptr, struct OUTPUT **outputptr)
{
  FILE *fpwgi;
  char strbuffer[MAXCHARS];
  char wgifile[FILESIZE];
  char *strptr;
  char *matchcharptr;
  char *modconflagptr;
  char *structflagptr;
  char designname[MAXCHARS - 1];
  char keyword[MAXKEYWORD];

  struct LAYER *layerptr;
  struct GAMOUT *gamoutptr;

  int num;
  
  layerptr = NULL;

  *caseptr = (struct CASE *)malloc(sizeof(struct CASE));
  *structptr = (struct STRUCT *)malloc(sizeof(struct STRUCT));
  *modconptr = (struct MODCON *)malloc(sizeof(struct MODCON));
  *outputptr = (struct OUTPUT *)malloc(sizeof(struct OUTPUT));
  gamoutptr = (struct GAMOUT *)malloc(sizeof(struct GAMOUT));
  initialize(*caseptr, *structptr, *modconptr);
  strcpy(wgifile, filename);
  strcat(wgifile, ".wgi");
  fpwgi = fopen(wgifile,"r");
  if (fpwgi == NULL)
    {
      printf("Cannot open input file to read \n");
      exit(1);
    }
  
  strptr = fgets(strbuffer, MAXCHARS, fpwgi);
  if (strptr == NULL)
    {
      printf("Error reading input file \n");
      exit(1);
    }

  if (strbuffer[0] == '!')
    {
      matchcharptr = (char *)strchr(strbuffer, '!');
      if (matchcharptr == NULL)
	{
	  printf("character match not found \n");
	  printf("The first line and character of the input ");
	  printf(" file must start with '!' \n");
	  printf("to be followed by the design name \n");
	  exit(1);
	}
      strcpy(designname,(matchcharptr + 1));
      printf("Design Name: %s \n",designname);
    }
  else
    {
      printf("The first line of the input file must start with '!' \n");
      printf("to be followed by the design name \n");
      exit(1);
    }
  
  strptr = fgets(strbuffer, MAXCHARS, fpwgi); 
  while (strptr != NULL)
    {
      sscanf(strbuffer, "%s", keyword);
      if (!strcmp(keyword, "CASE"))
	{
	  readcase(strbuffer, caseptr);
	}
      else if (!strcmp(keyword, "MODCON"))
	{
	  modconflagptr = (char *) strstr(strbuffer, "KPOL=");
	  if (modconflagptr != NULL)
	    {
	      sscanf(modconflagptr + 5, "%d", &(*modconptr)->kpol);
	      printf("KPOL: calculate TE/TM mode: %d \n",
		     (*modconptr)->kpol);
	    }

	  modconflagptr = (char *) strstr(strbuffer, "APB1=");
	  if (modconflagptr != NULL)
	    {
	      sscanf(modconflagptr + 5, "%f", &(*modconptr)->apb1);
	      printf("APB1:  branch specification for branch cuts: %g \n",
		     (*modconptr)->apb1);
	    }

	  modconflagptr = (char *) strstr(strbuffer, "APB2=");
	  if (modconflagptr != NULL)
	    {
	      sscanf(modconflagptr + 5, "%f", &(*modconptr)->apb2);
	      printf("APB2: branch specification for branch cuts:  %g \n",
		     (*modconptr)->apb2);
	    }

	  modconflagptr = (char *) strstr(strbuffer, "KBC1=");
	  if (modconflagptr != NULL)
	    {
	      sscanf(modconflagptr + 5, "%d", &(*modconptr)->kbc1);
	      printf("KBC1: field solutions:  %d \n",
		     (*modconptr)->kbc1);
	    }

	  modconflagptr = (char *) strstr(strbuffer, "KBC2=");
	  if (modconflagptr != NULL)
	    {
	      sscanf(modconflagptr + 5, "%d", &(*modconptr)->kbc2);
	      printf("KBC2: field solutions:  %d \n",
		     (*modconptr)->kbc2);
	    }

	}
      else if (!strcmp(keyword, "STRUCT"))
	{
	  structflagptr = (char *) strstr(strbuffer, "WVL=");
	  if (structflagptr != NULL)
	    {
	      sscanf(structflagptr + 4, "%f", &(*structptr)->wvl);
	      printf("WVL: wavelength: %g \n", (*structptr)->wvl);
	    }	  
	}
      else if (!strcmp(keyword, "LAYER"))
	{
	  layerptr = (struct LAYER *)malloc(sizeof(struct LAYER));
	  readlayers(strbuffer, layerptr, headlayerptr); 
	}
      else if (!strcmp(keyword, "GLAYER"))
	{
	  readgradedlayer(strbuffer, headlayerptr, *structptr); 
	}
      else if (!strcmp(keyword, "OUTPUT"))
	{
	  readoutput(strbuffer, *outputptr);
	}
      else if (!strcmp(keyword, "GAMOUT"))
	{
	  readgamout(strbuffer, gamoutptr);
	}
      else if (!strncmp(keyword, "LOOPZ", 5))
	{
	  printf("Hello keyword \n");
	  readloopz(strbuffer, *structptr);
	}
      strcpy(keyword,"");
      strptr = fgets(strbuffer, MAXCHARS, fpwgi); 
    } 
  fclose(fpwgi);
  return 0;
}
