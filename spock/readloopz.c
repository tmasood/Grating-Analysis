#include <stdio.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int readloopz(char *strbuffer, struct STRUCT *structptr)
{
  int tmploop; /* loop number */
  int value;
  char *loopzptr;
  char loopvar[MAXSIMLOOP];
  
  loopzptr = (char *)strstr(strbuffer, "LOOPZ");
  tmploop = structptr->loopzv;

  if (loopzptr != NULL)
    {
      sscanf(loopzptr + 5, "%d", &tmploop);   
      /* increment number of simultaneous loops */
      structptr->loopzv = tmploop;
      structptr->lzcnt[tmploop] = structptr->lzcnt[tmploop] + 1;
    }
  
  loopzptr = (char *)strstr(strbuffer, "FINV=");
  if (loopzptr != NULL)
    {
      sscanf(loopzptr + 5, "%lf", 
	     &structptr->zfinv[tmploop][structptr->lzcnt[tmploop]]);
    }

  loopzptr = (char *)strstr(strbuffer, "ZINC=");
  if (loopzptr != NULL)
    {
      sscanf(loopzptr + 5, "%lf", 
	     &structptr->zinc[tmploop][structptr->lzcnt[tmploop]]);
    }

  loopzptr = (char *)strstr(strbuffer, "ILZ=");
  if (loopzptr != NULL)
    {
      sscanf(loopzptr + 4, "%s", loopvar);
    }
  if (!strcmp(loopvar, "QZMR"))
    {
      value = 1;
    }
  else if (!strcmp(loopvar, "QZMI"))
    {
      value = 2;
    }
  else
    {
      printf("ZLOOP: unknown parameter \n");
      exit(1);
    }

  structptr->ilz[tmploop][structptr->lzcnt[tmploop]] = value;
  structptr->loopzv = tmploop;

  return (0);
}
