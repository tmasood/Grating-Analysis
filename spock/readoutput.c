#include <stdio.h>
#include "util.h"
#include "output.h"

int readoutput(char *strbuffer, struct OUTPUT *outputptr)
{
  char *outputflagptr;

  outputflagptr = (char *)strstr(strbuffer, "PHMO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 5, "%d", &outputptr->phmo);
      printf("PHMO: phase integral (for root search PHM < 1): %d \n",
	     outputptr->phmo);
    }

  outputflagptr = (char *)strstr(strbuffer, "GAMMAO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->gammao);
      printf("GAMMAO: confinement factor: %d \n", outputptr->gammao);
    }

  outputflagptr = (char *)strstr(strbuffer, "WZRO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 5, "%d", &outputptr->wzro);
      printf("WZRO: eigenvalue (real root): %d \n", outputptr->wzro);
    }

  outputflagptr = (char *)strstr(strbuffer, "WZIO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 5, "%d", &outputptr->wzio);
      printf("WZIO: eigenvalue (imaginary root): %d \n", outputptr->wzio);
    }

  outputflagptr = (char *)strstr(strbuffer, "QZRO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 5, "%d", &outputptr->qzro);
      printf("WZIO: eigenvalue squared (real root squared): %d \n",
	     outputptr->qzro);
    }

  outputflagptr = (char *)strstr(strbuffer, "QZIO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 5, "%d", &outputptr->qzio);
      printf("WZIO: eigenvalue squared (imaginary root squared): %d \n",
	     outputptr->qzio);
    }

  outputflagptr = (char *)strstr(strbuffer, "FWHPNO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->fwhpno);
      printf("FWHPNO: near field: %d \n", outputptr->fwhpno);
    }

  outputflagptr = (char *)strstr(strbuffer, "FWHPFO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->fwhpfo);
      printf("FWHPFO: far field: %d \n", outputptr->fwhpfo);
    }

  outputflagptr = (char *)strstr(strbuffer, "KMO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 4, "%d", &outputptr->kmo);
      printf("KMO: quality of root guess: %d \n", outputptr->kmo);
    }

  outputflagptr = (char *)strstr(strbuffer, "ITO=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 4, "%d", &outputptr->ito);
      printf("ITO: number of iterations before reaching root: %d \n",
	     outputptr->ito);
    }

  outputflagptr = (char *)strstr(strbuffer, "MODOUT=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->modout);
      printf("MODOUT: running account of calculations: %d \n",
	     outputptr->modout);
    }

  outputflagptr = (char *)strstr(strbuffer, "LYROUT=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->lyrout);
      printf("LYROUT: layer index info for each layer: %d \n",
	     outputptr->lyrout);
    }

  outputflagptr = (char *)strstr(strbuffer, "SPLTFL=");
  if (outputflagptr != NULL)
    {
      sscanf(outputflagptr + 7, "%d", &outputptr->spltfl);
      printf("SPLTFL: split the output file- one for each loop: %d \n",
	     outputptr->spltfl);
    }

  return 0;
}
