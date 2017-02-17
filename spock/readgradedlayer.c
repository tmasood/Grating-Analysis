#include <stdio.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int readgradedlayer(char *strbuffer, struct LAYER **headlayerptr,
		    struct STRUCT *structptr)
{
  char *layerflagptr;
  char *matsysflagptr;

  struct GLAYER glayer;
  struct LAYER layerstart, layerend;
  struct LAYER *layerptr, *tmplayerptr;

  double indexstep, widthstep;
  int i;

  layerflagptr = (char *)strstr(strbuffer, "MATSYS=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 7, "%d", &glayer.gradedmatsys);
      printf("GRADED MATSYS: material system: %d \t", glayer.gradedmatsys);
      matsysflagptr = (char *)strstr(strbuffer, "XPERCI=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 7, "%f", &glayer.xperci);
	  printf("XPERCI- Initial X mole fraction: %g \t", glayer.xperci);
	}
      matsysflagptr = (char *)strstr(strbuffer, "YPERCI=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 7, "%f", &glayer.yperci);
	  printf("YPERCI- Initial Y mole fraction: %g \t", glayer.yperci);
	}
      matsysflagptr = (char *)strstr(strbuffer, "XPERCF=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 7, "%f", &glayer.xpercf);
	  printf("XPERCF- final X mole fraction: %g \t", glayer.xpercf);  
	}
      matsysflagptr = (char *)strstr(strbuffer, "YPERCF=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 7, "%f", &glayer.ypercf);
	  printf("YPERCF- final Y mole fraction: %g \t", glayer.ypercf);	  
	}
      matsysflagptr = (char *)strstr(strbuffer, "NSLC=");
      if (matsysflagptr != NULL)
	{
	  sscanf(matsysflagptr + 5, "%d", &glayer.nslc);
	  printf("NSLC- number of slices: %d \t", glayer.nslc);	  
	}
    }

  layerflagptr = (char *)strstr(strbuffer, "NREALI=");
  if (layerflagptr != NULL)
    {
      glayer.gradedmatsys = 0;
      sscanf(layerflagptr + 7, "%lf", &glayer.nreali);
      printf("NREALI: real index: %lf  \t", glayer.nreali);
    }

  layerflagptr = (char *)strstr(strbuffer, "NREALF=");
  if (layerflagptr != NULL)
    {
      glayer.gradedmatsys = 0;
      sscanf(layerflagptr + 7, "%lf", &glayer.nrealf);
      printf("NREALF: real index: %lf  \t", glayer.nrealf);
    }

  layerflagptr = (char *)strstr(strbuffer, "NLOSS=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 6, "%lf", &glayer.nloss);
      printf("NLOSS: imag index: %lf  \t", glayer.nloss);
    }
  else
    {
      glayer.nloss = 0.0;
    }

  layerflagptr = (char *)strstr(strbuffer, "TL=");
  if (layerflagptr != NULL)
    {
      sscanf(layerflagptr + 3, "%lf", &glayer.gtl);
      printf("TL: layer thickness: %lf  \n", glayer.gtl);
    }

  layerstart.matsys = glayer.gradedmatsys;
  layerstart.xperc =  glayer.xperci;
  layerstart.yperc =  glayer.yperci;
  layerstart.nloss =  glayer.nloss;
  layerstart.nextptr = NULL;
  calcindex(&layerstart, structptr);

  layerend.matsys = glayer.gradedmatsys;
  layerend.xperc =  glayer.xpercf;
  layerend.yperc =  glayer.ypercf;
  layerend.nloss =  glayer.nloss;
  layerend.nextptr = NULL;
  calcindex(&layerend, structptr);

  indexstep = (layerend.nreal - layerstart.nreal) / glayer.nslc;
  widthstep = glayer.gtl / glayer.nslc;

  for (i=0; i<=glayer.nslc; i++)
    {
      layerptr = (struct LAYER *)malloc(sizeof(struct LAYER));

      layerptr->nextptr = NULL;
      layerptr->matsys = 0;
      layerptr->nreal = layerstart.nreal + i*indexstep;
      layerptr->nloss = glayer.nloss;
      if ((i == 0) || (i == glayer.nslc))
	{
	  layerptr->tl = widthstep/2.0;
	}
      else
	{
	  layerptr->tl = widthstep;
	}

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
    }

  return 0;
}
