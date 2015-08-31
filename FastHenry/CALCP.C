/* this calls the routine to calculate the filament-filament 
   interaction exactly */

#include "induct.h"
#include "../FHWindow.h"

static int num2nd=0, num4th=0, numexact=0;
static int num2ndsav=0, num4thsav=0, numexactsav=0;

double calcp(pchg1, pchg2, pfd)
charge *pchg1, *pchg2;
double *pfd;   /* left over from fastcap */
{
  double mutual(), selfterm();

  if (pfd != NULL)
    viewprintf(stderr, "calcp: I don't know what to do with pfd!=NULL\n");

  if (pchg1->fil->filnumber == pchg2->fil->filnumber) 
    /* self term */
    return selfterm(pchg1->fil);
  else
    /* calculate mutual inductance of the two filaments */
    return mutual(pchg1->fil, pchg2->fil);
} 

/* from the fastcap calcp */
void dumpnums(flag, size)
int flag, size;
{
  double total;

  if(flag == ON) {		/* if first call */
    num2ndsav = num2nd;
    num4thsav = num4th;
    numexactsav = numexact;
  }
  else {
    total = num2ndsav + num4thsav + numexactsav;
#if MULDAT == ON
    viewprintf(stdout, "Potential coefficient counts\n multipole only:\n");
    viewprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2ndsav/total), num4th, 100*(num4thsav/total), 
	    numexact, 100*(numexactsav/total));
#endif
    total = num2nd + num4th + numexact;
#if MULDAT == ON
    viewprintf(stdout, " multipole plus adaptive:\n");
    viewprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2nd/total), num4th, 100*(num4th/total), 
	    numexact, 100*(numexact/total));
#endif
    viewprintf(stdout, "Percentage of multiplies done by multipole: %.3g%%\n",
	    100*(size*size - total)/(size*size));
    if(size*size == total) 
	viewprintf(stdout, "Warning: no multipole acceleration\n");
  }
}

double tilelength(nq)
charge *nq;
{
  return nq->max_diag;
}

void InitCalcpVars(void)
{
	num2nd=0;
	num4th=0;
	numexact=0;
	num2ndsav=0;
	num4thsav=0;
	numexactsav=0;
}
