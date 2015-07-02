/*!\page LICENSE LICENSE
 
Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of
Technology, hereafter designated as the Copyright Owners.
 
License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:
 
1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.
 
2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.
 
3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.
 
LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.
 
*//* this calls the routine to calculate the filament-filament 
   interaction exactly */

#include "induct.h"

static int num2nd=0, num4th=0, numexact=0;
static int num2ndsav=0, num4thsav=0, numexactsav=0;

double calcp(pchg1, pchg2, pfd)
charge *pchg1, *pchg2;
double *pfd;   /* left over from fastcap */
{
  double mutual(), selfterm();

  if (pfd != NULL)
    fprintf(stderr, "calcp: I don't know what to do with pfd!=NULL\n");

  if (pchg1->fil->filnumber == pchg2->fil->filnumber) 
    /* self term */
    return selfterm(pchg1->fil);
  else
    /* calculate mutual inductance of the two filaments */
    return mutual(pchg1->fil, pchg2->fil);
} 

/* from the fastcap calcp */
dumpnums(flag, size)
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
    fprintf(stdout, "Potential coefficient counts\n multipole only:\n");
    fprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2ndsav/total), num4th, 100*(num4thsav/total), 
	    numexact, 100*(numexactsav/total));
#endif
    total = num2nd + num4th + numexact;
#if MULDAT == ON
    fprintf(stdout, " multipole plus adaptive:\n");
    fprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2nd/total), num4th, 100*(num4th/total), 
	    numexact, 100*(numexact/total));
#endif
    fprintf(stdout, "Percentage of multiplies done by multipole: %.3g%%\n",
	    100*(size*size - total)/(size*size));
    if(size*size == total) 
	fprintf(stdout, "Warning: no multipole acceleration\n");
  }
}

double tilelength(nq)
charge *nq;
{
  return nq->max_diag;
}
