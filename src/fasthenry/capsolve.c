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
 
*//* # ***** sort to /src/main
   # ***** */
#include "mulGlobal.h"

/* 
ComputePsi computes the potential from the charge vector, or may
include a preconditioner.  It is assumed that the vectors for the
charge and potential have already been set up and that the potential
vector has been zeroed.  ARBITRARY VECTORS CAN NOT BE USED.
*/

computePsi(sys, q, p, size, chglist)
ssystem *sys;
double *q, *p;
int size;
charge *chglist;
{
  extern double dirtime, uptime, downtime, evaltime, prectime;
  extern int real_size;
  int i;

  ASSERT(p == sys->p);
  ASSERT(q == sys->q);

  for(i=1; i <= size; i++) p[i] = 0;

#if PRECOND != NONE
  starttimer;
  mulPrecond(sys, PRECOND);
  stoptimer;
  prectime += dtime;
#endif

#if EXPGCR == ON
  blkCompressVector(q+1, size, real_size, sys->is_dummy+1);
  blkAqprod(p+1, q+1, real_size, sqrmat);	/* offset since index from 1 */
  blkExpandVector(p+1, size, real_size); /* ap changed to p, r chged to q */
  blkExpandVector(q+1, size, real_size); /*    7 Oct 91 */
#else
/*  moved into SetupComputePsi since it only is done once 
  starttimer;
  mulDirect(sys);
  stoptimer;
  dirtime += dtime;
*/

  starttimer;
  mulUp(sys);
  stoptimer;
  uptime += dtime;

#if DUPVEC == ON
  dumpLevOneUpVecs(sys);
#endif

#if DNTYPE == NOSHFT
  mulDown(sys);		/* do downward pass without local exp shifts */
#endif

#if DNTYPE == GRENGD
  mulDown(sys);	       	/* do heirarchical local shift dwnwd pass */
#endif
  stoptimer;
  downtime += dtime;

  starttimer;
#if MULTI == ON
  mulEval(sys);		/* evaluate either locals or multis or both */
#endif
  stoptimer;
  evaltime += dtime;

#if DMPCHG == LAST
  fprintf(stdout, "\nPanel potentials divided by areas\n");
  dumpChgDen(stdout, p, chglist);
  fprintf(stdout, "End panel potentials\n");
#endif

  /* convert the voltage vec entries on dielectric i/f's into eps1E1-eps2E2 */
/*  compute_electric_fields(sys, chglist); */

#if 1==0
#if OPCNT == ON
  printops();
  exit(0);
#endif				/* OPCNT == ON */
#endif

#endif				/* EXPGCR == ON */
}
