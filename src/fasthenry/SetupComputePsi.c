/* this sets up the vectors to call Fastcap's ComputePsi */
/* It will be called twice for each coordinate direction.  Once for real
   and once for imaginary */

#include "induct.h"

/* Vs will contain the result,  Im is the 'q',  Size is the size of vectors. */
/* This will alter Im.  Im = Precond*Im */
SetupComputePsi(Vs, sys, Im, size, chglist, w, R, indsys)
CX *Vs, *Im;
int size;
ssystem *sys;
charge *chglist;
double w;  /* radian frequency */
double *R;  /* resistance vector */
SYS *indsys;
{
  extern double dirtime;
  double *q, *p;
  static CX *Ib = NULL, *Vb = NULL, *Vdirect = NULL, *ctemp;
  int branches;
  CX temp;
  MELEMENT *mtemp;
  charge *chg;
  int i, j;
  double rtemp;
  MELEMENT **Mtrans, **Mlist;
  double maxdiff,pdiff;
  int maxindx;
  int ind_opcnt_mult = 0, ind_opcnt_real = 0;

  branches = indsys->num_fils;
  Mtrans = indsys->Mtrans;
  Mlist = indsys->Mlist;

  if (Ib == NULL) {
    Ib = (CX *)MattAlloc(branches, sizeof(CX));
    Vb = (CX *)MattAlloc(branches, sizeof(CX));
    ctemp = (CX *)MattAlloc(size, sizeof(CX));
#ifndef NODEBUG
    Vdirect = (CX *)MattAlloc(branches, sizeof(CX));
#endif
  }

  for(i = 0; i < branches; i++)
    Vb[i] = CXZERO;

  q = sys->q;
  p = sys->p;
  ASSERT(size == indsys->num_mesh);

  if (indsys->precond_type == LOC) {
    multPrecond(indsys->Precond, Im, ctemp, size);
    for(i = 0; i < size; i++)
      Im[i] = ctemp[i];
  }
  else if (indsys->precond_type == SPARSE) 
    spSolve(indsys->sparMatrix, Im, Im);

  /* do  Ib = Mtrans*Im */
  for(i = 0; i < branches; i++) {
    Ib[i] = CXZERO;
    for(mtemp = Mtrans[i]; mtemp != NULL; mtemp = mtemp->mnext) {
      if (mtemp->sign == 1) 
	cx_add(Ib[i], Ib[i], Im[mtemp->filindex]);
      else
	cx_sub(Ib[i], Ib[i], Im[mtemp->filindex]);
    }
  }

  /* Evaluate M*L*Mt*Im = M*L*Ib using the multipole algorithm */

  /* Do all of the non-direct parts first */
  sys->DirectEval = FALSE;
  for(i = 0; i < 3; i++) {  /* for each of the coordinate directions */

    /* do the real part */
    for(chg = chglist; chg != NULL; chg = chg->next) {
      /* fill the pseudo-charge vector */
      q[chg->index] = Ib[chg->fil->filnumber].real*chg->fil->lenvect[i];
#if OPCNT == ON
      ind_opcnt_mult++;
#endif

    }
    computePsi(sys, q, p, branches, chglist);
    for(chg = chglist; chg != NULL; chg = chg->next) {
      /* add potential due to i direction */
      Vb[chg->fil->filnumber].real += p[chg->index]*chg->fil->lenvect[i]*MUOVER4PI;
#if OPCNT == ON
      ind_opcnt_mult++;
#endif
    }

    /* do the imaginary part */
    for(chg = chglist; chg != NULL; chg = chg->next) {
      /* fill the pseudo-charge vector */
      q[chg->index] = Ib[chg->fil->filnumber].imag*chg->fil->lenvect[i];
#if OPCNT == ON
      ind_opcnt_mult++;
#endif
    }
    computePsi(sys, q, p, branches, chglist);
    for(chg = chglist; chg != NULL; chg = chg->next) {
      /* add potential due to i direction */
      Vb[chg->fil->filnumber].imag += p[chg->index]*chg->fil->lenvect[i]*MUOVER4PI;
#if OPCNT == ON
      ind_opcnt_mult++;
#endif
    }
    
  }

  /* do the direct parts */
  sys->DirectEval = TRUE;

  /* do the real part of the Direct part */
  for(i = 1; i <= branches; i++)
    p[i] = 0;
  for(chg = chglist; chg != NULL; chg = chg->next) 
    /* fill the pseudo-charge vector */
    q[chg->index] = Ib[chg->fil->filnumber].real;

  /* starttimer; */
  mulDirect(sys);
  mulEval(sys);
  /* stoptimer; */
  dirtime += dtime;

  for(chg = chglist; chg != NULL; chg = chg->next) {
    /* add potential due to i direction */
    Vb[chg->fil->filnumber].real += p[chg->index];
  }
  
  /* do the imaginary part of the Direct part */
  for(i = 1; i <= branches; i++)
    p[i] = 0;
  for(chg = chglist; chg != NULL; chg = chg->next)
    /* fill the pseudo-charge vector */
    q[chg->index] = Ib[chg->fil->filnumber].imag;

  /* starttimer; */
  mulDirect(sys);
  mulEval(sys);
  /* stoptimer; */
  dirtime += dtime;

  for(chg = chglist; chg != NULL; chg = chg->next) { 
    /* add potential due to i direction */
    Vb[chg->fil->filnumber].imag += p[chg->index];
  }

  /* do Vs = M*Vb*jw */
  for(i = 0; i < size; i++) {
    Vs[i] = CXZERO;
    for(mtemp = Mlist[i]; mtemp != NULL; mtemp = mtemp->mnext)
      if (mtemp->sign == 1) 
	cx_add(Vs[i], Vs[i], Vb[mtemp->filindex]);
      else
	cx_sub(Vs[i], Vs[i], Vb[mtemp->filindex]);

    /* multiply by jw */
    rtemp = -Vs[i].imag*w;
    Vs[i].imag = Vs[i].real*w;
    Vs[i].real = rtemp;
  }

  /* add in M*R*Mt*Im = M*R*Ib */
  for(i = 0; i < size; i++) {
    for(mtemp = Mlist[i]; mtemp != NULL; mtemp = mtemp->mnext) {
     cx_scalar_mult(temp, mtemp->sign*R[mtemp->filindex], Ib[mtemp->filindex]);
     cx_add(Vs[i], Vs[i], temp);
#if OPCNT == ON
      ind_opcnt_mult+=2;
      ind_opcnt_real+=2;
#endif
    }
  }

#if OPCNT == ON
  printf("Inductance (mesh to branch) mults: %d\n",ind_opcnt_mult);
  printf("Just doing MRMtIm: %d\n",ind_opcnt_real);
   printops();
  exit(0);
#endif

#ifdef NODEBUG
  /* for debugging, compare to direct Vb = ZM Ib */
  realmatCXvec(Vdirect, indsys->Z, Ib, branches);
  maxdiff = 0; 
  maxindx = 0;
  for(i = 0; i < branches; i++) {
    if (cx_abs(Vb[i]) > 1e-23) {
      cx_sub(temp, Vdirect[i], Vb[i]);
      pdiff = cx_abs(temp )/cx_abs(Vb[i]) ;
    }
    else
      pdiff = cx_abs(Vb[i]);

    if (pdiff > maxdiff) {
      maxdiff = pdiff;
      maxindx = i;
    }
  }
  if (maxdiff < .3)
    printf("maxdiff: %g  Vb[%d]=%g  Vdirect[%d]=%g\n",
	   maxdiff,maxindx,cx_abs(Vb[maxindx]),maxindx,cx_abs(Vdirect[maxindx]));
  else
    printf("***maxdiff: %g  Vb[%d]=%g  Vdirect[%d]=%g***\n",
	   maxdiff,maxindx,cx_abs(Vb[maxindx]),maxindx,cx_abs(Vdirect[maxindx]));    


#endif
}

realmatCXvec(y, A, x, size)
CX *y, *x;
double **A;
int size;
{
  int i, j;
  CX temp;

  for (i = 0; i < size; i++) {
    y[i] = CXZERO;
    for(j = 0; j < size; j++) {
      cx_scalar_mult(temp, A[i][j], x[j]);
      cx_add(y[i], y[i], temp);
    }
  }
}

/* this function fixes Eval matrices which are computed directly */
/* This is necessary since direct mutual terms are not componentwise,
   but the multipole routines are called once for each component direction.
   Basically, componentwise multiplication will cause the elements
   to be multiplied by the dot product of the fil->lenvect vectors of 
   the two filaments.  This will divide that product out.  Also, MUOVER4PI
   must also be divided out 
*/

fixEvalDirect(qchgs, numqchgs, is_dummy, pchgs, numpchgs, mat)
charge **qchgs, **pchgs;
int numqchgs, numpchgs;
int *is_dummy;
double **mat;
{
  int i,j, k;
  double dotprod, magi, magj;
  double *lenvecti, *lenvectj;
  static double eps = EPS;

  for(i = 0; i < numpchgs; i++) {
    lenvecti = pchgs[i]->fil->lenvect;
    magi = 0;
    for(k = 0; k < 3; k++)
      magi += lenvecti[k]*lenvecti[k];
    for(j = 0; j < numqchgs; j++) {
      lenvectj = qchgs[j]->fil->lenvect;
      magj = dotprod = 0;
      for(k = 0; k < 3; k++) {
	magj += lenvectj[k]*lenvectj[k];
	dotprod += lenvecti[k]*lenvectj[k];
      }
      if (fabs(dotprod)/sqrt(magi*magj) > EPS) /* filaments aren't perpendicular */
	mat[i][j] = mat[i][j]/(dotprod*MUOVER4PI);
      else { /* if they are, mat[i][j] == 0.0, hopefully */
	if (mat[i][j] != 0.0)
	  printf("Warning: dot product = %lg < EPS, but mat[i][j] = %lg\n",dotprod, mat[i][j]);
      }
    }
  }
}
