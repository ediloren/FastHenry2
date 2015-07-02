/* This is a complex Generalized Minimal residual algorithm */

#include "induct.h"
#define EPSGMRES 1e-13    /* machine precision * 10 */

/*
mat - dense matrix to block precondition. Mat is replaced with  D^-1 Mat, 
except that D^-1 is on the block diagonal instead of the identity.
rhs - rhs to block precond.
size - matrix size.
vect - temporary vector of complex numbers (of length size).
pvect - temporary vector of ptrs to complex numbers (of length size).
numblocks - number of diag blocks.
startrows - beginning row of each diag block, 
startrows[0] = 0, startrows[numblocks] = size;
*/

int blockPrecond(mat, rhs, size, vect, pvect, numblocks, startrows)
CX **mat, *rhs;
CX *vect, **pvect;
int size, numblocks;
int *startrows;
{
int i, j;
CX **diagBlock;
  if(startrows[0] != 0) { printf("startrow not zero\n"); exit(0);}
  if(startrows[numblocks] != size) { printf("startrow not size\n"); exit(0);}

  for(i=0; i < numblocks; i++) {
    /* First invert the diagonal block. */
    for(j=startrows[i]; j < startrows[i+1]; j++) 
	pvect[j] = &(mat[j][startrows[i]]);
    diagBlock = &(pvect[startrows[i]]);
    cx_invert(diagBlock, startrows[i+1] - startrows[i]);

    /* Fix the startrows[i] thru startrows[i+1] elements of the rhs. */
    matMultVect(diagBlock, rhs, vect, startrows[i+1] - startrows[i]);
    for(j=startrows[i]; j < startrows[i+1]; j++) 
	rhs[j] = vect[j];

    /* form the startrows[i] thru startrows[i+1] rows of D^-1 A. */
    matMult(diagBlock, &(mat[startrows[i]]),startrows[i+1] - startrows[i], 
	    0, startrows[i], vect);
    matMult(diagBlock, &(mat[startrows[i]]), startrows[i+1] - startrows[i], 
	    startrows[i+1], startrows[numblocks], vect);

  }

}
      

/* 
  In-place inverts a matrix using guass-jordan.
*/
cx_invert(mat, size)
CX **mat;
int size;
{
  int i, j, k;
  CX normal, multiplier, tmp;

  for(i=0; i < size; i++) {
    /* First i^{th} column of A. */
    cx_div(normal, CXONE, mat[i][i]);
    for(j=0; j < size; j++) {
      cx_mul(tmp, mat[j][i], normal);
      mat[j][i] = tmp;
    }
    mat[i][i] = normal;

    /* Fix the backward columns. */
    for(j=0; j < size; j++) {
      if(j != i) {
	cx_mul(multiplier, CXMONE, mat[i][j]);
	for(k=0; k < size; k++) {
	  cx_mul(tmp, mat[k][i], multiplier);
	  if(k != i) cx_add(mat[k][j], mat[k][j], tmp);
	  else mat[k][j] = tmp;
	}
      }
    }
  }
}

/*
Multiplies mat1 * mat2, replacing mat2 with the product. 
*/
matMult(mat1, mat2, rows1, firstcol, lastcol, vtemp)
CX **mat1, **mat2, *vtemp;
int rows1, firstcol, lastcol;
{
int i, j, k;
CX *row, tmp;

  for(i=firstcol; i < lastcol; i++) {  /* Get the next column of mat2. */
    for(j=0; j < rows1; j++) {
      row = mat1[j];
      vtemp[j] = CXZERO;
      for(k=0; k < rows1; k++) {
	cx_mul(tmp, mat2[k][i], row[k]);
	cx_add(vtemp[j], vtemp[j], tmp);
      }
    }
    for(j=0; j < rows1; j++) mat2[j][i] = vtemp[j];
  }
}


/*
dvect = mat1 * svect.
*/
matMultVect(mat1, svect, dvect, rows)
CX **mat1, *svect, *dvect;
int rows;
{
  int i, j, k;
  CX *row, tmp, sum;

  for(j=0; j < rows; j++) {
    row = mat1[j];
    sum = CXZERO;
    for(k=0; k < rows; k++) {
      cx_mul(tmp, svect[k], row[k]);
      cx_add(sum, sum, tmp);
    }
    dvect[j] = sum;
  }
}
    
int
gmres(A, b, x0, inner, matvec, size, maxiters, tol, sys, chglist, w, R, indsys,
      cond)
  CX **A;
  CX *b;
  CX *x0;
  CX (*inner)();
  int (*matvec)();
  int size;
  int maxiters;
  double tol;
  SYS *indsys;
  ssystem *sys;
  charge *chglist;
  double w;  /* frequency */
  double *R;
  int cond;   /* conductor number for saving residual data */
{
  register int i, j, k;
  double rnorm, norm, length, blahnorm;
  int retval = -maxiters;
  CX hi, hip1;
  CX tmp1, tmp2;
  double zeta, rtmp1, rtmp2, r_real, r_imag, max_real, max_imag, absolute;
  static int lastmax=0, lastsize=0;
  static CX *c = NULL, *s = NULL, *g = NULL, *y = NULL;
  static CX *Rw, *Pw, *APw, *xlast, *x, *temp;
  static CX **bh, **bv;
  static CX *temp1, *temp2;

  if((size > lastsize) || (maxiters > lastmax)) {
    c = (CX *) MattAlloc(maxiters+1, sizeof(CX));
    s = (CX *) MattAlloc(maxiters+1, sizeof(CX));
    g = (CX *) MattAlloc(maxiters+1, sizeof(CX));
    y = (CX *) MattAlloc(maxiters+1, sizeof(CX));
    if (y == NULL) { printf("no space \n"); exit(1); }
    Rw = (CX *) MattAlloc(size, sizeof(CX));
    x = (CX *) MattAlloc(size, sizeof(CX));
    xlast = (CX *) MattAlloc(size, sizeof(CX));
    APw = (CX *) MattAlloc(size, sizeof(CX));
    Pw = (CX *) MattAlloc(size, sizeof(CX));
    temp1 = (CX *) MattAlloc(size, sizeof(CX));
    if (Pw == NULL){ printf("no space \n"); exit(1); }
  
    bh = (CX **) MattAlloc(maxiters+1, sizeof(CX*));
    for (i = 0; i <= maxiters; i++){
      bh[i] = (CX *) MattAlloc(maxiters+1, sizeof(CX));
    }
    if (bh[maxiters] == NULL){ printf("no space \n"); exit(1); }
    
    bv = (CX **) MattAlloc(maxiters + 1, sizeof(CX*));
    for (i = 0; i <= maxiters; i++){
      bv[i] = (CX *) MattAlloc(size, sizeof(CX));
    }
    if (bv[maxiters] == NULL){ printf("no space \n"); exit(1); }

    lastsize = size;
    lastmax = maxiters;
  }

  matvec(Rw, sys, x0, size, chglist, w, R, indsys);
  sub_3(Rw, b, Rw, size);
    
  for (i = 0; i <= maxiters; i++) {
    c[i] = s[i] = g[i] = y[i] = CXZERO;
    for (j = 0; j <= maxiters; j++)
      bh[i][j] = CXZERO;
  }
  rnorm = cx_abs(inner(Rw, Rw, size));
  rnorm = sqrt(rnorm);
  zeta = rnorm;
  if (indsys->opts->debug == ON)
    printf("zeta=%lg (initial rnorm)\n",zeta);

  set_scalar_mult_2(Pw, Rw, 1./rnorm, size);
  
  g[0].real = rnorm;
  g[0].imag = 0.0;
  r_real = r_imag = tol*2;

  for (k = 0; (k < maxiters) && (rnorm > (zeta*EPSGMRES))
              && ((rnorm > (zeta * tol)) || (r_real > tol) 
		  || (r_imag > tol)); k++) {

    /* Save p as the v{iter}. */
    for(i=0; i < size; i++) 
      bv[k][i] = Pw[i];

    matvec(APw, sys, Pw, size, chglist, w, R, indsys);

    for(i=0; i < size; i++) 
      Pw[i] = APw[i];

    for (j = 0; j <= k; j++) {
/*    hi = inner(APw, bv[j], size); */       /* h_jk = (A*v_k, v_j) */
      hi = inner(Pw, bv[j], size);        /* modified gram-schmidt */
      sub_cx_mult_2(Pw, bv[j], hi, size);  /* vhat_k+1 = vhat_k+1 - h_jk*v_j */
      bh[k][j] = hi;                       /* save h_jk in Hbar */
    }
    
    norm = cx_abs(inner(Pw, Pw, size));   /* norm = norm(vhat_k+1) */
    norm = sqrt(norm);
    if (norm != 0.0) 
      set_scalar_mult_1(Pw, 1./norm, size); /* Pw = vhat_k+1/norm */
    else
      /* We have the exact solution, and we don't want to cause divide by 0*/
      /*  We got exact solution because RHS happens to be a sum of 
          <# of iters> eigenvectors.  Not unlikely for freq=0 and conductors
          are decoupled */
      set_scalar_mult_1(Pw, 0.0, size); /* Pw = vhat_k+1/norm */
      
    bh[k][k+1].real = norm;               /* h_k+1,k = norm(vhat_k+1) */
    bh[k][k+1].imag = 0.0;
    
    /* rotate bh[k] by previous rotations used to make bh upper triangular */
    for (i = 0; i < k; i++) {
      hi = bh[k][i];
      hip1 = bh[k][i+1];
/* This doesn't work anymore.  rotation matrix should have conjugate first row 
      cx_mul(tmp1, c[i], hi);
      cx_mul(tmp2, s[i], hip1);
*/
      cx_conj_mul(tmp1, hi, c[i]);
      cx_conj_mul(tmp2, hip1, s[i]);
      cx_sub(bh[k][i],tmp1,tmp2);

      cx_mul(tmp1, c[i], hip1);
      cx_mul(tmp2, s[i], hi);
      cx_add(bh[k][i+1],tmp1,tmp2);
    }
    
    hi = bh[k][k];         /* r in GMRES paper */
    hip1 = bh[k][k+1];     /* h in GMRES paper */

    rtmp1 = cx_abs(hi);    /* || r || */
    rtmp2 = cx_abs(hip1);  /* || h || */
    length = sqrt(rtmp1 * rtmp1 + rtmp2 * rtmp2);

    /* compute next rotation */
    c[k].real = hi.real / length;
    c[k].imag = hi.imag / length;
    s[k].real = -hip1.real / length;
    s[k].imag = -hip1.imag / length;

/*        need conjugate now
    cx_mul(tmp1, c[k], hi);
    cx_mul(tmp2, s[k], hip1);
*/
    cx_conj_mul(tmp1, hi, c[k]);
    cx_conj_mul(tmp2, hip1, s[k]);
    cx_sub(bh[k][k],tmp1,tmp2);

    cx_mul(tmp1, c[k], hip1);
    cx_mul(tmp2, s[k], hi);
    cx_add(bh[k][k+1],tmp1,tmp2);    /* this should be zero */


    hi = g[k];
/*    cx_mul(g[k], c[k], hi);  need conjugate */
    cx_conj_mul(g[k], hi, c[k]);
    cx_mul(g[k+1], s[k], hi);
    
    rnorm = cx_abs(g[k+1]);

    /* Compute solution, note, bh is bh[col][row]. */
    /* backsolve for y: g = bh*y */
    for(i=0; i <= k; i++) y[i] = g[i];
    for(i = k; i >= 0; i--) {
      cx_div(hi, y[i], bh[i][i]);
      y[i] = hi;
      for(j = i-1; j >= 0; j--) {
	cx_mul(hi, bh[i][j],y[i]);
	cx_sub(y[j], y[j], hi); 
      }
    }
    /* x = x0 + bv*y */
    for(i=0; i < size; i++) {
      xlast[i] = x[i];
      x[i] = x0[i];
      for(j=0; j <= k; j++) {
	cx_mul(hi, y[j], bv[j][i]);
	cx_add(x[i], x[i], hi);
      }
    }
    calc_other_norm(x, xlast, size, indsys->opts->abs_tol, 
		    &r_real, &r_imag, &max_real, &max_imag);

    if (indsys->opts->debug == ON) {
      printf("%d rnorm=%g real_part=%lg imag_part=%lg\n", 
	     k+1,rnorm, r_real, r_imag);
      fflush(stdout);
    }
    else {
      printf("%d ",k+1);
      fflush(stdout);
    }

    if (1 == 1) {  /* lets always dump residual data */
      indsys->resids[cond][k] = rnorm;
      indsys->resid_real[cond][k] = r_real;
      indsys->resid_imag[cond][k] = r_imag;
    }

#ifdef DEBUG
    matvec(APw, sys, x, size, chglist, w, R, indsys);
    for(i=0; i < size; i++) {
      temp1[i].real = b[i].real - APw[i].real;
      temp1[i].imag = b[i].imag - APw[i].imag;
    }
    blahnorm = cx_abs(inner(temp1, temp1, size));
    blahnorm = sqrt(blahnorm);
    printf("check norm =%g\n", blahnorm); 
#endif

  }

  /* Decrement from the last increment. */
  k--;

  if (1 == 1) {  /* lets always do it */
    /* save number of iters data */
    indsys->niters[cond] = k;
  }
  
  if (k == maxiters - 1) 
    printf("Warning: Convergence not reached. Stopping iteration...\n");
  if (indsys->opts->debug == ON)
    printf("Total GMRES iters = %d\n", k+1);
  else
    printf("\n");
  
  /* i did this already.  But I'll leave it in case I ever remove 
     the backsolving in the main loop */

  /* Compute solution, note, bh is bh[col][row]. */
/* backsolve g = bh*y */
  for(i=0; i <= k; i++) y[i] = g[i];
  for(i = k; i >= 0; i--) {
    cx_div(hi, y[i], bh[i][i]);
    y[i] = hi;
    for(j = i-1; j >= 0; j--) {
      cx_mul(hi, bh[i][j],y[i]);
      cx_sub(y[j], y[j], hi); 
    }
  }
/* x0 = x0 + bv*y */
  for(i=0; i < size; i++) {
    for(j=0; j <= k; j++) {
      cx_mul(hi, y[j], bv[j][i]);
      cx_add(x0[i], x0[i], hi);
    }
  }

  /* Note, x0 now has the solution. */
  return retval;
}


sub_3(z,x,y,size)
  CX *z, *x, *y;
  int size;
{
  int i;
  CX tmp;
  
  for (i = 0; i < size; i ++) {
    cx_sub(tmp,x[i],y[i]);
    z[i] = tmp;
  }
  
}


set_scalar_mult_1(x,alpha,size)
  CX *x;
  double alpha;
  int size;
{
  int i;
  
  for (i = 0; i < size; i ++) {
    x[i].real *= alpha;
    x[i].imag *= alpha;
  }
  
}


set_scalar_mult_2(x,y,alpha,size)
  CX *x, *y;
  double alpha;
  int size;
{
  int i;
  
  for (i = 0; i < size; i ++) {
    x[i].real = alpha * y[i].real;
    x[i].imag = alpha * y[i].imag;
  }
  
}


add_scalar_mult_2(x,y,alpha,size)
  CX *x, *y;
  double alpha;
  int size;
{
  int i;
  
  for (i = 0; i < size; i ++) {
    x[i].real += alpha * y[i].real;
    x[i].imag += alpha * y[i].imag;
  }
  
}


add_cx_mult_2(x,y,alpha,size)
  CX *x, *y;
  CX alpha;
  int size;
{
  int i;
  CX tmp1;
  CX tmp2;
  
  for (i = 0; i < size; i ++) {
    cx_mul(tmp1,alpha,y[i]);
    cx_add(tmp2,x[i],tmp1);
    x[i] = tmp2;
  }
  
}

sub_cx_mult_2(x,y,alpha,size)
  CX *x, *y;
  CX alpha;
  int size;
{
  int i;
  CX tmp1;
  CX tmp2;
  
  for (i = 0; i < size; i ++) {
    cx_mul(tmp1,alpha,y[i]);
    cx_sub(tmp2,x[i],tmp1);
    x[i] = tmp2;
  }
  
}

CX
inner(x, y, size)
  CX *x, *y;
  int size;
{
  int i;
  CX tmp1, tmp2;

  tmp2 = CXZERO;
  for (i = 0; i < size; i++) {
    cx_conj_mul(tmp1, x[i], y[i]);
    cx_add(tmp2,tmp2,tmp1);
  }

  return tmp2;
}

#if 1==0
/* this multiplies the preconditioned matrix by a vector */
#ifdef NOMULTI
/* it uses the global variables because I have no patience */
#endif NOMULTI
/* (The blockdiagonal pieces are the identity) */

int
#if !defined ( NOMULTI )
oldmatvec(y, A, x, size,startrows, numblocks)
#else NOMULTI
matvec(y, A, x, size,startrows, numblocks)
#endif NOMULTI
  CX *x, *y;
  CX **A;
  int *startrows, numblocks;
{
  int i,j, row;
  CX temp;

#if !defined(NOPRECOND)
  for(i = 0; i < numblocks; i++) {
    for(row = startrows[i]; row < startrows[i+1]; row++) {
      y[row] = CXZERO;
      for(j = 0; j < startrows[i]; j++) {
	cx_mul(temp, A[row][j], x[j]);
	cx_add(y[row], y[row], temp);
      }
      cx_add(y[row], y[row], x[row]);
      for(j = startrows[i+1]; j < size; j++) {
	cx_mul(temp, A[row][j], x[j]);
	cx_add(y[row], y[row], temp);
      }
    }
  }
#else
/*  The old way */
  for (i = 0; i < size; i++) {
    y[i] = CXZERO;
    for(j = 0; j < size; j++) {
      cx_mul(temp, A[i][j], x[j]);
      cx_add(y[i], y[i], temp);
    }
  }
#endif
  
}
#endif

veryoldmatvec(y, A, x, size)
  CX *x, *y;
  CX **A;
{
  int i, j;
  CX temp;

  for (i = 0; i < size; i++) {
    y[i] = CXZERO;
    for(j = 0; j < size; j++) {
      cx_mul(temp, A[i][j], x[j]);
      cx_add(y[i], y[i], temp);
    }
  }
}

disMat(A, size)
CX **A;
int size;
{
  int i, j;

  for(i=0; i < size; i++) {
    printf("\n");
    for(j=0; j < size; j++) {
      printf("r=%g i=%g ", A[i][j].real, A[i][j].imag);
    }
  }
  printf("\n");
}

/* multiplies the vector by just the diagonal blocks */
/* That is, it multiplies by the preconditioner */
multvec2(A, x, y, startrows, numblocks)
int *startrows, numblocks;
  CX *x, *y;
  CX **A;
{
  int i,j, row;
  CX temp;

  for(i = 0; i < numblocks; i++) {
    for(row = startrows[i]; row < startrows[i+1]; row++) {
      y[row] = CXZERO;
      for(j = startrows[i]; j < startrows[i+1]; j++) {
	cx_mul(temp, A[row][j], x[j]);
	cx_add(y[row], y[row], temp);
      }
    }
  }
}

/* Vs holds the result of (indsys->MtZM)*(indsys->Precond)*Im */
/* All other arguments are not used.   MK 10/92 */
directmatvec(Vs, sys, Im, size, chglist, w, R, indsys)
CX *Vs, *Im;
int size;
ssystem *sys;
charge *chglist;
double w;  /* radian frequency */
double *R;  /* resistance vector */
SYS *indsys;
{
  int i;

  if (indsys->precond_type == LOC) {
    multPrecond(indsys->Precond, Im, Vs, size);
    for(i = 0; i < size; i++)
      Im[i] = Vs[i];
  }
  else if (indsys->precond_type == SPARSE) 
    spSolve(indsys->sparMatrix, Im, Im);

  veryoldmatvec(Vs, indsys->MtZM, Im, size);
}

calc_other_norm(x, xlast, size, abs_tol, r_real, r_imag, max_real, max_imag)
     double abs_tol, *r_real, *r_imag, *max_real, *max_imag;
     CX *x, *xlast;
     int size;
{
  double r_max_diff, r_max_val, i_max_diff, i_max_val, r_abstol, i_abstol;
  int i;
  double temp, atol;

  r_max_diff = r_max_val = i_max_diff = i_max_val = 0;

  for(i = 0; i < size; i++) {
    if (fabs(x[i].real) > r_max_val) r_max_val = fabs(x[i].real);
    if (fabs(x[i].imag) > i_max_val) i_max_val = fabs(x[i].imag);
  }

  r_abstol = abs_tol*r_max_val;
  i_abstol = abs_tol*i_max_val;
/*   atol = EPSGMRES*MAX(r_max_val, i_max_val); */
  if (i_abstol == 0) i_abstol = 1e-200;
  if (r_abstol == 0) i_abstol = 1e-200;
  
  for(i = 0; i < size; i++) {
    temp = fabs(x[i].real - xlast[i].real)/(r_abstol + fabs(x[i].real));
    if (temp > r_max_diff) r_max_diff = temp;

    temp = fabs(x[i].imag - xlast[i].imag)/(i_abstol + fabs(x[i].imag));
    if (temp > i_max_diff) i_max_diff = temp;
  }

  *max_real = r_max_val;
  *r_real = r_max_diff;

  *max_imag = i_max_val;
  *r_imag = i_max_diff;
}
