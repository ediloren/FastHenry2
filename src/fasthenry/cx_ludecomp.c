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
 
*//* Complex LU decomposition routines */
/*
  - returned matrix has L below the diagonal, U above (GVL1 pg 58)
  - if allocate == TRUE ends up storing LU (could be a lot)
*/

#include "induct.h"

CX **cx_ludecomp(matin, size, allocate)
CX **matin;
int size;
int allocate;
{
  extern int fulldirops;
  CX factor, **mat, tmp;
  int i, j, k;

  if(allocate == TRUE) {
    /* allocate for LU matrix and copy A */
    MALLOC(mat, size, CX*, ON, IND);
    for(i = 0; i < size; i++) {
      MALLOC(mat[i], size, CX, ON, IND);
      for(j = 0; j < size; j++) mat[i][j] = matin[i][j];
    }
  }
  else mat = matin;

  for(k = 0; k < size-1; k++) {	/* loop on rows */
    if(mat[k][k].real == 0.0 && mat[k][k].imag == 0.0) {
      fprintf(stderr, "ludecomp: zero pivot\n");
      exit(0);
    }
    for(i = k+1; i < size; i++) { /* loop on remaining rows */
      /*factor = (mat[i][k] /= mat[k][k]);*/
      cx_div(tmp, mat[i][k], mat[k][k]);
      factor = mat[i][k] = tmp;
      fulldirops++;
      for(j = k+1; j < size; j++) { /* loop on remaining columns */
	/* mat[i][j] -= (factor*mat[k][j]);*/
	cx_mul(tmp, factor, mat[k][j]);
	cx_sub(mat[i][j], mat[i][j], tmp);
	fulldirops++;
      }
    }
  }
  return(mat);
}

/*
  For direct solution of Pq = psi, used if DIRSOL == ON or if preconditioning.
*/
void cx_lu_solve(mat, x, b, size)
CX **mat, *x, *b;
int size;
{
  extern int fulldirops;
  int i, j;
  CX tmp;

  /* copy rhs */
  if(x != b) for(i = 0; i < size; i++) x[i] = b[i];

  /* forward elimination */
  for(i = 0; i < size; i++) {	/* loop on pivot row */
    for(j = i+1; j < size; j++) { /* loop on elimnation row */
      /* x[j] -= mat[j][i]*x[i]; */
      cx_mul(tmp, mat[j][i], x[i]);
      cx_sub(x[j], x[j], tmp);
      fulldirops++;
    }
  }

  /* back substitution */
  for(i--; i > -1; i--) {		/* loop on rows */
    for(j = i+1; j < size; j++) { /* loop on columns */
      /* x[i] -= mat[i][j]*x[j]; */
      cx_mul(tmp, mat[i][j], x[j]);
      cx_sub(x[i], x[i], tmp);
      fulldirops++;
    }
    /* x[i] /= mat[i][i]; */
    cx_div(tmp, x[i], mat[i][i]);
    x[i] = tmp;
    fulldirops++;
  }
}

