/* Complex LU decomposition routines */
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

