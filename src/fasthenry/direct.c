/* # ***** sort to /src/direct
   # ***** */
#include "mulGlobal.h"

double **Q2PDiag(chgs, numchgs, is_dummy, calc)
charge **chgs;
int numchgs, *is_dummy;
int calc;
{
  double **mat;
  int i, j;
  double calcp();

  /* Allocate storage for the potential coefficients. */
  CALLOC(mat, numchgs, double*, ON, AQ2PD);
  for(i=0; i < numchgs; i++) CALLOC(mat[i], numchgs, double, ON, AQ2PD);

  if(calc) {
    /* Compute the potential coeffs. */
    /* - exclude dummy panels when they would need to carry charge
       - exclude dielec i/f panels when they would lead to evals at their
         centers (only if using two-point flux-den-diff evaluations) */
    for(i=0; i < numchgs; i++) { 
#if NUMDPT == 2
      if(chgs[i]->dummy);	/* don't check surface of a dummy */
      else if(chgs[i]->surf->type == DIELEC || chgs[i]->surf->type == BOTH)
	  continue;
#endif
      for(j=0; j < numchgs; j++) { /* need to have charge on them */
#if SKIPQD == ON
	if(chgs[j]->pos_dummy == chgs[i] || chgs[j]->neg_dummy == chgs[i])
	    continue;
#endif
	if(!is_dummy[j]) mat[i][j] = calcp(chgs[j], chgs[i], NULL);
      }
    }
  }

#if DSQ2PD == ON
  dispQ2PDiag(mat, chgs, numchgs, is_dummy);
#endif

  return(mat);
}

double **Q2P(qchgs, numqchgs, is_dummy, pchgs, numpchgs, calc)
charge **qchgs, **pchgs;
int numqchgs, numpchgs;
int *is_dummy;
int calc;
{
  double **mat;
  int i, j;
  double calcp();

  /* Allocate storage for the potential coefficients. P rows by Q cols. */
  CALLOC(mat, numpchgs, double*, ON, AQ2P);
  for(i=0; i < numpchgs; i++) {
    CALLOC(mat[i], numqchgs, double, ON, AQ2P);
    if(calc) {
      /* exclude:
         - dummy panels in the charge list
	 - dielectric i/f panels in the eval list (if doing 2-point E's)*/
#if NUMDPT == 2
      if(pchgs[i]->dummy);	/* don't check the surface of a dummy */
      else if(pchgs[i]->surf->type == DIELEC || pchgs[i]->surf->type == BOTH)
	  continue;
#endif
      for(j=0; j < numqchgs; j++) { /* only dummy panels in the charge list */
	if(!is_dummy[j])          /* (not the eval list) are excluded */
	    mat[i][j] = calcp(qchgs[j], pchgs[i], NULL);
	                                 /* old: pchgs[i],qchgs[j] */
      }
    }
  }

#if DISQ2P == ON
  dispQ2P(mat, qchgs, numqchgs, is_dummy, pchgs, numpchgs);
#endif

  return(mat);
}

/*
  used only in conjunction with DMPMAT == ON  and DIRSOL == ON
  to make 1st directlist mat = full P mat
*/
double **Q2Pfull(directlist, numchgs)
int numchgs;
cube *directlist;
{
  int i, j, fromp, fromq, top, toq;
  double **mat, calcp();
  cube *pq, *pp;
  charge **pchgs, **qchgs, *eval;

  /* allocate room for full P matrix */
  MALLOC(mat, numchgs, double*, ON, AQ2P);
  for(i = 0; i < numchgs; i++) MALLOC(mat[i], numchgs, double, ON, AQ2P);

  /* load the matrix in the style of Q2P() - no attempt to exploit symmetry */
  /* calculate interaction between every direct list entry and entire dlist */
  for(pp = directlist; pp != NULL; pp = pp->dnext) {
    pchgs = pp->chgs;
    fromp = pchgs[0]->index - 1; /* row index range */
    top = fromp + pp->upnumeles[0];
    for(pq = directlist; pq != NULL; pq = pq->dnext) {
      qchgs = pq->chgs;
      fromq = qchgs[0]->index - 1; /* column index range */
      toq = fromq + pq->upnumeles[0];

      for(i = fromp; i < top; i++) {
	for(j = fromq; j < toq; j++) { 
	  eval = qchgs[j-fromq];
	  mat[i][j] = calcp(pchgs[i-fromp],eval, NULL);
	}
      }

    }
  }
  return(mat);
}

/*
  - returned matrix has L below the diagonal, U above (GVL1 pg 58)
  - if allocate == TRUE ends up storing P and LU (could be a lot)
*/
double **ludecomp(matin, size, allocate)
double **matin;
int size;
int allocate;
{
  extern int fulldirops;
  double factor, **mat;
  int i, j, k;

  if(allocate == TRUE) {
    /* allocate for LU matrix and copy A */
    MALLOC(mat, size, double*, ON, AMSC);
    for(i = 0; i < size; i++) {
      MALLOC(mat[i], size, double, ON, AMSC);
      for(j = 0; j < size; j++) mat[i][j] = matin[i][j];
    }
  }
  else mat = matin;

  for(k = 0; k < size-1; k++) {	/* loop on rows */
    if(mat[k][k] == 0.0) {
      fprintf(stderr, "ludecomp: zero piovt\n");
      exit(0);
    }
    for(i = k+1; i < size; i++) { /* loop on remaining rows */
      factor = (mat[i][k] /= mat[k][k]);
      fulldirops++;
      for(j = k+1; j < size; j++) { /* loop on remaining columns */
	mat[i][j] -= (factor*mat[k][j]);
	fulldirops++;
      }
    }
  }
  return(mat);
}

/*
  For direct solution of Pq = psi, used if DIRSOL == ON or if preconditioning.
*/
void solve(mat, x, b, size)
double **mat, *x, *b;
int size;
{
  extern int fulldirops;
  int i, j;

  /* copy rhs */
  if(x != b) for(i = 0; i < size; i++) x[i] = b[i];

  /* forward elimination */
  for(i = 0; i < size; i++) {	/* loop on pivot row */
    for(j = i+1; j < size; j++) { /* loop on elimnation row */
      x[j] -= mat[j][i]*x[i];
      fulldirops++;
    }
  }

  /* back substitution */
  for(i--; i > -1; i--) {		/* loop on rows */
    for(j = i+1; j < size; j++) { /* loop on columns */
      x[i] -= mat[i][j]*x[j];
      fulldirops++;
    }
    x[i] /= mat[i][i];
    fulldirops++;
  }
}

/* 
  In-place inverts a matrix using guass-jordan.
  - is_dummy[i] = 0 => ignore row/col i
*/
void invert(mat, size, reorder)
double **mat;
int size, *reorder;
{
  int i, j, k, best;
  double normal, multiplier, bestval, nextbest;
/*
  matlabDump(mat,size,"p");
*/
  for(i=0; i < size; i++) {
    best = i;
    bestval = ABS(mat[i][i]);
    for(j = i+1; j < size; j++) {
      nextbest = ABS(mat[i][j]);
      if(nextbest > bestval) {
	best = j;
	bestval = nextbest;
      }
    }

    /* If reordering, find the best pivot. */
    if(reorder != NULL) {
      reorder[i] = best;
      if(best != i) {
	for(k=0; k < size; k++) {
	  bestval = mat[k][best];
	  mat[k][best] = mat[k][i];
	  mat[k][i] = bestval;
	}
      }
    }

    /* First i^{th} column of A. */
    normal = 1.0 / mat[i][i];
    for(j=0; j < size; j++) {
      mat[j][i] *= normal;
    }
    mat[i][i] = normal;

    /* Fix the backward columns. */
    for(j=0; j < size; j++) {
      if(j != i) {
	multiplier = -mat[i][j];
	for(k=0; k < size; k++) {
	  if(k != i) mat[k][j] += mat[k][i] * multiplier;
	  else mat[k][j] = mat[k][i] * multiplier;
	}
      }
    }
  }

  /* Unravel the reordering, starting with the last column. */
  if(reorder != NULL) {
    for(i=size-2; i >= 0; i--) {
      if(reorder[i] != i) {
	for(k=0; k < size; k++) {
	  bestval = mat[k][i];
	  mat[k][reorder[i]] = mat[k][i];
	  mat[k][i] = bestval;
	}
      }
    }
  }
/*
  matlabDump(mat,size,"c");
*/

}

/*
  Used in conjuction with invert() to remove dummy row/columns
   comp_rows = TRUE => remove rows corresponding to ones in is_dummy[]
   comp_rows = FALSE => remove cols corresponding to ones in is_dummy[]
   comp_rows = BOTH => remove both rows and columns
   returns number of rows/cols in compressed matrix
*/
int compressMat(mat, size, is_dummy, comp_rows)
double **mat;
int size, *is_dummy, comp_rows;
{
  static int *cur_order;
  static int cur_order_array_size = 0;
  int cur_order_size, i, j, k;
  
  if(cur_order_array_size < size) {
    CALLOC(cur_order, size, int, ON, AMSC);
  }
  
  /* figure the new order vector (cur_order[i] = index of ith row/col) */
  for(i = cur_order_size = 0; i < size; i++) {
    if(!is_dummy[i]) cur_order[cur_order_size++] = i;
  }

  if(comp_rows == TRUE || comp_rows == BOTH) {
    /* compress by removing rows from the matrix */
    for(i = 0; i < cur_order_size; i++) {
      if((k = cur_order[i]) == i) continue; /* if not reordered */
      for(j = 0; j < size; j++) { /* copy the row to its new location */
	mat[i][j] = mat[k][j];
      }
    }
  }
  if(comp_rows == FALSE || comp_rows == BOTH) {
    /* compress by removing columns from the matrix */
    for(j = 0; j < cur_order_size; j++) {
      if((k = cur_order[j]) == j) continue; /* if not reordered */
      for(i = 0; i < size; i++) { /* copy the col to its new location */
	mat[i][j] = mat[i][k];
      }
    }
  }
  return(cur_order_size);
}

/*
  Used in conjuction with invert() to add dummy row/columns
   exp_rows = TRUE => add rows corresponding to ones in is_dummy[]
   exp_rows = FALSE => add cols corresponding to ones in is_dummy[]
   exp_rows = BOTH => add rows and columns
*/
void expandMat(mat, size, comp_size, is_dummy, exp_rows)
double **mat;
int size, *is_dummy, exp_rows, comp_size;
{
  int i, j, k, next_rc;

  if(exp_rows == TRUE || exp_rows == BOTH) {
    next_rc = comp_size - 1;
    /* add rows to the matrix starting from the bottom */
    for(i = size - 1; i >= 0; i--) {
      if(is_dummy[i]) {		/* zero out dummy row */
	for(j = 0; j < size; j++) mat[i][j] = 0.0;
      }
      else {			/* copy the row from its compressed location */
	for(j = 0; j < size; j++) mat[i][j] = mat[next_rc][j];
	next_rc--;
      }
    }
  }
  if(exp_rows == FALSE || exp_rows == BOTH) {
    next_rc = comp_size - 1;
    /* add columns to the matrix starting from the right */
    for(j = size - 1; j >= 0; j--) {
      if(is_dummy[j]) {		/* zero out dummy column */
	for(i = 0; i < size; i++) mat[i][j] = 0.0;
      }
      else {			/* copy the col from its compressed location */
	for(i = 0; i < size; i++) mat[i][j] = mat[i][next_rc];
	next_rc--;
      }
    }
  }
}

/* 
Checks to see if the matrix has the M-matrix sign pattern and if
it is diagonally dominant. 
*/
matcheck(mat, rows, size)
double **mat;
int rows, size;
{
  double rowsum;
  int i, j;

  for(i = rows - 1; i >= 0; i--) {
    for(rowsum = 0.0, j = size - 1; j >= 0; j--) {
      if((i != j)  && (mat[i][j] > 0.0)) {
	printf("violation mat[%d][%d] =%g\n", i, j, mat[i][j]);
      }
      if(i != j) rowsum += ABS(mat[i][j]);
    }
    printf("row %d diag=%g rowsum=%g\n", i, mat[i][i], rowsum);
    if(rowsum > mat[i][i]) {
      for(j = size - 1; j >= 0; j--) {
	printf("col%d = %g ", j, mat[i][j]);
      }
      printf("\n");
    }
  }
}


matlabDump(mat, size, name)
double **mat;
int size;
char *name;
{
FILE *foo;
int i,j;
char fname[100];

  sprintf(fname, "%s.m", name);
  foo = fopen(fname, "w");
  fprintf(foo, "%s = [\n", name);
  for(i=0; i < size; i++) {
    for(j=0; j < size; j++) {
      fprintf(foo, "%.10e  ", mat[i][j]);
    }
    fprintf(foo, "\n");
  }
  fprintf(foo, "]\n");
}








