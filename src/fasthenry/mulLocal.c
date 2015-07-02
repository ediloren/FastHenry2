/* # ***** sort to /src/main
   # ***** */
#include "mulGlobal.h"

/*
  globals used for temporary storage
*/
double **facFrA;		/* array of factorial fractions */
double *cosmkB;			/* array used to look up cos[(m+-k)beta] */
double *sinmkB;			/* array used to look up sin[(m+-k)beta] */

/*
  initializes the factorial fraction array used in M2L, L2L matrix calculation
*/
void evalFacFra(array, order)
int order;			/* array is 2*order+1 x 2*order+1 */
double **array;			/* array[num][den] = num!/den! */
{
  int d, i;
  for(i = 0; i <= 2*order; i++) {
    array[i][i] = 1.0; /* do main diagonal */
    if(i > 0 && i < 2*order) array[i+1][i] = i+1; /* do first sub diagonal */
  }
  for(d = 3; d <= 2*order; d++) { /* loop on lower triangular rows */
    for(i = 1; i < d-1; i++) {	/* loop on columns */
      array[d][i] = array[d-1][i]*array[d][d-1];
    }
  }
  /* invert lower part entries and copy to top */
  for(d = 2; d <= 2*order; d++) {
    for(i = 1; i <= d-1; i++) {
      array[i][d] = 1/array[d][i];
    }
  }
  /* copy 1st row and column from computed values */
  for(d = 1; d <= 2*order; d++) {
    array[0][d] = array[1][d];
    array[d][0] = array[d][1];
  }

#if DISFAC == ON
  fprintf(stdout, "FACTORIAL FRACTION ARRAY:\n");
  dumpMat(array, 2*order+1, 2*order+1);
#endif

}

/*
  initializes sqrt((m+n)!/(n-m)!) lookup table (for L2L)
*/
void evalSqrtFac(arrayout, arrayin, order)
int order;
double **arrayout, **arrayin;
{
  int n, m;			/* arrayout[n][m] = sqrt((m+n)!/(n-m)!) */

  /* set up first column, always ones */
  for(n = 0; n < order+1; n++) arrayout[n][0] = 1.0;

  /* set up lower triangular (n+m)!/(n-m)! */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      arrayout[n][m] = sqrt(arrayin[n][m]);
    }
  }

#if DISSFA == ON
  fprintf(stdout, "SQUARE ROOT FACTORIAL ARRAY:\n");
  dumpMat(arrayout, order+1, order+1);
#endif

}


/*
  initializes cos[(m+-k)beta] and sin[(m+-k)beta] lookup tables (M2L and L2L)
*/
void evalSinCos(beta, order)
int order;
double beta;
{
  int i;
  double temp = beta;

  for(i = 1; i <= 2*order; beta += temp, i++) {
    sinmkB[i] = sin(beta);
    cosmkB[i] = cos(beta);
  }
}

/*
  looks up sin[(m+-k)beta]
*/
double sinB(sum)
int sum;
{
  if(sum < 0) return(-sinmkB[abs(sum)]);
  else return(sinmkB[sum]);
}

/*
  looks up cos[(m+-k)beta]
*/
double cosB(sum)
int sum;
{
  return(cosmkB[abs(sum)]);
}

/* 
  Used for all but no local downward pass. 
*/
double **mulMulti2Local(x, y, z, xp, yp, zp, order)
int order;
double x, y, z, xp, yp, zp;	/* multipole and local cube centers */
{
  int i, j, k, n, m;
  int terms = multerms(order);	/* the number of non-zero moments */
  int ct = costerms(order);	/* the number of non-zero cos (bar) moments */
  double **mat;			/* the transformation matrix */
  double rho, cosA, beta;	/* spher. position of multi rel to local */
  double rhoJ, rhoN;		/* rho^j and (-1)^n*rho^(n+1) in main loop */
  double rhoFac;		/* = rhoJ*rhoN intermediate storage */
  double temp1, temp2, temp3;
  double iPwr(), sinB(), cosB();
  extern double *tleg, *Ir, *Irn, *phi, *Mphi; /* external temporary storage */

  /* allocate the multi to local transformation matrix */
  CALLOC(mat, terms, double*, ON, AM2L);
  for(i = 0; i < terms; i++)
      CALLOC(mat[i], terms, double, ON, AM2L);

  /* find relative spherical coordinates */
  xyz2sphere(x, y, z, xp, yp, zp, &rho, &cosA, &beta);

  /* generate legendre function evaluations */
  evalLegendre(cosA, tleg, 2*order); /* multi->loc needs 2x legendres */

  /* generate sin[(m+-k)beta] and cos[(m+-k)beta] look up arrays */
  /*  other lookup arrays generated in mulMultiAlloc() */
  evalSinCos(beta, order);

  /* generate multi to local transformation matrix; uses NB12 pg30 */
  /*  rhoFac factor divides could be reduced to once per loop */
  for(j = 0, rhoJ = 1.0; j <= order; rhoJ *= rho, j++) {
    for(k = 0; k <= j; k++) {	/* loop on Nj^k's, local exp moments */
      for(n = 0, rhoN = rho; n <= order; rhoN *= (-rho), n++) {
	for(m = 0; m <= n; m++) { /* loop on On^m's, multipole moments */

	  /* generate a bar(N)j^k and dblbar(N)j^k entry */
	  rhoFac = rhoJ*rhoN;	/* divisor to give (-1)^n/rho^(j+n+1) factor */
	  if(k == 0) {		/* use abbreviated formulae in this case */

	    /* generate only bar(N)j^0 entry (dblbar(N)j^0 = 0 always) */
	    if(m != 0) {
	      temp1 = tleg[CINDEX(j+n, m)]*facFrA[j+n-m][n+m];
	      mat[CINDEX(j, 0)][CINDEX(n, m)] += temp1*cosB(m)/rhoFac;
	      mat[CINDEX(j, 0)][SINDEX(n, m, ct)] += temp1*sinB(m)/rhoFac;
	    }
	    else mat[CINDEX(j, 0)][CINDEX(n, 0)] 
		+= tleg[CINDEX(j+n, 0)]*facFrA[j+n][n]/rhoFac;
	  }
	  else {
	    temp1 = tleg[CINDEX(j+n, abs(m-k))]
		*facFrA[j+n-abs(m-k)][n+m]*iPwr(abs(k-m)-k-m);
	    temp2 = tleg[CINDEX(j+n, m+k)]*facFrA[j+n-m-k][n+m];
	    temp3 = tleg[CINDEX(j+n, k)]*facFrA[j+n-k][n]*2;

	    /* generate bar(N)j^k entry */
	    if(m != 0) {
	      mat[CINDEX(j, k)][CINDEX(n, m)] 
		  += (temp1*cosB(m-k)+temp2*cosB(m+k))/rhoFac;
	      mat[CINDEX(j, k)][SINDEX(n, m, ct)] 
		  += (temp1*sinB(m-k)+temp2*sinB(m+k))/rhoFac;
	    }
	    else mat[CINDEX(j, k)][CINDEX(n, 0)] += temp3*cosB(k)/rhoFac;

	    /* generate dblbar(N)j^k entry */
	    if(m != 0) {
	      mat[SINDEX(j, k, ct)][CINDEX(n, m)] 
		  += (-temp1*sinB(m-k)+temp2*sinB(m+k))/rhoFac;
	      mat[SINDEX(j, k, ct)][SINDEX(n, m, ct)] 
		  += (temp1*cosB(m-k)-temp2*cosB(m+k))/rhoFac;
	    }
	    else mat[SINDEX(j, k, ct)][CINDEX(n, 0)] += temp3*sinB(k)/rhoFac;
	  }
	}
      }
    }
  }

#if DISM2L == ON
  dispM2L(mat, x, y, z, xp, yp, zp, order);
#endif

  return(mat);
}

/* 
  Used only for true (Greengard) downward pass - similar to Multi2Local
*/
double **mulLocal2Local(x, y, z, xc, yc, zc, order)
int order;
double x, y, z, xc, yc, zc;	/* parent and child cube centers */
{
  int i, j, k, n, m;
  int terms = multerms(order);	/* the number of non-zero moments */
  int ct = costerms(order);	/* the number of non-zero cos (bar) moments */
  double **mat;			/* the transformation matrix */
  double rho, cosA, beta;	/* spher. position of multi rel to local */
  double rhoJ, rhoN;		/* rho^j and (-1)^n*rho^(n+1) in main loop */
  double rhoFac;		/* = rhoJ*rhoN intermediate storage */
  double temp1, temp2, temp3;
  double iPwr(), sinB(), cosB();
  extern double *tleg, *Ir, *Irn, *phi, *Mphi; /* external temporary storage */

  /* allocate the local to local transformation matrix */
  CALLOC(mat, terms, double*, ON, AL2L);
  for(i = 0; i < terms; i++)
      CALLOC(mat[i], terms, double, ON, AL2L);

  /* find relative spherical coordinates */
  xyz2sphere(x, y, z, xc, yc, zc, &rho, &cosA, &beta);

  /* generate legendre function evaluations */
  evalLegendre(cosA, tleg, 2*order); /* local->local needs 2x legendres */

  /* generate sin[(m+-k)beta] and cos[(m+-k)beta] look up arrays */
  /*  other lookup arrays generated in mulMultiAlloc() */
  evalSinCos(beta, order);

  /* generate local to local transformation matrix; uses NB12 pg36Y */
  /*  rhoFac factor divides could be reduced to once per loop */
  for(j = 0, rhoJ = 1.0; j <= order; rhoJ *= (-rho), j++) {
    for(k = 0; k <= j; k++) {	/* loop on Nj^k's, local exp moments */
      for(n = j, rhoN = rhoJ; n <= order; rhoN *= (-rho), n++) {
	for(m = 0; m <= n; m++) { /* loop on On^m's, old local moments */

	  /* generate a bar(N)j^k and dblbar(N)j^k entry */
	  rhoFac = rhoN/rhoJ;	/* divide to give (-rho)^(n-j) factor */
	  if(k == 0 && n-j >= m) {  /* use abbreviated formulae in this case */

	    /* generate only bar(N)j^0 entry (dblbar(N)j^0 = 0 always) */
	    if(m != 0) {
	      temp1 = tleg[CINDEX(n-j, m)]*facFrA[0][n-j+m]*rhoFac;
	      mat[CINDEX(j, 0)][CINDEX(n, m)] += temp1*cosB(m);
	      mat[CINDEX(j, 0)][SINDEX(n, m, ct)] += temp1*sinB(m);
	    }
	    else mat[CINDEX(j, 0)][CINDEX(n, 0)] += tleg[CINDEX(n-j, 0)]
		    *facFrA[0][n-j]*rhoFac;
	  }
	  else {
	    if(n-j >= abs(m-k)) temp1 = tleg[CINDEX(n-j, abs(m-k))]
		*facFrA[0][n-j+abs(m-k)]*iPwr(m-k-abs(m-k))*rhoFac;
	    if(n-j >= m+k) temp2 = tleg[CINDEX(n-j, m+k)]
		*facFrA[0][n-j+m+k]*iPwr(2*k)*rhoFac;
	    if(n-j >= k) temp3 = 2*tleg[CINDEX(n-j, k)]
		*facFrA[0][n-j+k]*iPwr(2*k)*rhoFac;

	    /* generate bar(N)j^k entry */
	    if(m != 0) {
	      if(n-j >= abs(m-k)) {
		mat[CINDEX(j, k)][CINDEX(n, m)] += temp1*cosB(m-k);
		mat[CINDEX(j, k)][SINDEX(n, m, ct)] += temp1*sinB(m-k);
	      }
	      if(n-j >= m+k) {
		mat[CINDEX(j, k)][CINDEX(n, m)] += temp2*cosB(m+k);
		mat[CINDEX(j, k)][SINDEX(n, m, ct)] += temp2*sinB(m+k);
	      }
	    }
	    else if(n-j >= k) mat[CINDEX(j, k)][CINDEX(n, 0)] += temp3*cosB(k);

	    /* generate dblbar(N)j^k entry */
	    if(m != 0) {
	      if(n-j >= abs(m-k)) {
		mat[SINDEX(j, k, ct)][CINDEX(n, m)] += (-temp1*sinB(m-k));
		mat[SINDEX(j, k, ct)][SINDEX(n, m, ct)] += temp1*cosB(m-k);
	      }
	      if(n-j >= m+k) {
		mat[SINDEX(j, k, ct)][CINDEX(n, m)] += (-temp2*sinB(m+k));
		mat[SINDEX(j, k, ct)][SINDEX(n, m, ct)] += temp2*cosB(m+k);
	      }
	    }
	    else if(n-j >= k) 
		mat[SINDEX(j, k, ct)][CINDEX(n, 0)] += temp3*sinB(k);
	  }
	}
      }
    }
  }

#if DISL2L == ON
  dispL2L(mat, x, y, z, xc, yc, zc, order);
#endif

  return(mat);
}

/*
  sets up xformation for distant cube charges to local expansion
  form almost identical to mulQ2Multi - follows NB12 pg 32 w/m,n replacing k,j
  OPTIMIZATIONS INVOLVING is_dummy HAVE NOT BEEN COMPLETELY IMPLEMENTED
*/
double **mulQ2Local(chgs, numchgs, is_dummy, x, y, z, order)
double x, y, z;
charge **chgs;
int numchgs, order, *is_dummy;
{
  int i, j, k, kold, n, m, start;
  int cterms = costerms(order), terms = multerms(order);
  double **mat, temp, fact();
  double cosA;			/* cosine of elevation coordinate */
  extern double *Rhon, *Rho, *Betam, *Beta, *tleg;

  /* Allocate the matrix. */
  CALLOC(mat, terms, double*, ON, AQ2L);
  for(i=0; i < terms; i++) 
      CALLOC(mat[i], numchgs, double, ON, AQ2L);

  /* get Legendre function evaluations, one set for each charge */
  /*  also get charge coordinates, set up for subsequent evals */
  for(j = 0; j < numchgs; j++) { /* for each charge */

    /* get cosA for eval; save rho, beta in rho^n and cos/sin(m*beta) arrays */
    xyz2sphere(chgs[j]->x, chgs[j]->y, chgs[j]->z,
	       x, y, z, &(Rho[j]), &cosA, &(Beta[j]));
    Rhon[j] = Rho[j]; /* init powers of rho_i's */
    Betam[j] = Beta[j];		/* init multiples of beta */
    evalLegendre(cosA, tleg, order);	/* write moments to temporary array */
    /* write a column of the matrix with each set of legendre evaluations */
    for(i = 0; i < cterms; i++) mat[i][j] = tleg[i]; /* copy for cos terms */
  }

#if DALQ2L == ON
  fprintf(stdout,
	  "\nQ2L MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, terms, numchgs);
#endif

  /* add the rho^n+1 factors to the cos matrix entries. */
  for(n = 0; n <= order; n++) {	/* loop on rows of matrix */
    for(m = 0; m <= n; m++) {
      for(j = 0; j < numchgs; j++) 
	  mat[CINDEX(n, m)][j] /= Rhon[j]; /* divide by factor */
    }
    for(j = 0; j < numchgs; j++) Rhon[j] *= Rho[j];	/* rho^n -> rho^n+1 */
  }

#if DALQ2L == ON
  fprintf(stdout,"    AFTER ADDITION OF (1/RHO)^N+1 FACTORS\n");
  dumpMat(mat, terms, numchgs);
#endif

  /* copy result to lower (sine) part of matrix */
  for(n = 1; n <= order; n++) { /* loop on rows of matrix */
    for(m = 1; m <= n; m++) {
      for(j = 0; j < numchgs; j++) { /* copy a row */
	mat[SINDEX(n, m, cterms)][j] = mat[CINDEX(n, m)][j];
      }
    }
  }

#if DALQ2L == ON
  fprintf(stdout,"    AFTER COPYING SINE (LOWER) HALF\n");
  dumpMat(mat, terms, numchgs);
#endif

  /* add factors of cos(m*beta) and sin(m*beta) to matrix entries */
  for(m = 0; m <= order; m++) {	/* lp on m in Mn^m */
    for(n = m; n <= order; n++) { /* loop over rows with same m */
      for(j = 0; j < numchgs; j++) { /* add factors to a row */
	if(m == 0)  mat[CINDEX(n, m)][j] *= fact(n); /* j! part of bar(N)j^0 */
	else {			/* for Nj^k, k != 0 */
	  temp = 2.0*fact(n-m);    	/* find the factorial for moment */
	  mat[CINDEX(n, m)][j] *= (temp*cos(Betam[j]));   /* note mul by 2 */
	  mat[SINDEX(n, m, cterms)][j] *= (temp*sin(Betam[j]));
	}
      }
    }
    if(m > 0) {
      for(j = 0; j < numchgs; j++) Betam[j] += Beta[j];/* (m-1)*beta->m*beta */
    }
  }

  /* THIS IS NOT VERY GOOD: zero out columns corresponding to dummy panels */
  for(j = 0; j < numchgs; j++) {
    if(is_dummy[j]) {
      for(i = 0; i < terms; i++) {
	mat[i][j] = 0.0;
      }
    }
  }

#if DISQ2L == ON
   dispQ2L(mat, chgs, numchgs, x, y, z, order);
#endif

  return(mat);
}

/*
  builds local expansion evaluation matrix; not used for fake dwnwd pass
  follows NB10 equation marked circle(2A) except roles of j,k and n,m switched
  very similar to mulMulti2P()
*/
double **mulLocal2P(x, y, z, chgs, numchgs, order)
double x, y, z;
charge **chgs;
int numchgs, order;
{
  double **mat;
  double cosTh;			/* cosine of elevation coordinate */
  double fact();
  extern double *Irn, *Mphi, *phi, *Ir;
  int i, j, k, m, n, kold, start;
  int cterms = costerms(order), terms = multerms(order);

  CALLOC(mat, numchgs, double*, ON, AL2P);
  for(i = 0; i < numchgs; i++) 
      CALLOC(mat[i], terms, double, ON, AL2P);

  /* get Legendre function evaluations, one set for each charge */
  /*   also get charge coordinates to set up rest of matrix */
  for(i = 0; i < numchgs; i++) { /* for each charge, do a legendre eval set */
    xyz2sphere(chgs[i]->x, chgs[i]->y, chgs[i]->z,
	       x, y, z, &(Ir[i]), &cosTh, &(phi[i]));
    Irn[i] = 1.0; /* initialize r^n vec. */
    Mphi[i] = phi[i];		/* initialize m*phi vector */
    evalLegendre(cosTh, mat[i], order);	/* wr moms to 1st (cos) half of row */
  }

#if DALL2P == ON
  fprintf(stdout,
	  "\nL2P MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* add the r^n factors to the left (cos(m*phi)) half of the matrix */
  for(j = 0, k = kold = 1; j < cterms; j++) { /* loop on columns of matrix */
    for(i = 0; i < numchgs; i++) mat[i][j] *= Irn[i]; /* multiply by r^n */
    k -= 1;
    if(k == 0) {		/* so that n changes as appropriate */
      kold = k = kold + 1;
      for(i = 0; i < numchgs; i++) Irn[i] *= Ir[i]; /* r^n -> r^n+1 */
    }
  }

#if DALL2P == ON
  fprintf(stdout,
	  "    AFTER ADDITION OF R^N FACTORS\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* add the factorial fraction factors to the left (cos(m*phi)) part of mat */
  for(n = 0; n <= order; n++) {
    for(m = 0; m <= n; m++) {
      for(i = 0; i < numchgs; i++) mat[i][CINDEX(n, m)] /= fact(n+m);
    }
  }

#if DALL2P == ON
  fprintf(stdout,
	  "    AFTER ADDITION OF FACTORIAL FRACTION FACTORS\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* copy left half of matrix to right half for sin(m*phi) terms */
  for(i = 0; i < numchgs; i++) { /* loop on rows of matrix */
    for(n = 1; n <= order; n++) { 
      for(m = 1; m <= n; m++) {	/* copy a row */
	mat[i][SINDEX(n, m, cterms)] = mat[i][CINDEX(n, m)];
      }
    }
  }

#if DALL2P == ON
  fprintf(stdout,
	  "    AFTER COPYING SINE (RIGHT) HALF\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* add factors of cos(m*phi) and sin(m*phi) to left and right halves resp. */
  for(m = 1; m <= order; m++) {	/* lp on m in Mn^m (no m=0 since cos(0)=1) */
    for(n = m; n <= order; n++) { /* loop over cols with same m */
      for(i = 0; i < numchgs; i++) { /* add factors to a column */
	mat[i][CINDEX(n, m)] *= cos(Mphi[i]);
	mat[i][SINDEX(n, m, cterms)] *= sin(Mphi[i]);
      }
    }
    for(i = 0; i < numchgs; i++) Mphi[i] += phi[i]; /* (m-1)*phi->m*phi */
  }

#if DISL2P == ON
  dispL2P(mat, x, y, z, chgs, numchgs, order);
#endif

  return(mat);
}
