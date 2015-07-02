/* # ***** sort to /src/main
   # ***** */
#include "mulGlobal.h"

/*
  Globals used for temporary storage.
*/
double *Irn, *Mphi;		/* (1/r)^n+1, m*phi vect's */
double *Ir, *phi;		/* 1/r and phi arrays, used to update above */
double *Rho, *Rhon;		/* rho and rho^n array */
double *Beta, *Betam;		/* beta and beta*m array */
double *tleg;		/* Temporary Legendre storage. */
double **factFac;		/* factorial factor array: (n-m+1)...(n+m) */

/* 
   Used various places.  Returns number of coefficients in the multipole 
   expansion. 
*/
int multerms(order)
int order;
{
  return(costerms(order) + sinterms(order));
}

/*
  returns number of cos(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int costerms(order)
int order;
{
  return(((order + 1) * (order + 2)) / 2);
}

/*
  returns number of sin(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int sinterms(order)
int order;
{
  return((((order + 1) * (order + 2)) / 2) - (order+1));
}


/*
  takes two sets of cartesian absolute coordinates; finds rel. spherical coor.
*/
void xyz2sphere(x, y, z, x0, y0, z0, rho, cosA, beta)
double x, y, z, x0, y0, z0, *rho, *cosA, *beta;
{
  /* get relative coordinates */
  x -= x0;			/* "0" coordinates play the role of origin */
  y -= y0;
  z -= z0;
  /* get spherical coordinates */
  *rho = sqrt(x*x + y*y + z*z);

  if(*rho == 0.0) *cosA = 1.0;
  else *cosA = z/(*rho);

  if(x == 0.0 && y == 0.0) *beta = 0.0;
  else *beta = atan2(y, x);

}

/*
  gives the linear index into vector from n and m used by all routines dealing
  with moments (cosine parts) and Leg. function evals, e.g. Mn^m and Pn^m(cosA)
  used for all cases except for the sine part of arrays (use sindex() instead)
  assumed entry order: (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0)...
*/
/* REPLACED BY MACRO CINDEX(N, M) 24July91 */

#if TRUE == FALSE

int index(n, m)
int n, m;
{
  if(m > n) {
    fprintf(stderr, "index: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0) {
    fprintf(stderr, "index: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
  return(m + (n*(n+1))/2);
}

/*
  gives the linear index into vector from n and m used by all routines dealing
  with moments (sine parts), e.g. Mn^m 
  assumes an array with all m = 0 (Mn^0) entries omitted to save space
  assumed entry order: (n,m) = (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1)...
*/
/* REPLACED BY MACRO SINDEX(N, M, CTERMS) 24July91 */
int sindex(n, m, cterms)
int n, m, cterms;		/* cterms is costerms(order) */
{
  if(m > n) {
    fprintf(stderr, "sindex: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0) {
    fprintf(stderr, "sindex: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
  if(m == 0) {
    fprintf(stderr, "sindex: attempt to index M%d^0\n", n);
    exit(0);
  }
  return(cterms + m + (n*(n+1))/2 - (n+1));
}

#endif				/* end of index(), sindex() comment out */

/*
  returns i = sqrt(-1) to the power of the argument
*/
double iPwr(e)
int e;				/* exponent, computes i^e */
{
  if(e == 0) return(1.0);
  if(e % 2 != 0) {
    fprintf(stderr, "iPwr: odd exponent %d\n", e);
    exit(0);
  }
  else {
    e = e/2;			/* get power of negative 1 */
    if(e % 2 == 0) return(1.0);
    else return(-1.0);
  }
}

/*
  returns factorial of the argument (x!)
*/
double fact(x)
int x;
{
  double ret = 1.0;
  if(x == 0 || x == 1) return(1.0);
  else if(x < 0) {
    fprintf(stderr, "fact: attempt to take factorial of neg no. %d\n", x);
    exit(0);
  }
  else {
    while(x > 1) {
      ret *= x;
      x--;
    }
    return(ret);
  }
}

/*
  produces factorial factor array for mulMulti2P
*/
void evalFactFac(array, order)
int order;
double **array;
{
  int n, m;			/* array[n][m] = (m+n)!/(n-m)! */

  /* do first column of lower triangular part - always 1's */
  for(n = 0; n < order+1; n++) array[n][0] = 1.0;

  /* do remaining columns of lower triangular part */
  /*  use (n-m)!/(n+m)! = 1/(n-m+1)...(n+m) since m \leq n */
  /*  (array entry is divided into number to effect mul by (n-m)!/(n+m)!) */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      array[n][m] = (n-(m-1)) * array[n][m-1] * (n+m);
    }
  }

#if DISFAF == ON
  fprintf(stdout, "FACTORIAL FACTOR ARRAY:\n");
  dumpMat(array, order+1, order+1);
#endif

}

/*
  Allocates space for temporary vectors.
*/
void mulMultiAlloc(maxchgs, order, depth)
int maxchgs, order, depth;
{
  int x;

#if DISSYN == ON
  extern int *multicnt, *localcnt, *evalcnt;
#endif
#if DMTCNT == ON
  extern int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt;
  extern int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt, **Q2PDcnt;
#endif
  extern double *sinmkB, *cosmkB, **facFrA;

  if(maxchgs > 0) {
    CALLOC(Rho, maxchgs, double, ON, AMSC); /* rho array */
    CALLOC(Rhon, maxchgs, double, ON, AMSC); /* rho^n array */
    CALLOC(Beta, maxchgs, double, ON, AMSC); /* beta array */
    CALLOC(Betam, maxchgs, double, ON, AMSC); /* beta*m array */
    CALLOC(Irn, maxchgs, double, ON, AMSC);	/* (1/r)^n+1 vector */
    CALLOC(Ir, maxchgs, double, ON, AMSC); /* 1/r vector */
    CALLOC(Mphi, maxchgs, double, ON, AMSC); /* m*phi vector */
    CALLOC(phi, maxchgs, double, ON, AMSC); /* phi vector */
  }
  CALLOC(tleg, costerms(2*order), double, ON, AMSC);
	       	/* temp legendre storage (2*order needed for local exp) */
  CALLOC(factFac, order+1, double*, ON, AMSC);
  for(x = 0; x < order+1; x++) {
    CALLOC(factFac[x], order+1, double, ON, AMSC);
  }
  evalFactFac(factFac, order);	/* get factorial factors for mulMulti2P */

#if DISSYN == ON
  /* for counts of local/multipole expansions and eval mat builds by level */
  CALLOC(localcnt, depth+1, int, ON, AMSC);
  CALLOC(multicnt, depth+1, int, ON, AMSC);
  CALLOC(evalcnt, depth+1, int, ON, AMSC);
#endif

#if DMTCNT == ON
  /* for counts of transformation matrices by level  */
  CALLOC(Q2Mcnt, depth+1, int*, ON, AMSC);
  CALLOC(Q2Lcnt, depth+1, int*, ON, AMSC);
  CALLOC(Q2Pcnt, depth+1, int*, ON, AMSC);
  CALLOC(L2Lcnt, depth+1, int*, ON, AMSC);
  CALLOC(M2Mcnt, depth+1, int*, ON, AMSC); 
  CALLOC(M2Lcnt, depth+1, int*, ON, AMSC); 
  CALLOC(M2Pcnt, depth+1, int*, ON, AMSC); 
  CALLOC(L2Pcnt, depth+1, int*, ON, AMSC);
  CALLOC(Q2PDcnt, depth+1, int*, ON, AMSC);
  for(x = 0; x < depth+1; x++) {
    CALLOC(Q2Mcnt[x], depth+1, int, ON, AMSC);
    CALLOC(Q2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(Q2Pcnt[x], depth+1, int, ON, AMSC);
    CALLOC(L2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Mcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Pcnt[x], depth+1, int, ON, AMSC);
    CALLOC(L2Pcnt[x], depth+1, int, ON, AMSC);
    CALLOC(Q2PDcnt[x], depth+1, int, ON, AMSC);
  }
#endif

  /* from here down could be switched out when the fake dwnwd pass is used */
  CALLOC(facFrA, 2*order+1, double*, ON, AMSC);
  for(x = 0; x < 2*order+1; x++) {
    CALLOC(facFrA[x], 2*order+1, double, ON, AMSC);
  }
  /* generate table of factorial fraction evaluations (for M2L and L2L) */
  evalFacFra(facFrA, order);
  CALLOC(sinmkB, 2*order+1, double, ON, AMSC); /* sin[(m+-k)beta] */
  CALLOC(cosmkB, 2*order+1, double, ON, AMSC); /* cos[(m+-k)beta] */
  cosmkB[0] = 1.0;		/* look up arrays used for local exp */
  /* generate array of sqrt((n+m)!(n-m)!)'s for L2L
  evalSqrtFac(sqrtFac, factFac, order); */
}

/*
  returns an vector of Legendre function evaluations of form Pn^m(cosA)
  n and m have maximum value order
  vector entries correspond to (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1)...
*/
void evalLegendre(cosA, vector, order)
double cosA, *vector;
int order;
{
  int x;
  int n, m;			/* as in Pn^m, both <= order */
  double sinMA;			/* becomes sin^m(alpha) in higher order P's */
  double fact;			/* factorial factor */

  /* do evaluations of first four functions separately w/o recursions */
  vector[CINDEX(0, 0)] = 1.0;	/* P0^0 */
  if(order > 0) {
    vector[CINDEX(1, 0)] = cosA;	/* P1^0 */
    vector[CINDEX(1, 1)] = sinMA = -sqrt(1-cosA*cosA); /* P1^1 = -sin(alpha) */
  }
  if(order > 1) vector[CINDEX(2, 1)] = 3*sinMA*cosA; /* P2^1 = -3sin()cos() */

  /* generate remaining evaluations by recursion on lower triangular array */
  fact = 1.0;
  for(m = 0; m < order+1; m++) {
    if(m != 0 && m != 1) {	/* set up first two evaluations in row */
      fact *= (2*m - 1); /* (2(m-1)-1)!! -> (2m-1)!! */
      /* use recursion on m */
      if(vector[CINDEX(1, 1)] == 0.0) {
	vector[CINDEX(m, m)] = 0.0;
	if(m != order) vector[CINDEX(m+1, m)] = 0.0;	/* if not last row */
      }
      else {
	cosA = vector[CINDEX(1,0)]/vector[CINDEX(1,1)]; /* cosA= -cot(theta) */
	sinMA *= vector[CINDEX(1,1)]; /*(-sin(alpha))^(m-1)->(-sin(alpha))^m*/
	vector[CINDEX(m, m)] = fact * sinMA;
	if(m != order) {		/* do if not on last row */
	  vector[CINDEX(m+1, m)] = vector[CINDEX(1, 0)]*(2*m+1)
	      *vector[CINDEX(m, m)];
	}
      }
    }
    for(x = 2; x < order-m+1; x++) { /* generate row of evals recursively */
      vector[CINDEX(x+m, m)] = 
	  ((2*(x+m)-1)*vector[CINDEX(1, 0)]*vector[CINDEX(x+m-1, m)]
	   - (x + 2*m - 1)*vector[CINDEX(x+m-2, m)])/x;
    }
  }
}
    
  
/* 
  Returns a matrix which gives a cube's multipole expansion when *'d by chg vec
  OPTIMIZATIONS USING is_dummy HAVE NOT BEEN COMPLETELY IMPLEMENTED
*/
double **mulQ2Multi(chgs, is_dummy, numchgs, x, y, z, order)
double x, y, z; 
charge **chgs;
int numchgs, order, *is_dummy;
{
  double **mat;
  double cosA;			/* cosine of elevation coordinate */
  int i, j, k, kold, n, m, start;
  int cterms = costerms(order), terms = multerms(order);

  /* Allocate the matrix. */
  CALLOC(mat, terms, double*, ON, AQ2M);
  for(i=0; i < terms; i++) 
      CALLOC(mat[i], numchgs, double, ON, AQ2M);

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

#if DALQ2M == ON
  fprintf(stdout,
	  "\nQ2M MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, terms, numchgs);
#endif

  /* some of this can be avoided using is_dummy to skip unused columns */
  /* add the rho^n factors to the cos matrix entries. */
  for(i = 1, k = kold = 2; i < cterms; i++) { /* loop on rows of matrix */
    for(j = 0; j < numchgs; j++) mat[i][j] *= Rhon[j]; /* mul in factor */
    k -= 1;
    if(k == 0) {		/* so that effective n varys appropriately */
      kold = k = kold + 1;
      for(j = 0; j < numchgs; j++) Rhon[j] *= Rho[j];	/* r^n-1 -> r^n */
    }
  }

#if DALQ2M == ON
  fprintf(stdout,"    AFTER ADDITION OF RHO^N FACTORS\n");
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

#if DALQ2M == ON
  fprintf(stdout,"    AFTER COPYING SINE (LOWER) HALF\n");
  dumpMat(mat, terms, numchgs);
#endif

  /* add factors of cos(m*beta) and sin(m*beta) to matrix entries */
  for(m = 1; m <= order; m++) {	/* lp on m in Mn^m (no m=0 since cos(0)=1) */
    for(n = m; n <= order; n++) { /* loop over rows with same m */
      for(j = 0; j < numchgs; j++) { /* add factors to a row */
	mat[CINDEX(n, m)][j] *= (2.0*cos(Betam[j]));   /* note factors of 2 */
	mat[SINDEX(n, m, cterms)][j] *= (2.0*sin(Betam[j]));
      }
    }
    for(j = 0; j < numchgs; j++) Betam[j] += Beta[j]; /* (m-1)*beta->m*beta */
  }

  /* THIS IS NOT VERY GOOD: zero out columns corresponding to dummy panels */
  for(j = 0; j < numchgs; j++) {
    if(is_dummy[j]) {
      for(i = 0; i < terms; i++) {
	mat[i][j] = 0.0;
      }
    }
  }

#if DISQ2M == ON
  dispQ2M(mat, chgs, numchgs, x, y, z, order);
#endif

  return(mat);
}

double **mulMulti2Multi(x, y, z, xp, yp, zp, order)
double x, y, z, xp, yp, zp;	/* cube center, parent cube center */
int order;
{
  double **mat, rho, rhoPwr, cosA, beta, mBeta, temp1, temp2; 
  double iPwr(), fact();
  int r, j, k, m, n, c;
  int cterms = costerms(order), sterms = sinterms(order);
  int terms = cterms + sterms;

  /* Allocate the matrix (terms x terms ) */
  CALLOC(mat, terms, double*, ON, AM2M);
  for(r=0; r < terms; r++) CALLOC(mat[r], terms, double, ON, AM2M);
  for(r = 0; r < terms; r++)
      for(c = 0; c < terms; c++) mat[r][c] = 0.0;

  /* get relative distance in spherical coordinates */
  xyz2sphere(x, y, z, xp, yp, zp, &rho, &cosA, &beta);

  /* get the requisite Legendre function evaluations */
  evalLegendre(cosA, tleg, order);

  /* for each new moment (Nj^k) stuff the appropriate matrix entries */
  /* done completely brute force, one term at a time; uses exp in nb 12, p29 */
  for(j = 0; j <= order; j++) {
    for(k = 0; k <= j; k++) {
      for(n = 0, rhoPwr = 1.0; n <= j; n++, rhoPwr *= rho) {
	for(m = 0, mBeta = 0.0; m <= n; m++, mBeta += beta) {

	  if(k == 0) {		/* figure terms for Nj^0, ie k = 0 */
	    if(m <= j-n) {	/* if O moments are nonzero */
	      temp1 = fact(j)*rhoPwr*iPwr(2*m)*tleg[CINDEX(n, m)];
	      temp1 /= (fact(j-n+m)*fact(n+m));
	      mat[CINDEX(j, k)][CINDEX(j-n, m)] += temp1*cos(mBeta);
	      if(m != 0) {		/* if sin term is non-zero */
		mat[CINDEX(j, k)][SINDEX(j-n, m, cterms)] += temp1*sin(mBeta);
	      }
	    }
	  }
	  else {		/* figure terms for Nj^k, k != 0 */
	    temp1 = fact(j+k)*rhoPwr*tleg[CINDEX(n, m)]/fact(n+m);
	    temp2 = temp1*iPwr(2*m)/fact(j-n+k+m);
	    temp1 = temp1*iPwr(k-m-abs(k-m))/fact(j-n+abs(k-m));

	    /* write the cos(kPhi) coeff, bar(N)j^k */
	    if(m != 0) {
	      if(k-m < 0 && abs(k-m) <= j-n) {	/* use conjugates here */
		mat[CINDEX(j, k)][CINDEX(j-n, m-k)] += temp1*cos(mBeta);
		mat[CINDEX(j, k)][SINDEX(j-n, m-k,cterms)] += temp1*sin(mBeta);
	      }
	      else if(k-m == 0) {	/* double to compensate for 2Re sub. */
		mat[CINDEX(j, k)][CINDEX(j-n, k-m)] += 2*temp1*cos(mBeta);
		/* sin term is always zero */
	      }
	      else if(k-m > 0 && k-m <= j-n) {
		mat[CINDEX(j, k)][CINDEX(j-n, k-m)] += temp1*cos(mBeta);
		mat[CINDEX(j, k)][SINDEX(j-n, k-m,cterms)] -= temp1*sin(mBeta);
	      }
	      if(k+m <= j-n) {
		mat[CINDEX(j, k)][CINDEX(j-n, k+m)] += temp2*cos(mBeta);
		mat[CINDEX(j, k)][SINDEX(j-n, k+m,cterms)] += temp2*sin(mBeta);
	      }
	    }			/* do if m = 0 and O moments not zero */
	    else if(k <= j-n) mat[CINDEX(j, k)][CINDEX(j-n, k)] += temp2;

	    /* write the sin(kPhi) coeff, dblbar(N)j^k, if it is non-zero */
	    if(m != 0) {
	      if(k-m < 0 && abs(k-m) <= j-n) {	/* use conjugates here */
		mat[SINDEX(j, k,cterms)][CINDEX(j-n, m-k)] += temp1*sin(mBeta);
		mat[SINDEX(j, k, cterms)][SINDEX(j-n, m-k, cterms)] 
		    -= temp1*cos(mBeta);
	      }
	      else if(k-m == 0) {/* double to compensate for 2Re sub */
		mat[SINDEX(j, k, cterms)][CINDEX(j-n, k-m)] 
		    += 2*temp1*sin(mBeta);
		/* sine term is always zero */
	      }
	      else if(k-m > 0 && k-m <= j-n) {
		mat[SINDEX(j, k,cterms)][CINDEX(j-n, k-m)] += temp1*sin(mBeta);
		mat[SINDEX(j, k, cterms)][SINDEX(j-n, k-m, cterms)] 
		    += temp1*cos(mBeta);
	      }
	      if(k+m <= j-n) {
		mat[SINDEX(j, k,cterms)][CINDEX(j-n, k+m)] -= temp2*sin(mBeta);
		mat[SINDEX(j, k, cterms)][SINDEX(j-n, k+m, cterms)] 
		    += temp2*cos(mBeta);
	      }
	    }			/* do if m = 0 and moments not zero */
	    else if(k <= j-n) mat[SINDEX(j, k,cterms)][SINDEX(j-n, k, cterms)] 
		+= temp2;
	  }
	}
      }
    }
  }
#if DISM2M == ON
  dispM2M(mat, x, y, z, xp, yp, zp, order);
#endif
  return(mat);
}

/* 
  builds multipole evaluation matrix; used only for fake downward pass 
*/
double **mulMulti2P(x, y, z, chgs, numchgs, order)
double x, y, z;			/* multipole expansion origin */
charge **chgs;
int numchgs, order;
{
  double **mat;
  double cosTh;			/* cosine of elevation coordinate */
  double factorial;		/* 1/factorial = (n-m)!/(n+m)! */
  int i, j, k, m, n, kold, start;
  int cterms = costerms(order), sterms = sinterms(order);
  int terms = cterms + sterms;

  CALLOC(mat, numchgs, double*, ON, AM2P);
  for(i=0; i < numchgs; i++) 
      CALLOC(mat[i], terms, double, ON, AM2P);

  /* get Legendre function evaluations, one set for each charge */
  /*   also get charge coordinates to set up rest of matrix */
  for(i = 0; i < numchgs; i++) { /* for each charge, do a legendre eval set */
    xyz2sphere(chgs[i]->x, chgs[i]->y, chgs[i]->z,
	       x, y, z, &(Ir[i]), &cosTh, &(phi[i]));

    Irn[i] = Ir[i]; /* initialize (1/r)^n+1 vec. */
    Mphi[i] = phi[i];		/* initialize m*phi vector */

    evalLegendre(cosTh, mat[i], order);	/* wr moms to 1st (cos) half of row */

  }

#if DALM2P == ON
  fprintf(stdout,
	  "\nM2P MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* add the (1/r)^n+1 factors to the left (cos(m*phi)) half of the matrix */
  for(j = 0, k = kold = 1; j < cterms; j++) { /* loop on columns of matrix */
    for(i = 0; i < numchgs; i++) mat[i][j] /= Irn[i]; /* divide by r^n+1 */
    k -= 1;
    if(k == 0) {		/* so that n changes as appropriate */
      kold = k = kold + 1;
      for(i = 0; i < numchgs; i++) Irn[i] *= Ir[i]; /* r^n -> r^n+1 */
    }
  }

#if DALM2P == ON
  fprintf(stdout,
	  "    AFTER ADDITION OF (1/R)^N+1 FACTORS\n");
  dumpMat(mat, numchgs, terms);
#endif

  /* add the factorial fraction factors to the left (cos(m*phi)) part of mat */
  /*  note that (n-m)!/(n+m)! = 1/(n-m+1)...(n+m) since m \leq n */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      for(i = 0; i < numchgs; i++) mat[i][CINDEX(n, m)] /= factFac[n][m];
    }
  }

#if DALM2P == ON
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

#if DALM2P == ON
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

#if DISM2P == ON
  dispM2P(mat, x, y, z, chgs, numchgs, order);
#endif

  return(mat);
}







