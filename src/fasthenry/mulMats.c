/* # ***** sort to /src/main
   # ***** */
#include "mulGlobal.h"

double **Q2P(), **Q2PDiag();
double **mulMulti2P(), **mulQ2Multi(), **mulMulti2Multi();
double **mulLocal2Local(), **mulLocal2P(), **mulQ2Local(), **mulMulti2Local();

int *localcnt, *multicnt, *evalcnt;	/* counts of builds done by level */

int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt; /* counts of xformation mats */
int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt, **Q2PDcnt;

/*
MulMatDirect creates the matrices for the piece of the problem that is done
directly exactly.
*/
mulMatDirect(sys)
ssystem *sys;
{
  cube *nextc, *nextnbr;
  int i, nummats, **temp;
  extern double lutime, dirtime;

#if DIRSOL == ON || EXPGCR == ON
  double **ludecomp();
  extern double *trimat, *sqrmat; /* flattened triangular, square matrices */
  extern int up_size, eval_size;
  extern int *real_index;	/* for map btwn condensed/expanded vectors */
#endif

  /* First count the number of matrices to be done directly. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    for(nummats=1, i=0; i < nextc->numnbrs; i++) {
      nextnbr = nextc->nbrs[i];
      ASSERT(nextnbr->upnumvects > 0);
      nummats++;
    }

  /* Allocate space for the vects and mats. */
    nextc->directnumvects = nummats;
    if(nummats > 0) {
      CALLOC(nextc->directq, nummats, double*, ON, AMSC);
      CALLOC(temp, nummats, int*, ON, AMSC);
      CALLOC(nextc->directnumeles, nummats, int, ON, AMSC);
      CALLOC(nextc->directmats, nummats, double**, ON, AMSC);
/*      CALLOC(nextc->precondmats, nummats, double**, ON, AMSC); */
    }

    /* initialize the pointer from this cube to its part of dummy vector
       - save the self part found in indexkid() */
    temp[0] = nextc->nbr_is_dummy[0];
    nextc->nbr_is_dummy = temp;
  }

/* Now place in the matrices. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    nextc->directq[0] = nextc->upvects[0];
    nextc->directnumeles[0] = nextc->upnumeles[0];

    /* starttimer; */
#if DIRSOL == ON || EXPGCR == ON
    if(nextc == sys->directlist) {
      if(eval_size < MAXSIZ) {
	fprintf(stderr, 
		"mulMatDirect: non-block direct methods not supported\n");
	exit(0);
	/* if this is going to work, need a special, condensing Q2P
	   as well as some way to use it in the framework of the GCR loop */
	nextc->directmats[0] = Q2P(nextc->chgs, eval_size,
				   nextc->nbr_is_dummy[0], nextc->chgs, 
				   eval_size, TRUE);
      }
      else blkQ2Pfull(sys->directlist, up_size, eval_size,
		      &trimat, &sqrmat, &real_index, sys->is_dummy);
    }
    else nextc->directmats[0] 
	= Q2PDiag(nextc->chgs, nextc->upnumeles[0], nextc->nbr_is_dummy[0],
		  TRUE);
#else
    nextc->directmats[0] 
	= Q2PDiag(nextc->chgs, nextc->upnumeles[0], nextc->nbr_is_dummy[0],
		  TRUE);
/*
    nextc->precondmats[0] 
	= Q2PDiag(nextc->chgs, nextc->upnumeles[0], nextc->nbr_is_dummy[0],
		  FALSE);
*/
    /*dumpMatCor(nextc->directmats[0], (double *)NULL, nextc->upnumeles[0]);*/

#endif
    /* stoptimer; */
    dirtime += dtime;

#if DSQ2PD == ON
    dumpQ2PDiag(nextc);
#endif

#if DMTCNT == ON
    Q2PDcnt[nextc->level][nextc->level]++;
#endif

#if DIRSOL == ON
    /* transform A into LU */
    if(eval_size > MAXSIZ) {
      blkLUdecomp(sqrmat, trimat, up_size);
    }
    else if(nextc == sys->directlist) {
      /* starttimer; */
      nextc->directlu = ludecomp(nextc->directmats[0], eval_size, TRUE);
      /* stoptimer; */
      lutime += dtime;
    }
#endif
    
    /* starttimer; */
    for(nummats=1, i=0; i < nextc->numnbrs; i++) {
      nextnbr = nextc->nbrs[i];
      ASSERT(nextnbr->upnumvects > 0);
      nextc->directq[nummats] = nextnbr->upvects[0];
      nextc->nbr_is_dummy[nummats] = nextnbr->nbr_is_dummy[0];
      nextc->directnumeles[nummats] = nextnbr->upnumeles[0];
      nextc->directmats[nummats] = Q2P(nextnbr->chgs, 
				       nextnbr->upnumeles[0], 
				       nextnbr->nbr_is_dummy[0],
				       nextc->chgs, nextc->upnumeles[0],
				       TRUE);
      nummats++;
/*
      nextc->precondmats[nummats++] = Q2P(nextnbr->chgs, 
					  nextnbr->upnumeles[0], 
					  nextnbr->nbr_is_dummy[0],
					  nextc->chgs, nextc->upnumeles[0],
					  FALSE);
*/
#if DMTCNT == ON
      Q2Pcnt[nextc->level][nextnbr->level]++;
#endif
    }
    /* stoptimer; */
    dirtime += dtime;
  }
}


/*
MulMatPrecond creates the preconditioner matrix
*/
bdmulMatPrecond(sys)
ssystem *sys;
{
  cube *nc, *kid, *kidnbr;
  double **mat, **nbrmat;
  int i, j, k, l, kidi;
  int kidsize, nbrsize, size, row, col, first, offset;
  double **ludecomp(), factor;
  charge *pc;
  surface *surf;

  printf("This Preconditioner is not used in FastHenry\n");
  exit(1);

  for(nc=sys->precondlist; nc != NULL; nc = nc->pnext) {

    /* find total number of charges in cube to dimension P. */
    for(size=0, i=0; i < nc->numkids; i++) {
      kid = nc->kids[i];
      if(kid != NULL) {
	ASSERT(kid->level == sys->depth);
	size += kid->directnumeles[0];  /* Equals number of charges. */
      }
    }

    /* allocate and zero a preconditioner matrix. */
    MALLOC(mat, size, double*, ON, AMSC);
    for(i = 0; i < size; i++) {
      MALLOC(mat[i], size, double, ON, AMSC);
    }
    for(i = 0; i < size; i++) {
      for(j = 0; j < size; j++) {
	mat[i][j] = 0.0;
      }
    }

    /* Chase through the kids to place in potential coeffs. */
    for(first = TRUE, row=0, kidi=0; kidi < nc->numkids; kidi++) {
      kid = nc->kids[kidi];
      if(kid != NULL) {
	/* Exploit the hierarchical charge numbering to get precond vector. */
	if(first == TRUE) {
	  first = FALSE;
	  nc->prevectq = kid->directq[0];
	  nc->prevectp = kid->eval;
	}
	/* Get the diagonal block of P^{-1}. */
	kidsize = kid->directnumeles[0];
	for(k = kidsize - 1; k >= 0; k--) {
	  for(l = kidsize - 1; l >= 0; l--) {
	    mat[row + k][row + l] = kid->directmats[0][k][l];
	  }
	}
	/* Get the off-diagonals of P^{-1}. */
	for(col = 0, i = 0; i < nc->numkids; i++) {
	  kidnbr = nc->kids[i];
	  if(kidnbr != NULL) {
	    if(kidnbr != kid) {
	      /* Chase thru list of nbrs to get matrix associated with this 
		 kidnbr.  Note, this is because the kid list and nbr matrix
		 list are in different orders, could be fixed. */
	      for(j = kid->numnbrs-1; j >= 0; j--) {
		if(kidnbr == kid->nbrs[j]) {
		  nbrmat = kid->directmats[j+1];
		  nbrsize = kidnbr->directnumeles[0];
		  for(k = kidsize - 1; k >= 0; k--) {
		    for(l = nbrsize - 1; l >= 0; l--) {
		      mat[row + k][col + l] = nbrmat[k][l];
		    }
		  }
		  break;
		}
	      }
	    }
	    col += kidnbr->directnumeles[0];
	  }
	}
	ASSERT(col == size);
	row += kidsize;
      }
    }    
    ASSERT(row == size);

    nc->precond = ludecomp(mat, size, FALSE);
    nc->presize = size;
  }
}
      

/* This near picks up only the hamming distance one cubes. */    
#define HNEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) + ABS((nbr)->k - (nk)) + ABS((nbr)->l - (nl))) <= 1)

/* This near picks up all 27 neighboring cubes. */
#define NEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) <= 1) && \
 (ABS((nbr)->k - (nk)) <= 1) && \
 (ABS((nbr)->l - (nl)) <= 1))

/* This near picks only the diagonal, for testing. */
#define DNEAR(nbr, nj, nk, nl) \
(((nbr)->j == (nj)) && \
 ((nbr)->k == (nk)) && \
 ((nbr)->l == (nl)) )
	
olmulMatPrecond(sys)
ssystem *sys;
{
  cube *nc, *nnbr, *nnnbr;
  double **mat, **nmat;
  int i, j, k, l, m;
  int maxsize, nsize, nnsize, nnnsize, *reorder;
  int nj, nk, nl, offset, noffset;
  int dindex, *nc_dummy, *nnbr_dummy, *nnnbr_dummy;
  static int *is_dummy;		/* local dummy flag vector, stays around */
  static int big_mat_size = 0;	/* size of previous mat */
  charge **nnnbr_pc, **nnbr_pc, **nc_pc, **mpc, *dp;
  surface *surf;
  double factor;

  printf("This Preconditioner is not used in FastHenry\n");
  exit(1);

/* Figure out the max number of elements in any set of near cubes. */
  for(maxsize=0, nc=sys->directlist; nc != NULL; nc = nc->dnext) {
    nsize = nc->directnumeles[0];
    nj = nc->j;
    nk = nc->k;
    nl = nc->l;
    for(i=0; i < nc->numnbrs; i++) {
      nnbr = nc->nbrs[i];
      if(NEAR(nnbr, nj, nk, nl)) nsize += nnbr->directnumeles[0];
    }
    maxsize = MAX(nsize, maxsize);
  }

/* Allocate a matrix big enough for any set of 7. */
#if JACDBG == ON
  printf("max direct size =%d\n", maxsize);
#endif
  MALLOC(reorder, maxsize,int, ON, AMSC);
  MALLOC(mat, maxsize, double*, ON, AMSC);
  for(i=0; i < maxsize; i++) {
    MALLOC(mat[i], maxsize, double, ON, AMSC);
  }

/* Now go fill-in a matrix. */
  for(maxsize=0, nc=sys->directlist; nc != NULL; nc = nc->dnext) {
    nsize = nc->directnumeles[0];
    nc_dummy = nc->nbr_is_dummy[0];
    nc_pc = nc->chgs;
#if CHKDUM == ON
    chkDummyList(nc_pc, nc_dummy, nsize);
#endif
    nj = nc->j;
    nk = nc->k;
    nl = nc->l;
    for(i = nsize - 1; i >= 0; i--) {
      if(nc_dummy[i]) continue;	/* dummy rows copied only in divided diff */
      if(nc_pc[i]->surf->type != DIELEC) {
	for(j = nsize - 1; j >= 0; j--) {
	  mat[i][j] = nc->directmats[0][i][j];
	}
      }
      else {
#if DPCOMP == ON
	fprintf(stdout, "Source mat, nc to nc\n");
	dumpMat(nc->directmats[0], nsize, nsize);
#endif
	find_flux_density_row(mat, nc->directmats[0], i, nsize, nsize, 0, 0,
			      nc_pc, nc_pc, nc_dummy, nc_dummy);
      }	
    }
    offset = nsize;
    for(k=0; k < nc->numnbrs; k++) { /* loop on neighbors of nc */
      nnbr = nc->nbrs[k];
      if(NEAR(nnbr, nj, nk, nl)) {
	nnsize = nc->directnumeles[k+1];
	nmat = nc->directmats[k+1];
	ASSERT(nc->directnumeles[k+1] == nnbr->directnumeles[0]);
	nnbr_dummy = nnbr->nbr_is_dummy[0];
	nnbr_pc = nnbr->chgs;
#if CHKDUM == ON
	chkDummyList(nnbr_pc, nnbr_dummy, nnsize);
#endif
	for(i = nsize - 1; i >= 0; i--) {
	  if(nc_dummy[i]) continue;
	  if(nc_pc[i]->surf->type != DIELEC) {
	    for(j = nnsize - 1; j >= 0; j--) {
	      mat[i][offset + j] = nmat[i][j];
	    }
	  }
	  else {
#if DPCOMP == ON
	    fprintf(stdout, "Source mat, nnbr to nc\n");
	    dumpMat(nmat, nsize, nnsize);
#endif
	    find_flux_density_row(mat, nmat, i, nnsize, nsize, 0, offset,
				  nc_pc, nnbr_pc, nc_dummy, nnbr_dummy);
	  }
	}
	/* Get the row of the big matrix associated with this nnbr. */
	for(noffset = 0, l = -1; l < nc->numnbrs; l++) { /* lp on nc's nbrs */
	  if(l < 0) nnnbr = nc;
	  else nnnbr = nc->nbrs[l];
	  if(NEAR(nnnbr, nj, nk, nl)) {  /* Note, near to nc!! */
	    if(nnbr == nnnbr) m = -1;
	    else { /* Find this nnnbr's position in nnbr's list */
	      for(m=0; m < nnbr->numnbrs; m++) {
		if(nnbr->nbrs[m] == nnnbr) break;
	      }
	      ASSERT(m < nnbr->numnbrs);
	    }
	    nnnsize = nnbr->directnumeles[m+1];
	    nmat = nnbr->directmats[m+1];
	    ASSERT(nnbr->directnumeles[m+1] == nnnbr->directnumeles[0]);
	    nnnbr_pc = nnnbr->chgs; /* panels in nnnbr */
	    nnnbr_dummy = nnnbr->nbr_is_dummy[0];
#if CHKDUM == ON
	    chkDummyList(nnnbr_pc, nnnbr_dummy, nnnsize);
#endif
	    for(i = nnsize - 1; i >= 0; i--) { /* loop on panels in nnbr */
	      if(nnbr_dummy[i]) continue;
	      if(nnbr_pc[i]->surf->type != DIELEC) {
		for(j = nnnsize - 1; j >= 0; j--) {
		  mat[offset + i][noffset+j] = nmat[i][j];
		}
	      }
	      else {
#if DPCOMP == ON
		fprintf(stdout, "Source mat, nnnbr to nnbr\n");
		dumpMat(nmat, nnsize, nnnsize);
#endif
		find_flux_density_row(mat, nmat, i, nnnsize, nnsize, offset,
				      noffset, nnbr_pc, nnnbr_pc, nnbr_dummy,
				      nnnbr_dummy);
	      }
	    }
	    noffset += nnnsize;
	  }
	}
	offset += nnsize;
      }
    }

    /* set up the local is_dummy vector for the rows/cols of mat */
    /* THIS COULD BE AVOIDED BY USING CUBE is_dummy's INSIDE invert() */
    if(big_mat_size < offset) {	/* allocate only if larger array needed */
      CALLOC(is_dummy, offset, int, ON, AMSC);
    }
    /* dump sections of the dummy vector in order cubes appear in nbr lst */
    /* (use fragment of Jacob's loop above) */
    nnnsize = noffset = nc->directnumeles[0];
    nc_dummy = nc->nbr_is_dummy[0];
    for(i = nnnsize - 1; i >= 0; i--) {
      is_dummy[i] = nc_dummy[i];
    }
    for(l = 0; l < nc->numnbrs; l++) {
      nnnbr = nc->nbrs[l];
      if(NEAR(nnnbr, nj, nk, nl)) {
	nnnsize = nnnbr->directnumeles[0];
	nc_dummy = nnnbr->nbr_is_dummy[0];
	for(i = nnnsize - 1; i >= 0; i--) {
	  is_dummy[i + noffset] = nc_dummy[i];
	}
	noffset += nnnsize;
      }
    }

    /* The big Matrix is filled in, invert it and get the preconditioner. */
#if DPCOMP == ON
    fprintf(stdout, "Before compression\n");
    dumpMat(mat, offset, offset);
#endif
    nnnsize = compressMat(mat, offset, is_dummy, BOTH);
#if DPCOMP == ON
    fprintf(stdout, "After compression\n");
    dumpMat(mat, nnnsize, nnnsize);
#endif
    invert(mat, nnnsize, NULL);
    expandMat(mat, offset, nnnsize, is_dummy, BOTH);
#if DPCOMP == ON
    fprintf(stdout, "After expansion\n");
    dumpMat(mat, offset, offset);
#endif

    /* Copy out the preconditioner to the saved matrices. */
    for(i = nsize - 1; i >= 0; i--) {
      for(j = nsize - 1; j >= 0; j--) {
	nc->precondmats[0][i][j] = mat[i][j];
      }
    }
    offset = nsize;
    for(k=0; k < nc->numnbrs; k++) {
      nnbr = nc->nbrs[k];
      if(NEAR(nnbr, nj, nk, nl)) {
	nnsize = nc->directnumeles[k+1];
	nmat = nc->precondmats[k+1];
	for(i = nsize - 1; i >= 0; i--) {
	  for(j = nnsize - 1; j >= 0; j--) {
	    nmat[i][j] = mat[i][offset + j];
	  }
	}
	offset += nnsize;
      }
      else nc->precondmats[k+1] = NULL;
    }
  }
}

/*
  finds a row of flux density coeffs from three potential coeff rows
  - to_mat[eval_row][] is the destination row; from_mat[eval_row][]
    initially contains the potential coefficients for evals at the 
    center of eval_panels[eval_row] (unless NUMDPT == 2, is garbage then)
  - the eval panels are scaned until eval_panels[eval_row]'s
    dummies are found and the corresponding two rows are identified
  - the divided differences built with entries in the same columns in
    these three rows replace the to_mat[eval_row][] entries:
    to_mat[eval_row][j] = a1*from_mat[eval_row][j]
              + a2*from_mat[pos_dum_row][j] + a3*from_mat[neg_dum_row][j]
  - if a dummy panel is not found in the panel list, its row is generated
    using explicit calcp() calls (shouldn't happen much)
  - global flags used here
    NUMDPT = number of divided diff points, 2 or 3
    SKIPDQ = ON=>don't do cancellation-prone add-subtract of identical
      influence of DIELEC/BOTH panels' charges on dummy panel pot. evals
*/
find_flux_density_row(to_mat, from_mat, eval_row, n_chg, n_eval, row_offset,
		      col_offset, eval_panels, chg_panels, eval_is_dummy, 
		      chg_is_dummy)
double **to_mat, **from_mat;
int eval_row, n_chg, n_eval, row_offset, col_offset;
charge **eval_panels, **chg_panels;
int *eval_is_dummy, *chg_is_dummy;
{
  int dindex, j;
  double factor, calcp();
  charge *dp;
  surface *surf = eval_panels[eval_row]->surf;

  /* do divided difference w/ three rows to get dielectric row */
#if NUMDPT == 3
  /* - dielectric panel row first */
  factor = -(surf->outer_perm + surf->inner_perm)/
      (eval_panels[eval_row]->pos_dummy->area);
#if DPDDIF == ON
  fprintf(stdout, "Center row, factor = %g\n", factor);
#endif
  for(j = n_chg - 1; j >= 0; j--) { /* loop on columns */
    if(!chg_is_dummy[j]) 
	to_mat[row_offset + eval_row][col_offset + j] 
	    = from_mat[eval_row][j]*factor;
#if DPDDIF == ON
    fprintf(stdout, " %.16e", from_mat[eval_row][j]);
#endif				/* #if DPDDIF == ON */
  }
#endif				/* #if NUMDPT == 3 */
  /* - do positive dummy row */
  /*   first find the dummy row */
  dindex = -1;
  dp = eval_panels[eval_row]->pos_dummy; /* get dummy panel from eval panel */
  for(j = n_eval - 1; j >= 0; j--) {
    if(!eval_is_dummy[j]) continue;
    if(dp == eval_panels[j]) {
      dindex = j;
      break;
    }
  }
  if(dindex != -1) { /* dummy row found */
#if NUMDPT == 3
    factor = surf->outer_perm/eval_panels[dindex]->area;
#else
    /* this is the only factor required for two dummy rows in two point case */
    factor = (surf->inner_perm - surf->outer_perm)
	/(eval_panels[eval_row]->neg_dummy->area 
	  + eval_panels[eval_row]->pos_dummy->area);
#endif
#if DPDDIF == ON
    fprintf(stdout, "\nPos dummy row, factor = %g\n", factor);
#endif
    for(j = n_chg - 1; j >= 0; j--) {
#if SKIPDQ == ON
      if(chg_panels[j]->index == eval_panels[eval_row]->index) {
	to_mat[row_offset + eval_row][col_offset + j] = 0.0;
	continue;
      }
#endif
      if(!chg_is_dummy[j])
#if NUMDPT == 3
	  to_mat[row_offset + eval_row][col_offset + j] 
	      += from_mat[dindex][j]*factor;
#else				/* make sure to overwrite possible garbage */
	  to_mat[row_offset + eval_row][col_offset + j] 
	      = -from_mat[dindex][j]*factor;
#endif
#if DPDDIF == ON
      fprintf(stdout, " %.16e (%d)", from_mat[dindex][j],chg_panels[j]->index);
#endif
    }
  }
  else {		/* dummy row out of cube => build it w/calcp */
#if NUMDPT == 3
    factor = surf->outer_perm/dp->area;
#else
    /* this is the only factor required for two dummy rows in two point case */
    factor = (surf->inner_perm - surf->outer_perm)
	/(eval_panels[eval_row]->neg_dummy->area 
	  + eval_panels[eval_row]->pos_dummy->area);
#endif
#if DPDDIF == ON
    fprintf(stdout, "\nPos dummy calcp row, factor = %g\n", factor);
#else
    fprintf(stderr, "\nolmulMatPrecond: building pos. dummy row\n");
#endif
    for(j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
      if(chg_panels[j]->index == eval_panels[eval_row]->index) {
	to_mat[row_offset + eval_row][col_offset + j] = 0.0;
	continue;
      }
#endif
      if(!chg_is_dummy[j]) {
#if NUMDPT == 3
	to_mat[row_offset + eval_row][col_offset + j] 
	    += calcp(chg_panels[j], dp, NULL)*factor;
#else
	to_mat[row_offset + eval_row][col_offset + j] 
	    = -calcp(chg_panels[j], dp, NULL)*factor;
#endif
#if DPDDIF == ON
	fprintf(stdout, " %.16e (%d)",
		calcp(chg_panels[j], dp, NULL),
		chg_panels[j]->index);
      }
      else {
	fprintf(stdout, " dummy");
#endif
      }
    }
  }
  /* - do negative dummy row */
  /*   first find the dummy row */
  dindex = -1;
  dp = eval_panels[eval_row]->neg_dummy; /* get dummy panel from eval panel */
  for(j = n_eval - 1; j >= 0; j--) {
    if(!eval_is_dummy[j]) continue;
    if(dp == eval_panels[j]) {
      dindex = j;
      break;
    }
  }
  if(dindex != -1) { /* dummy row found */
#if NUMDPT == 3
    factor = surf->inner_perm/eval_panels[dindex]->area;
#endif
#if DPDDIF == ON
    fprintf(stdout, "\nNeg dummy row, factor = %g\n", factor);
#endif
    for(j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
      if(chg_panels[j]->index == eval_panels[eval_row]->index) continue;
#endif
      if(!chg_is_dummy[j])
	  to_mat[row_offset + eval_row][col_offset + j] 
	      += from_mat[dindex][j]*factor;
#if DPDDIF == ON
      fprintf(stdout, " %.16e (%d)", from_mat[dindex][j],chg_panels[j]->index);
#endif
    }
  }
  else {		/* dummy row out of cube => build it w/calcp */
    factor = surf->inner_perm/dp->area;
#if DPDDIF == ON
    fprintf(stdout, "\nNeg dummy calcp row, factor = %g\n", factor);
#else
    fprintf(stderr, "olmulMatPrecond: building neg. dummy row\n");
#endif
    for(j = n_chg - 1; j >= 0; j--) {
#if SKIPQD == ON
      if(chg_panels[j]->index == eval_panels[eval_row]->index) continue;
#endif
      if(!chg_is_dummy[j]) {
	to_mat[row_offset + eval_row][col_offset + j] 
	    += calcp(chg_panels[j], dp, NULL)*factor;
#if DPDDIF == ON
	fprintf(stdout, " %.16e (%d)",
		calcp(chg_panels[j], dp, NULL),
		chg_panels[j]->index);
      }
      else {
	fprintf(stdout, " dummy");
#endif
      }
    }
  }
#if NUMDPT == 2
  /* - do row entry due to panel contribution 
     - entry only necessary if eval panel is in chg panel list */
  /*   search for the eval panel in the charge panel list */
  dp = NULL;
  for(j = n_chg - 1; j >= 0; j--) {
    if(!chg_is_dummy[j]) {
      if(eval_panels[eval_row] == chg_panels[j]) {
	dp = eval_panels[eval_row];
	break;
      }
    }
  }
  /*   set entry if eval panel found in chg panel list 
       - this is an overwrite; contributions of other rows should cancel */
  if(dp != NULL) {
    to_mat[row_offset + eval_row][col_offset + j]
	= -(2*M_PI*(surf->inner_perm + surf->outer_perm)
	    /eval_panels[eval_row]->area);
  }
#endif

#if DPDDIF == ON
  fprintf(stdout, "\nDivided difference row (%d)\n", 
	  eval_panels[eval_row]->index);
  for(j = n_chg - 1; j >= 0; j--) {
    fprintf(stdout, " %.16e (%d)", 
	    to_mat[row_offset + eval_row][col_offset + j], 
	    chg_panels[j]->index);
  }
  fprintf(stdout, "\n\n");
#endif
}


/* 
MulMatUp computes the multipole to multipole or charge to
multipole matrices that map to a parent's multipole coeffs from its
children's multipoles or charges. Note that only one set of
multipole to multipole matrices is computed per level by exploiting the
uniform break-up of three-space (ie many shifts have similar geometries).  
*/
mulMatUp(sys) 
ssystem *sys; 
{
cube *nextc, *kid;
int i, j, numterms, depth, order = sys->order;
double **multimats[8];

#if OFF == ON			/* OFF == OFF produces the M2M error!! */
    for(i=0; i < 8; i++) multimats[i] = NULL;
#endif

  numterms = multerms(order);

  if(sys->depth < 2) {
    /* fprintf(stdout, "\nWarning: no multipole acceleration\n");*/
    return;	/* return if upward pass not possible */
  }

/* Handle the lowest level cubes first (set up Q2M's). */
  for(nextc=sys->multilist[sys->depth]; nextc != NULL; nextc = nextc->mnext) {
    nextc->multisize = numterms;
    CALLOC(nextc->multi, numterms, double, ON, AMSC);
    CALLOC(nextc->upmats, 1, double**, ON, AMSC);
    nextc->upmats[0] = mulQ2Multi(nextc->chgs, nextc->nbr_is_dummy[0],
				  nextc->upnumeles[0],
				  nextc->x, nextc->y, nextc->z, order);

#if DISSYN == ON
    multicnt[nextc->level]++;
#endif

#if DMTCNT == ON
    Q2Mcnt[nextc->level][nextc->level]++;
#endif

  }

  if(sys->locallist[sys->depth] == NULL
     && sys->multilist[sys->depth] == NULL) {
    fprintf(stdout, "No expansions at level %d (lowest)\n", sys->depth);
    /*if(sys->depth < 3) 
	fprintf(stdout, " (Warning: no multipole acceleration)\n");*/
  }
  else if(sys->locallist[sys->depth] == NULL) {
    fprintf(stdout, "No local expansions at level %d (lowest)\n", sys->depth);
  }
  else if(sys->multilist[sys->depth] == NULL) {
    fprintf(stdout, "No multipole expansions at level %d (lowest)\n", 
	    sys->depth); 
    /*if(sys->depth < 3) 
	fprintf(stdout, " (Warning: no multipole acceleration)\n");*/
  }

/* Allocate the vectors and matrices for the cubes. */
  /* no multipoles over root cube or its kids (would not be used if made) */
  for(depth = (sys->depth - 1); depth > 1; depth--) {
  /* set up M2M's and Q2M's to compute the multipoles needed for this level */
    if(sys->locallist[depth] == NULL && sys->multilist[depth] == NULL) {
      fprintf(stdout, "No expansions at level %d\n", depth);
      /*if(depth < 3) fprintf(stdout, " (Warning: no multipole acceleration)\n");
      else fprintf(stdout, "\n");*/
    }
    else if(sys->locallist[depth] == NULL) {
      fprintf(stdout, "No local expansions at level %d\n", depth);
    }
    else if(sys->multilist[depth] == NULL) {
      fprintf(stdout, "No multipole expansions at level %d\n", depth); 
      /*if(depth < 3) fprintf(stdout, " (Warning: no multipole acceleration)\n");
      else fprintf(stdout, "\n");*/
    }

#if ON == ON			/* ON == OFF produces the M2M error!! */
    /* NULL out pointers to same-geometry M2M mats for this level */
    for(i=0; i < 8; i++) multimats[i] = NULL;
#endif

  /* Hit nonempty cubes at this level assigning ptrs to precomputed   */
  /* M2M mats (for this lev), or if a kid is exact, computing Q2M matrices. */
    for(nextc=sys->multilist[depth]; nextc != NULL; nextc = nextc->mnext) {
      
#if DISSYN == ON
      multicnt[nextc->level]++;
#endif

    /* Save space for upvector sizes, upvect ptrs, and upmats. */
      nextc->multisize = numterms;
      if(numterms > 0) {
	CALLOC(nextc->multi, numterms, double, ON, AMSC);
      }
      if(nextc->upnumvects) {
	CALLOC(nextc->upnumeles, nextc->upnumvects, int, ON, AMSC);
	CALLOC(nextc->upvects, nextc->upnumvects, double*, ON, AMSC);
	CALLOC(nextc->upmats, nextc->upnumvects, double**, ON, AMSC);
      }
    /* Go through nonempty kids and fill in upvectors and upmats. */
      for(i=0, j=0; j < nextc->numkids; j++) {
	if((kid = nextc->kids[j]) != NULL) { /* NULL => empty kid cube */
	  if(kid->mul_exact == FALSE) { /* if kid has a multi */
	    nextc->upvects[i] = kid->multi;
	    nextc->upnumeles[i] = kid->multisize;
	    if(multimats[j] == NULL) { /* Build the needed matrix only once. */
	      multimats[j] = mulMulti2Multi(kid->x, kid->y, kid->z, nextc->x, 
					    nextc->y, nextc->z, order);
	    }
	    nextc->upmats[i] = multimats[j];

#if DMTCNT == ON
	    M2Mcnt[kid->level][nextc->level]++; /* cnts usage, ~computation */
#endif				/* only at most 8 mats really built/level */

	  }
	  else {		/* if kid is exact, has no multi */
	    nextc->upvects[i] = kid->upvects[0];
	    nextc->upnumeles[i] = kid->upnumeles[0];
	    nextc->upmats[i] = mulQ2Multi(kid->chgs, kid->nbr_is_dummy[0],
					  kid->upnumeles[0],
					  nextc->x, nextc->y, nextc->z, order);

#if DMTCNT == ON
	    Q2Mcnt[kid->level][nextc->level]++;
#endif

	  }
	  i++;			/* only increments if kid is not empty */
	}
      }
    }
  }
}

/*
  builds the transformation matrices for the final evaluation pass (M2P, L2P)
  for all the direct list (generated by linkcubes(), non-empty
  lowest level cubes) cubes:

  for each cube A in the direct list:
  1) if the cube is not exact (always the case if ADAPT = OFF)
     a) and if DNTYPE = GRENGD build an L2P matrix from A to A 
     b) and if DNTYPE = NOSHFT build an L2P matrix from each of A's ancestors
        with level > 1 (including A) to A
     c) and if DNTYPE = NOLOCL build an M2P matrix from each of A's fake 
        ilist entries to A (same action as 2b)
  2) if the cube is exact, find the 1st ancestor of A, cube B, 
     which either is not exact and is at level 2,3,4... or is at level 1
     a) if B is at level 2,3,4... 
        i) if DNTYPE = GRENGD, construct an L2P from B to A and M2P's
           from the cubes in the true interaction lists of A and all its
	   ancestors up to and including B (a partial fake interaction list)
	j) if DNTYPE = NOSHFT, find cube C, the ancestor of B at level 1;
	   construct L2P's from the ancestors of B (including B but not C)
	   to A and Q- or M2P's from the cubes in the true interaction lists 
	   of A and all its ancestors up to and including B (a partial fake 
	   interaction list)
	k) if DNTYPE = NOLOCL, do 2b
     b) if B is at level 1 construct M2P's for all the cubes in A's
        fake interaction list

  true interaction list - RADINTER = OFF, those sibling
  (same level) cubes of a given cube who are children of the neighbors
  of the given cube's parent and are not neighbors of the given cube 
  - ie those cubes required to cover charges well separated from the given
  cube but not accounted for in the parent's local expansion 
  - the flag NNBRS is the number of sibling cube "shells" taken as neighbors 
  
  fake interaction list - RADINTER = OFF, the combined true interaction lists
  of a given cube and all its ancestors at levels 2,3,4...

  if RADINTER = ON, any 8 siblings of the given cube which form a well 
  separated cube one level up are included in the lists as a single higher
  level cube
  
  if ADAPT = OFF, no cube is exact so step 2 is never done

  this routine is used alone if compiled with DNTYPE = NOLOCL or after
  mulMatDown, which produces M2L and L2L matrices (DNTYPE = GRENGD) or
  just M2L matrices (DNTYPE = NOSHFT) --  DNTYPE = GRENGD does the full
  Greengard hiearchical downward pass

*/
void mulMatEval(sys)
ssystem *sys;
{
  int i, j, k, ttlvects, vects;
  cube *na, *nc, *nexti;

  if(sys->depth < 2) return;	/* ret if upward pass not possible/worth it */

  for(nc = sys->directlist; nc != NULL; nc = nc->dnext) {

    ASSERT(nc->level == sys->depth);
    ASSERT(nc->upnumvects > 0);

    /* allocate space for evaluation pass vectors; check nc's ancestors */
    /* First count the number of transformations to do. */
    for(na = nc, ttlvects = 0; na->level > 1; na = na->parent) { 
      if(na->loc_exact == FALSE && DNTYPE != NOLOCL) {
	ttlvects++;  /* allow for na to na local expansion (L2P) */
	if(DNTYPE == GRENGD) break; /* Only one local expansion if shifting. */
      }
      else {
	ttlvects += na->interSize; /* room for Q2P and M2P xformations */
      }
    }
    nc->evalnumvects = ttlvects; /* save ttl # of transformations to do */
    if(ttlvects > 0) {
      CALLOC(nc->evalvects, ttlvects, double*, ON, AMSC);
      CALLOC(nc->evalnumeles, ttlvects, int, ON, AMSC);
      CALLOC(nc->evalmats, ttlvects, double**, ON, AMSC);
      CALLOC(nc->eval_isQ2P, ttlvects, int, ON, AMSC);
    }
    
#if DILIST == ON
    fprintf(stdout, "\nInteraction list (%d entries) for ", ttlvects);
    disExParsimpcube(nc);
#endif
    
    /* set up exp/charge vectors and L2P, Q2P and/or M2P matrices as req'd */
    for(j=0, na = nc, ttlvects = 0; na->level > 1; na = na->parent) { 
      if(na->loc_exact == FALSE && DNTYPE != NOLOCL) {  
	/* build matrices for local expansion evaluation */
	nc->evalmats[j] = mulLocal2P(na->x, na->y, na->z, nc->chgs,
				     nc->upnumeles[0], sys->order);
	nc->evalnumeles[j] = na->localsize;
	nc->evalvects[j] = na->local;
/*	add_to_counts(nc, na->localsize, sys->evalL2Ps, sys->cntL2Ps); */
	nc->eval_isQ2P[j] = FALSE;
	j++; 
	
#if DMTCNT == ON
	L2Pcnt[na->level][nc->level]++;
#endif
	
#if DILIST == ON
	fprintf(stdout, "L2P: ");
	disExtrasimpcube(na);
#endif
	if(DNTYPE == GRENGD) break; /* Only one local expansion if shifting. */
      }
      else { /* build matrices for ancestor's (or cube's if 1st time) ilist */
	for(i=0; i < na->interSize; i++) {
	  nexti = na->interList[i];
	  if(nexti->mul_exact == TRUE) {
	    nc->evalvects[j] = nexti->upvects[0];
	    nc->evalmats[j] = Q2P(nexti->chgs, nexti->upnumeles[0], 
				  nexti->nbr_is_dummy[0], nc->chgs, 
				  nc->upnumeles[0], TRUE);
	    /* this is a hack to fix the fact that direct stuff is */
	    /* don't done componentwise */

	    /* obsolete as of 12/92  due to eval_isQ2P stuff
	    fixEvalDirect(nexti->chgs, nexti->upnumeles[0], 
				  nexti->nbr_is_dummy[0], nc->chgs, 
				  nc->upnumeles[0], nc->evalmats[j]);
	    */

	    nc->evalnumeles[j] = nexti->upnumeles[0];
	    nc->eval_isQ2P[j] = TRUE;
/*	    add_to_counts(nc, nexti->upnumeles[0], sys->evalQ2Ps, sys->cntQ2Ps); */
	    j++;

#if DMTCNT == ON
	    Q2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
	    fprintf(stdout, "Q2P: ");
	    disExtrasimpcube(nexti);
#endif
	  }
	  else {
	    nc->evalvects[j] = nexti->multi;
	    nc->evalmats[j] = mulMulti2P(nexti->x, nexti->y, nexti->z, 
					 nc->chgs, nc->upnumeles[0], 
					 sys->order);
	    nc->evalnumeles[j] = nexti->multisize;
	    nc->eval_isQ2P[j] = FALSE;
/*	    add_to_counts(nc, nexti->multisize, sys->evalM2Ps, sys->cntM2Ps);*/
	    
	    j++;
	    
#if DMTCNT == ON
	    M2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
	    fprintf(stdout, "M2P: ");
	    disExtrasimpcube(nexti);
#endif
	  }
	}
      }
    }
  }
}


/* 
  sets up matrices for the downward pass
  For each cube in local list (parents always in list before kids):
  1) parent's local to child's local unless DNTYPE = NOSHFT or no parent local
  2) multipoles for (Parent+parent's nbrs - child nbrs) to child's local
  -eval is sum of ancestral local evals for each lowest lev cube if NOSHFT
    otherwise only lowest level local is evaluated (see mulMatEval)
  -with ADAPT = OFF no cube is exact so local list is all non-empty cube lev>1
  -mats that give potentials (M2P, L2P, Q2P) are calculated in mulMatEval()
  -this routine makes only L2L, M2L and Q2L matrices
*/
mulMatDown(sys)
ssystem *sys;
{
  int i, j, vects;
  cube *nc, *parent, *ni;
  int depth;

  ASSERT(DNTYPE != NOLOCL);	/* use mulMatEval() alone if NOLOCL */

  for(depth = 2; depth <= sys->depth; depth++) { /* no locals before level 2 */
    for(nc=sys->locallist[depth]; nc != NULL; nc = nc->lnext) {

      /* Allocate for interaction list, include one for parent if needed */
      if((depth <= 2) || (DNTYPE == NOSHFT)) vects = nc->interSize;
      else vects = nc->interSize + 1;
      nc->downnumvects = vects;
      if(vects > 0) {
	CALLOC(nc->downvects, vects, double*, ON, AMSC);
	CALLOC(nc->downnumeles, vects, int, ON, AMSC);
	CALLOC(nc->downmats, vects, double**, ON, AMSC);
      }

      parent = nc->parent;
      ASSERT(parent->loc_exact == FALSE); /* has >= #evals of any of its kids*/

#if DISSYN == ON
      localcnt[nc->level]++;
#endif

      if((depth <= 2) || (DNTYPE == NOSHFT)) i = 0; /* No parent local. */
      else { /* Create the mapping matrix for the parent to kid. */
	i = 1;

	nc->downmats[0] = mulLocal2Local(parent->x, parent->y, parent->z,
					 nc->x, nc->y, nc->z, sys->order);
	nc->downnumeles[0] = parent->localsize;
	nc->downvects[0] = parent->local;

#if DMTCNT == ON
	L2Lcnt[parent->level][nc->level]++;
#endif
      }

      /* Go through the interaction list and create mapping matrices. */
      for(j = 0; j < nc->interSize; j++, i++) {
	ni = nc->interList[j];
	if(ni->mul_exact == TRUE) {	/* ex->ex (Q2P) xforms in mulMatEval */
	  nc->downvects[i] = ni->upvects[0];
	  nc->downmats[i] = mulQ2Local(ni->chgs, ni->upnumeles[0],
				       ni->nbr_is_dummy[0],
				       nc->x, nc->y, nc->z, sys->order);
	  nc->downnumeles[i] = ni->upnumeles[0];
#if DMTCNT == ON
	  Q2Lcnt[ni->level][nc->level]++;
#endif
	}
	else {
	  nc->downvects[i] = ni->multi;
	  nc->downmats[i] = mulMulti2Local(ni->x, ni->y, ni->z, nc->x, 
					   nc->y, nc->z, sys->order);
	  nc->downnumeles[i] = ni->multisize;
#if DMTCNT == ON
	  M2Lcnt[ni->level][nc->level]++;
#endif
	}
      }
    }
  }
}




