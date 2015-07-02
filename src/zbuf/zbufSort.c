/*
Copyright (c) 1990 Massachusetts Institute of Technology, Cambridge, MA.
All rights reserved.

This Agreement gives you, the LICENSEE, certain rights and obligations.
By using the software, you indicate that you have read, understood, and
will comply with the terms.

Permission to use, copy and modify for internal, noncommercial purposes
is hereby granted.  Any distribution of this program or any part thereof
is strictly prohibited without prior written consent of M.I.T.

Title to copyright to this software and to any associated documentation
shall at all times remain with M.I.T. and LICENSEE agrees to preserve
same.  LICENSEE agrees not to make any copies except for LICENSEE'S
internal noncommercial use, or to use separately any portion of this
software without prior written consent of M.I.T.  LICENSEE agrees to
place the appropriate copyright notice on any such copies.

Nothing in this Agreement shall be construed as conferring rights to use
in advertising, publicity or otherwise any trademark or the name of
"Massachusetts Institute of Technology" or "M.I.T."

M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By
way of example, but not limitation, M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS OR DOCUMENTATION WILL
NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
M.I.T. shall not be held liable for any liability nor for any direct,
indirect or consequential damages with respect to any claim by LICENSEE
or any third party on account of or arising from this Agreement or use
of this software.
*/

#include "mulGlobal.h"
#include "zbufGlobal.h"

int cnt;			/* used in setting up the depth graph */

/*
  returns TRUE if difference of two doubles is within machine precision
*/
int diff_is_zero(num1, num2, bias)
double num1, num2, bias;
{
  double margin;

  margin = MAX(fabs(num1), fabs(num2));
  
  if(fabs(num1-num2) <= margin*MARGIN*bias + MARGIN*bias) return(TRUE);
  return(FALSE);
}
    

/*
  returns TRUE if diff. of two doubles is negative within machine precision
  - note that argument order is relevant
*/
int diff_is_negative(num1, num2, bias)
double num1, num2, bias;
{
  double margin;

  margin = MAX(fabs(num1), fabs(num2));
  
  if(num1-num2 < -margin*MARGIN*bias - MARGIN*bias) return(TRUE);
  return(FALSE);
}
    
/*
  returns the dot product
*/
double dot(vec1, vec2)
double *vec1, *vec2;
{
  int i;
  double ret = 0.0;
  
  for(i = 0; i < 3; i++) ret += vec1[i]*vec2[i];
  return(ret);
}

/*
  calculate the cross product out = in1Xin2
*/
void crossProd(out, in1, in2)
double *out, *in1, *in2;
{
  out[0] = in1[1]*in2[2] - in1[2]*in2[1];
  out[1] = in1[2]*in2[0] - in1[0]*in2[2];
  out[2] = in1[0]*in2[1] - in1[1]*in2[0];
}

/*
  returns a normal and rhs for a plane given by three points
*/
double getPlane(normal, p1, p2, p3)
double *normal, *p1, *p2, *p3;
{
  int i;
  double dot(), v12[3], v13[3], norm;

  /* set up side vectors */
  for(i = 0; i < 3; i++) {
    v12[i] = p2[i]-p1[i];
    v13[i] = p3[i]-p1[i];
  }

  /* get the normal using a cross product */
  crossProd(normal, v12, v13);
  norm = sqrt(dot(normal, normal));
  for(i = 0; i < 3; i++) normal[i] = normal[i]/norm;

  /* set up the rhs using p1 as a point on the plane */
  return(dot(normal, p1));
}

/*
  returns POS, NEG or SPLIT for which side face corners are rel to plane
*/
int whichSide(fac, facplane)
face *fac, *facplane;
{
  int i, neg, pos, zero;
  double value[MAXSIDES];      	/* holds values when subbed in plane equ */
  double dot(), margin[MAXSIDES], temp[MAXSIDES];

  /* if planes are identical, return SAME right away */
  if(fac->rhs < facplane->rhs+MARGIN && fac->rhs > facplane->rhs-MARGIN && 
     fac->normal[0] < facplane->normal[0]+MARGIN &&
     fac->normal[0] > facplane->normal[0]-MARGIN &&
     fac->normal[1] < facplane->normal[1]+MARGIN && 
     fac->normal[1] > facplane->normal[1]-MARGIN && 
     fac->normal[2] < facplane->normal[2]+MARGIN &&
     fac->normal[2] > facplane->normal[2]-MARGIN) return(SAME);

  /* sub all the points of fac into the equation for the plane of other
     - POS if result of all subs is positive or zero (normals pnt to view)
     - NEG if result of all subs is negative or zero
     - SPLIT whenever at least one positive and at least one negative
     - SAME if faces are in the same plane */
  for(i = 0; i < fac->numsides; i++) {
    temp[i] = dot(facplane->normal, fac->c[i]);
    value[i] = temp[i] - facplane->rhs;
    margin[i] = MARGIN*MAX(fabs(facplane->rhs), fabs(temp[i]));
  }
  /* count the number of each type of corner */
  for(i = neg = pos = zero = 0; i < fac->numsides; i++) {
    if(value[i] >= -margin[i] && value[i] <= margin[i]) zero++;
    else if(value[i] > 0.0) pos++;
    else if(value[i] < 0.0) neg++;
    else {
      fprintf(stderr, 
	      "\nwhichSide: strange corner, value = %g, rhs = %g, dot = %g\n",
	      value[i], facplane->rhs, temp[i]);
      exit(0);
    }
  }
  if(neg > 0 && pos == 0) return(NEG);
  else if(pos > 0 && neg == 0) return(POS);
  else if(pos > 0 && neg > 0) return(SPLIT);
  else if(zero == fac->numsides) return(SAME);
  else {
    fprintf(stderr,"\nwhichSide: face has %d corners, %d pos %d neg %d zero\n",
	    fac->numsides, pos, neg, zero);
    exit(0);
  }
}

/*
  returns TRUE if lines from1-to1 and from2-to2 intersect
  - if lines are parallel up to MARGIN, FALSE is returned
  - lines than hit at segment endpoints up to MARGIN return FALSE
*/
int doLinesIntersect(isect, from1, to1, from2, to2)
double *from1, *to1, *from2, *to2, *isect;
{
  double A[2][2], b[2];
  double det, margin, temp1, temp2, margin1, margin2;

#if DEBUGX == ON
  fprintf(stdout, "seg1 = (%g %g) (%g %g) seg2 = (%g %g) (%g %g)\n",
	  from1[0], from1[1], to1[0], to1[1], from2[0], from2[1], 
	  to2[0], to2[1]);
#endif

  /* build the 2x2 system to solve for alpha1, alpha2 in
     from1 + alpha1*(to1 - from1) = from2 + alpha2*(to2 - from2) */
  A[0][0] = to1[0]-from1[0];
  A[1][0] = to1[1]-from1[1];
  A[0][1] = from2[0]-to2[0];
  A[1][1] = from2[1]-to2[1];
  b[0] = from2[0]-from1[0];
  b[1] = from2[1]-from1[1];

  /* attempt to solve, checking carefully for cancellation in the det */
  temp1 = A[0][0]*A[1][1];
  temp2 = A[1][0]*A[0][1];
  /* margin = MARGIN*MAX(temp1, temp2); */
  det = temp1 - temp2;
  if(diff_is_zero(temp1, temp2, 1.0)) return(FALSE);	/* check for || */

  /* solve */
  /* margin1 = MARGIN*MAX(A[1][1]*b[0], A[0][1]*b[1]);
  margin2 = MARGIN*MAX(A[0][0]*b[1], A[1][0]*b[0]); */
  temp1 = A[1][1]*b[0]-A[0][1]*b[1]; /* temp = alpha*det */
  temp2 = A[0][0]*b[1]-A[1][0]*b[0];

  /* check for end hits => no crossing (just barely touching) */
  if(diff_is_zero(A[1][1]*b[0], A[0][1]*b[1], 1.0) ||
     diff_is_zero(A[1][1]*b[0], A[0][1]*b[1]+det, 1.0)) return(FALSE);
  if(diff_is_zero(A[0][0]*b[1], A[1][0]*b[0], 1.0) ||
     diff_is_zero(A[0][0]*b[1], A[1][0]*b[0]+det, 1.0)) return(FALSE);

  /* set up intersection point using isect = from1 + alpha1*(to1 - from1) */
  isect[0] = from1[0] + temp1*(to1[0] - from1[0])/det;
  isect[1] = from1[1] + temp1*(to1[1] - from1[1])/det;
  isect[2] = 0.0;		/* not really necessary */

  /* check for alphas between 0 and 1 (temps btwn 0 and det) */
  if(det > 0.0) {
    if(0.0 < temp1 && temp1 < det && 0.0 < temp2 && temp2 < det) return(TRUE);
    else return(FALSE);
  }
  else {
    if(det < temp1 && temp1 < 0.0 && det < temp2 && temp2 < 0.0) return(TRUE);
    else return(FALSE);
  }
}
  
/*
  returns TRUE if one panel is completely inside the other
  - does careful checks when vertices on one face are on verts, sides of other
*/
int face_is_inside(corners1, ccnt1, corners2, ccnt2, com_pnt)
double **corners1, **corners2, *com_pnt;
int ccnt1, ccnt2;
{
  int i, j, k, n, ccnt, zeros, pos, neg, ncnt;
  double *sfrom, *sto, margin1, margin2, **refcor, **curcor;
  double innerpnt[2], side[2], pnt[2];
  double temp;

  /* do cross products between sides of one face and vectors to first point
     of the other---must get same sign or zero for all if one inside other 
     - assumes faces have no sides that intersect 
     - do assuming both sides are outside to be sure to catch all cases */
  for(refcor = corners1, curcor = corners2, ccnt = ccnt1, ncnt = ccnt2, n = 0;
      n < 2; 
      refcor = corners2, curcor = corners1, ccnt = ccnt2, ncnt = ccnt1, n++) {
    zeros = pos = neg = 0;	/* cross-product value counts */
    for(i = 0; i < ccnt; i++) {	/* loop on "outer" face sides */
      if(i == ccnt-1) j = 0;
      else j = i+1;
      /* figure the center point of the non-traversed panel */
      innerpnt[0] = innerpnt[1] = 0.0;
      for(k = 0; k < ncnt; k++) {
	innerpnt[0] += curcor[k][0];
	innerpnt[1] += curcor[k][1];
      }
      innerpnt[0] /= (double)ncnt;
      innerpnt[1] /= (double)ncnt;
      /* figure side and point vector */
      side[0] = refcor[j][0] - refcor[i][0];
      side[1] = refcor[j][1] - refcor[i][1];
      pnt[0] = innerpnt[0] - refcor[i][0];
      pnt[1] = innerpnt[1] - refcor[i][1];
      margin1 = MARGIN*MAX(fabs(innerpnt[0]), fabs(refcor[i][0]));
      margin2 = MARGIN*MAX(fabs(innerpnt[1]), fabs(refcor[i][1]));
      /* figure cross-product---if a cross-product would be zero, skip it */
      if(diff_is_zero(innerpnt[0], refcor[i][0], 1.0)
	 && diff_is_zero(innerpnt[1], refcor[i][1], 1.0)) {
	zeros++;
	continue;
      }
      /* now check if vertex of "inside" face is on a side, not on corner */
      if(diff_is_zero(side[0]*pnt[1], pnt[0]*side[1], 1.0)) { /* is X-prod 0 */
	zeros++;		/*       due to vertex on side (not corner) */
      }
      else if(diff_is_negative(side[0]*pnt[1], pnt[0]*side[1], 1.0)) {
	neg++;
      }
      else pos++;		/* solid positive cross-product */
    }
    /* see if negative and positive are not mixed together => full overlap */
    if(((neg == 0 && pos != 0) || (neg != 0 && pos == 0)) && zeros == 0) {
      for(k = 0; k < 3; k++) com_pnt[k] = innerpnt[k];
      return(TRUE);
    }
  }

  /* wasn't any full overlap in either case => no overlap */
  return(FALSE);
	
}
  
#if 1 == 0			/* old is1stFaceDeeper */
/*
  returns TRUE if fac is deeper than facref and facref overlaps fac
  returns FALSE if facref has no overlap with fac
  returns REVERSE if facref is deeper than fac and fac overlaps facref
  - checks for intersections between each line of fac and all sides of facref
  - also checks for complete overlap (one face inside the other)
*/
int is1stFaceDeeper(fac, facref, view)
face *fac, *facref;
double *view;
{
  int i, j, k, olap[2], is_overlap;
  static double **cproj = NULL;	/* corners of fac in facref's plane */
  static double **cref = NULL;	/* corners of facref in facref plane */
  double alpha[MAXSIDES];	/* 1 => view point 0 => corner */
  double minref[2], maxref[2];	/* bounding box coordinates */
  double minfac[2], maxfac[2];
  double x[3], y[3];		/* coordinates of x and y in facref plane */
  double dot(), temp, tvec[3], tvec1[3], margin, ovrlapmgn = 0.0, temp1;
  double margin1;
  double *cfr, *ctr, *cff, *ctf, avg[3];
  int all_pos, all_neg, intersect;

  if(fac->index == 5 && facref->index == 14
     || fac->index == 14 && facref->index == 5) {
    i = i + 0;
  }

  /* allocate for local arrays on first call */
  if(cproj == NULL) {
    CALLOC(cproj, MAXSIDES, double *, ON, AMSC);
    CALLOC(cref, MAXSIDES, double *, ON, AMSC);
    for(i = 0; i < MAXSIDES; i++) {
      CALLOC(cproj[i], 3, double, ON, AMSC);
      CALLOC(cref[i], 3, double, ON, AMSC);
    }
  }

  /* figure if panels are in the same plane
     - if they are, they cant overlap if this is a legal discretization
       so return false */
  temp = temp1 = margin = 0.0;
  for(i = 0; i < 3; i++) {
    temp1 = fac->normal[i]-facref->normal[i];
    margin = MAX(margin, fabs(fac->normal[i]));
    margin = MAX(margin, fabs(facref->normal[i]));
    temp += (temp1*temp1);
  }
  temp = sqrt(temp);		/* get norm of difference of normals */
  margin *= MARGIN;		/* get norm of difference margin */
  /* check rhs and normal equivalence (panels in same plane) */
  if(diff_is_zero(fac->rhs, facref->rhs, 1.0) 
     && -margin <= temp && temp <= margin) {
    return(FALSE);
  }

  /* figure projections of fac's corners back to facref's plane rel to view */
  for(i = 0; i < fac->numsides; i++) {
    for(j = 0; j < 3; j++) tvec[j] = view[j] - fac->c[i][j]; /* get v-c */
    temp = dot(facref->normal, tvec); /* get n.(v-c) */
    margin = sqrt(dot(tvec, tvec))*MARGIN;
    /* test fails if v-c is perpendicular to n */
    if(temp > -margin && temp < margin) {
      fprintf(stderr, 
	      "is1stFaceDeeper: Warning, view-corner in view plane (?)\n");
      return(FALSE);
    }
    /* get alpha as in c + alpha*(v-c) = c' */
    alpha[i] = (facref->rhs - dot(facref->normal, fac->c[i]))/temp;
    for(j = 0; j < 3; j++)	/* get c' */
	cproj[i][j] = (1.0-alpha[i])*fac->c[i][j]+alpha[i]*view[j];
  }

  /* figure x and y coordinates in facref plane (normal is z coordinate) */
  /* x = c0-c1 always */
  for(j = 0; j < 3; j++) x[j] = facref->c[0][j] - facref->c[1][j];
  temp = sqrt(dot(x, x));
  for(j = 0; j < 3; j++) x[j] /= temp; /* normalize */
  /* y = zXx */
  crossProd(y, facref->normal, x);
  
  /* project all fac corner projections onto new x and y coordinates
     - facref->c[0] plays the role of origin */
  for(i = 0; i < fac->numsides; i++) {
    for(j = 0; j < 3; j++) tvec1[j] = cproj[i][j] - facref->c[0][j];
    tvec[0] = dot(x, tvec1);	/* get weight in x direction */
    tvec[1] = dot(y, tvec1);	/* get weight in y direction */
    tvec[2] = 0.0;		/* all z weights must = facref->rhs */
    for(j = 0; j < 3; j++) cproj[i][j] = tvec[j]; /* xfer */
  }
  for(j = 0; j < 3; j++) cref[0][j] = 0.0;
  for(i = 1; i < facref->numsides; i++) {
    for(j = 0; j < 3; j++) tvec1[j] = facref->c[i][j] - facref->c[0][j];
    cref[i][0] = dot(x, tvec1);	/* get weight in x direction */
    cref[i][1] = dot(y, tvec1);	/* get weight in y direction */
    cref[i][2] = 0.0;		/* all z weights must = facref->rhs */
  }

  /* up to here is identical to isThereBoxOverlap() */
  /* for each side of facref, see if there is an intersect. w/ all fac sides */
#if DEBUGX == ON
  fprintf(stdout, "Is face %d behind face %d?\n", fac->index, facref->index);
#endif
  is_overlap = FALSE;
  for(i = 0; i < facref->numsides && !is_overlap; i++) {
    cfr = cref[i];
    if(i == facref->numsides - 1) ctr = cref[0];
    else ctr = cref[i+1];
    for(j = 0; j < fac->numsides && !is_overlap; j++) {
      cff = cproj[j];
      if(j == fac->numsides - 1) ctf = cproj[0];
      else ctf = cproj[j+1];
      if((intersect = doLinesIntersect(cfr, ctr, cff, ctf)) == TRUE)
	  is_overlap = TRUE;
#if DEBUGX == ON
      fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
    }
  }
#if XOVTST == ON		/* do either this or face_is_inside() below */
  /* check for overlap with lines across face */
  if(facref->numsides == 4) {
    for(i = 0; i < 2; i++) {
      if(i == 0) {
	cfr = cref[0];
	ctr = cref[2];
      }
      else {
	cfr = cref[1];
	ctr = cref[3];
      }
      for(j = 0; j < fac->numsides; j++) {
	cff = cproj[j];
	if(j == fac->numsides - 1) ctf = cproj[0];
	else ctf = cproj[j+1];
	if((intersect = doLinesIntersect(cfr, ctr, cff, ctf)) == TRUE)
	    is_overlap = TRUE;
#if DEBUGX == ON
	fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
      }
    }
  }
  else if(facref->numsides == 3) {
    avg[0] = (cref[0][0] + cref[1][0] + cref[2][0])/3.0;
    avg[1] = (cref[0][1] + cref[1][1] + cref[2][1])/3.0;
    avg[2] = (cref[0][2] + cref[1][2] + cref[2][2])/3.0;
    cfr = avg;
    for(i = 0; i < 3; i++) {
      ctr = cref[i];
      for(j = 0; j < fac->numsides; j++) {
	cff = cproj[j];
	if(j == fac->numsides - 1) ctf = cproj[0];
	else ctf = cproj[j+1];
	if((intersect = doLinesIntersect(cfr, ctr, cff, ctf)) == TRUE)
	    is_overlap = TRUE;
#if DEBUGX == ON
	fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
      }
    }
  }
  else {
    fprintf(stderr, "isThereProjOverlap: can't handle %d side panel\n",
	    facref->numsides);
    exit(0);
  }
#else
  /* no sides overlap---check if one face completely obscures the other */
  if(!is_overlap) {
    is_overlap = face_is_inside(cref, facref->numsides, cproj, fac->numsides);
  }
#endif				/* XOVTST == ON */

  /* if no overlap, no edge in graph in any case */
  if(!is_overlap) return(FALSE);

  /* return TRUE only if fac is deeper than facref */
  /* return REVERSE only if facref is deeper than fac */
  /* if all alphas are positive or zero, fac is deeper (note use of margin) */
  /* if all alphas are negative or zero, facref is deeper */
  /* otherwise the test is inconclusive */
  all_pos = all_neg = TRUE;
  for(i = 0; i < fac->numsides; i++) {
    if(alpha[i] < -margin) all_pos = FALSE;
    else if(alpha[i] > margin) all_neg = FALSE;
  }

  if(all_pos && all_neg) {
    fprintf(stderr, 
	 "\nis1stFaceDeeper: Warning, all projection coefficients are zero\n");
    return(FALSE);
  }
  else if(all_pos) return(TRUE);
  else if(all_neg) return(REVERSE);
  else return(FALSE);
}

#else				/* new is1stFaceDeeper--proj to view plane */

/*
  returns TRUE if fac is deeper than facref and facref overlaps fac
  returns FALSE if facref has no overlap with fac
  returns REVERSE if facref is deeper than fac and fac overlaps facref
  - checks for intersections between each line of fac and all sides of facref
  - also checks for complete overlap (one face inside the other)
*/
int is1stFaceDeeper(fac, facref, view, rhs, normal)
face *fac, *facref;
double *view, rhs, *normal;
{
  int i, j, k, olap[2], is_overlap, isect_cnt;
  static double ***cproj = NULL;	/* corners of faces in view plane */
  double alpha[2][MAXSIDES];	/* 1 => view point 0 => corner */
  double minref[2], maxref[2];	/* bounding box coordinates */
  double minfac[2], maxfac[2];
  double x[3], y[3];		/* coordinates of x and y in facref plane */
  double dot(), temp, tvec[3], tvec1[3], margin, ovrlapmgn = 0.0, temp1;
  double margin1;
  double *cfr, *ctr, *cff, *ctf, avg[3], origin[3];
  double isect_avg[3], isect[3]; /* intersection points */
  double alpha_fac, alpha_facref, alpha_origin, alpha_x;
  int all_pos, all_neg, intersect, same_normal;
  face *curf;

  /* i don't know why this is here, but it does nothing, so i'll leave it */
  /*
  if(fac->index == 8 && facref->index == 20
     || fac->index == 20 && facref->index == 8) {
    i = i + 0;
  }
  */

  /* allocate for local arrays on first call */
  if(cproj == NULL) {
    CALLOC(cproj, 2, double **, ON, AMSC);
    for(k = 0; k < 2; k++) {
      CALLOC(cproj[k], MAXSIDES, double *, ON, AMSC);
      for(i = 0; i < MAXSIDES; i++) {
	CALLOC(cproj[k][i], 3, double, ON, AMSC);
      }
    }
  }

  /* figure if panels are in the same plane
     - if they are, they cant overlap if this is a legal discretization
       so return false */
#if 1 == 0
  temp = temp1 = margin = 0.0;
  for(i = 0; i < 3; i++) {
    temp1 = fac->normal[i]-facref->normal[i];
    temp1 = temp1*temp1;
    margin = MAX(margin, temp1);
    temp += temp1;
  }
  temp = sqrt(temp);		/* get norm of difference of normals */
  margin *= (MARGIN*PARMGN);	/* get norm of difference margin */
#endif
  same_normal = TRUE;
  for(i = 0; i < 3 && same_normal; i++) {
    if(!diff_is_zero(fac->normal[i], facref->normal[i], PARMGN)) 
	same_normal = FALSE;
  }
    
  /* check rhs and normal equivalence (panels in same plane) */
  if(diff_is_zero(fac->rhs, facref->rhs, PARMGN) && same_normal) return(FALSE);

  /* find projections of fac and facref corners onto view plane rel to view */
  for(curf = fac, k = 0; k < 2; k++, curf = facref) {
    for(i = 0; i < curf->numsides; i++) {
      for(j = 0; j < 3; j++) tvec[j] = view[j] - curf->c[i][j]; /* get v-c */
      temp = dot(normal, tvec); /* get n.(v-c) */
      margin = sqrt(dot(tvec, tvec))*MARGIN;
      ovrlapmgn = MAX(margin, ovrlapmgn);	/* used below */
      /* test fails if v-c is perpendicular to n */
      if(temp > -margin && temp < margin) return(FALSE);
      /* get alpha as in c + alpha*(v-c) = c' */
      alpha[k][i] = (rhs - dot(normal, curf->c[i]))/temp;
      for(j = 0; j < 3; j++)	/* get c' */
     /*cproj[k][i][j] = (1.0-alpha[k][i])*curf->c[i][j]+alpha[k][i]*view[j];*/
	  cproj[k][i][j] = curf->c[i][j] + alpha[k][i]*tvec[j];
    }
  }

#if 1 == 0
  /* project origin and (1 0 0) into view plane */
  /* live dangerously---no perpendicular check */
  temp = dot(normal, view);
  if(temp == 0.0) {
    fprintf(stderr, 
	"is1stFaceDeeper: origin-view vector perpendicular to view plane\n");
  }
  /* get projection alpha */
  alpha_origin = rhs/temp;
  for(j = 0; j < 3; j++) origin[j] = alpha_origin*view[j];

  /* project (1 0 0) */
  for(j = 0; j < 3; j++) tvec[j] = view[j];
  tvec[0] -= 1.0;
  if((temp = dot(normal, tvec)) == 0.0) {
    fprintf(stderr,
	    "is1stFaceDeeper: x-view vector perpendicular to view plane\n");
  }
  /* get projection alpha */
  tvec1[0] = 1.0; tvec1[1] = tvec1[2] = 0.0;
  alpha_x = (rhs - dot(normal, tvec1))/temp;
  /* get x-axis projection */
  for(j = 0; j < 3; j++) x[j] = tvec1[j] + alpha_x*tvec[j] - origin[j];
  temp = sqrt(dot(x, x));
  for(j = 0; j < 3; j++) x[j] /= temp;
  /* y = zXx */
  crossProd(y, normal, x);
#endif

#if 1 == 1
  /* figure x and y coordinates in view plane (normal is z coordinate) */
  /* x = c0-c1 projections from fac proj. always (should never be 0 len) */
  for(j = 0; j < 3; j++) origin[j] = cproj[0][0][j]; /* to stop overwrites */
  for(j = 0; j < 3; j++) x[j] = cproj[0][1][j] - origin[j];
  temp = sqrt(dot(x, x));
  for(j = 0; j < 3; j++) x[j] /= temp; /* normalize */
  /* y = zXx */
  crossProd(y, normal, x);
#endif
  
  /* project all face corner projections onto new x and y coordinates
     - cproj[0][0] plays the role of origin */
  for(curf = fac, k = 0; k < 2; k++, curf = facref) {	/* loop on faces */
    for(i = 0; i < curf->numsides; i++) {
      for(j = 0; j < 3; j++) tvec1[j] = cproj[k][i][j] - origin[j];
      tvec[0] = dot(x, tvec1);	/* get weight in x direction */
      tvec[1] = dot(y, tvec1);	/* get weight in y direction */
      tvec[2] = 0.0;		/* all z weights must = rhs */
      for(j = 0; j < 3; j++) cproj[k][i][j] = tvec[j]; /* xfer */
    }
  }
  for(j = 0; j < 3; j++) cproj[0][0][j] = 0.0; /* set origin explicitly zero */

  /* for each side of facref, see if there is an intersect. w/ all fac sides */
#if DEBUGX == ON
  fprintf(stdout, "Is face %d behind face %d?\n", fac->index, facref->index);
#endif
  is_overlap = FALSE;
  isect_cnt = 0;
  isect_avg[0] = isect_avg[1] = isect_avg[2] = 0.0;
  for(i = 0; i < facref->numsides; i++) {
    cfr = cproj[1][i];
    if(i == facref->numsides - 1) ctr = cproj[1][0];
    else ctr = cproj[1][i+1];
    for(j = 0; j < fac->numsides; j++) {
      cff = cproj[0][j];
      if(j == fac->numsides - 1) ctf = cproj[0][0];
      else ctf = cproj[0][j+1];
      if((intersect = doLinesIntersect(isect, cfr, ctr, cff, ctf)) == TRUE) {
	isect_cnt++;
	for(k = 0; k < 3; k++) isect_avg[k] += isect[k];
	is_overlap = TRUE;
      }
#if DEBUGX == ON
      fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
    }
  }
#if XOVTST == ON		/* do either this or face_is_inside() below */
  /* check for overlap with lines across face */
  if(facref->numsides == 4) {
    for(i = 0; i < 2; i++) {
      if(i == 0) {
	cfr = cproj[1][0];
	ctr = cproj[1][2];
      }
      else {
	cfr = cproj[1][1];
	ctr = cproj[1][3];
      }
      for(j = 0; j < fac->numsides; j++) {
	cff = cproj[0][j];
	if(j == fac->numsides - 1) ctf = cproj[0][0];
	else ctf = cproj[0][j+1];
	if((intersect = doLinesIntersect(cfr, ctr, cff, ctf)) == TRUE)
	    is_overlap = TRUE;
#if DEBUGX == ON
	fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
      }
    }
  }
  else if(facref->numsides == 3) {
    avg[0] = (cproj[1][0][0] + cproj[1][1][0] + cproj[1][2][0])/3.0;
    avg[1] = (cproj[1][0][1] + cproj[1][1][1] + cproj[1][2][1])/3.0;
    avg[2] = (cproj[1][0][2] + cproj[1][1][2] + cproj[1][2][2])/3.0;
    cfr = avg;
    for(i = 0; i < 3; i++) {
      ctr = cproj[1][i];
      for(j = 0; j < fac->numsides; j++) {
	cff = cproj[0][j];
	if(j == fac->numsides - 1) ctf = cproj[0][0];
	else ctf = cproj[0][j+1];
	if((intersect = doLinesIntersect(isect, cfr, ctr, cff, ctf)) == TRUE)
	    is_overlap = TRUE;
#if DEBUGX == ON
	fprintf(stdout, "doLinesIntersect returned %d\n", intersect);
#endif
      }
    }
  }
  else {
    fprintf(stderr, "isThereProjOverlap: can't handle %d side panel\n",
	    facref->numsides);
    exit(0);
  }
#else
  /* no sides overlap---check if one face completely obscures the other */
  if(!is_overlap) {
    is_overlap 
	= face_is_inside(cproj[1], facref->numsides, cproj[0], fac->numsides,
			 isect_avg);
  }
  else {
    /* figure average intersect point---should be inside overlap region */
    for(k = 0; k < 3; k++) isect_avg[k] /= (double)isect_cnt;
    if(isect_cnt % 2 != 0) {
      k = k + 0;
    }
  }
#endif				/* XOVTST == ON */

  /* if no overlap, no edge in graph in any case */
  if(!is_overlap) return(FALSE);

  /* return TRUE only if fac is deeper than facref */
  /* return REVERSE only if facref is deeper than fac */
  /* project average point back along view direction to fac, facref planes */
  /*   first convert average point back to 3D */
  isect[0] = isect_avg[0]*x[0] + isect_avg[1]*y[0];
  isect[1] = isect_avg[0]*x[1] + isect_avg[1]*y[1];
  isect[2] = isect_avg[0]*x[2] + isect_avg[1]*y[2];
  /*   add in origin to get point in absolute coordinates */
  for(k = 0; k < 3; k++) isect[k] += origin[k];
  /*   figure projection back to fac */
  for(k = 0; k < 3; k++) isect[k] -= view[k];
  alpha_fac = (fac->rhs - dot(fac->normal, view))/dot(fac->normal, isect);
  alpha_facref 
      = (facref->rhs - dot(facref->normal, view))/dot(facref->normal, isect);

  /* a larger alpha => greater distance from view plane */
  /* => fac is deeper than facref if alpha_fac > alpha_facref */
  /*     (ie if alpha_facref - alpha_fac < 0) */
  /* => facref is deeper than fac if alpha_facref > alpha_fac */
  /* => inconclusive if equal */
  if(diff_is_negative(alpha_facref, alpha_fac, 1.0)) return(TRUE);
  else if(diff_is_negative(alpha_fac, alpha_facref, 1.0)) return(REVERSE);
  else {
    fprintf(stderr, 
	    "\nis1stFaceDeeper: Warning, face ordering test failure\n");
    fprintf(stderr, "  alpha_fac, face %d = %g alpha_facref, face %d = %g\n", 
	    fac->index, alpha_fac, facref->index, alpha_facref);
    dump_face(stderr, fac);
    fprintf(stderr, " Projected corners\n");
    dumpCorners(stderr, cproj[0], fac->numsides, 3);
    dump_face(stderr, facref);
    fprintf(stderr, " Projected corners\n");
    dumpCorners(stderr, cproj[1], facref->numsides, 3);

    return(FALSE);		/* inconclusive test */
  }

}
#endif

/*
  returns TRUE if bounding box of facref and fac (proj to facref's plane) insct
*/
int isThereBoxOverlap(fac, facref, view)
face *fac, *facref;
double *view;
{
  int i, j, olap[2];
  double cproj[MAXSIDES][3];	/* corners of fac in facref's plane */
  double cref[MAXSIDES][3];	/* corners of facref in facref plane */
  double alpha[MAXSIDES];	/* 1 => view point 0 => corner */
  double minref[2], maxref[2];	/* bounding box coordinates */
  double minfac[2], maxfac[2];
  double x[3], y[3];		/* coordinates of x and y in facref plane */
  double dot(), temp, tvec[3], tvec1[3], margin, ovrlapmgn = 0.0;

  /* figure projections of fac's corners back to facref's plane rel to view */
  for(i = 0; i < fac->numsides; i++) {
    for(j = 0; j < 3; j++) tvec[j] = view[j] - fac->c[i][j]; /* get v-c */
    temp = dot(facref->normal, tvec); /* get n.(v-c) */
    margin = sqrt(dot(tvec, tvec))*MARGIN;
    ovrlapmgn = MAX(margin, ovrlapmgn);	/* used below */
    /* test fails if v-c is perpendicular to n */
    if(temp > -margin && temp < margin) return(FALSE);
    /* get alpha as in c + alpha*(v-c) = c' */
    alpha[i] = (facref->rhs - dot(facref->normal, fac->c[i]))/temp;
    /* if(alpha[i] < -margin || alpha[i] > 1.0+margin) {
      fprintf(stderr, "\nisThereBoxOverlap: big X failure, alpha = %g\n",
	      alpha[i]);
      exit(0);
    } */
    for(j = 0; j < 3; j++)	/* get c' */
	cproj[i][j] = (1.0-alpha[i])*fac->c[i][j]+alpha[i]*view[j];
  }

  /* figure x and y coordinates in facref plane (normal is z coordinate) */
  /* x = c0-c1 always */
  for(j = 0; j < 3; j++) x[j] = facref->c[0][j] - facref->c[1][j];
  temp = sqrt(dot(x, x));
  for(j = 0; j < 3; j++) x[j] /= temp; /* normalize */
  /* y = zXx */
  crossProd(y, facref->normal, x);
  
  /* project all fac corner projections onto new x and y coordinates
     - facref->c[0] plays the role of origin */
  for(i = 0; i < fac->numsides; i++) {
    for(j = 0; j < 3; j++) tvec1[j] = cproj[i][j] - facref->c[0][j];
    tvec[0] = dot(x, tvec1);	/* get weight in x direction */
    tvec[1] = dot(y, tvec1);	/* get weight in y direction */
    tvec[2] = 0.0;		/* all z weights must = facref->rhs */
    for(j = 0; j < 3; j++) cproj[i][j] = tvec[j]; /* xfer */
  }
  for(j = 0; j < 3; j++) cref[0][j] = 0.0;
  for(i = 1; i < facref->numsides; i++) {
    for(j = 0; j < 3; j++) tvec1[j] = facref->c[i][j] - facref->c[0][j];
    cref[i][0] = dot(x, tvec1);	/* get weight in x direction */
    cref[i][1] = dot(y, tvec1);	/* get weight in y direction */
    cref[i][2] = 0.0;		/* all z weights must = facref->rhs */
  }

  /* figure bounding boxes in new coordinates */
  minfac[0] = maxfac[0] = cproj[0][0];
  minfac[1] = maxfac[1] = cproj[0][1];
  for(i = 1; i < fac->numsides; i++) { /* find max, min for fac */
    for(j = 0; j < 2; j++) {
      minfac[j] = MIN(minfac[j], cproj[i][j]);
      maxfac[j] = MAX(maxfac[j], cproj[i][j]);
    }
  }
  minref[0] = maxref[0] = cref[0][0];
  minref[1] = maxref[1] = cref[0][1];
  for(i = 1; i < facref->numsides; i++) { /* find max, min for facref */
    for(j = 0; j < 2; j++) {
      minref[j] = MIN(minref[j], cref[i][j]);
      maxref[j] = MAX(maxref[j], cref[i][j]);
    }
  }

  /* check for overlap - call things overlaped even if not to be sure
     - no overlap if either no overlap in x or none in y */
  /* check for x overlap */
  olap[0] = olap[1] = FALSE;
  for(j = 0; j < 2; j++) {
    if((minref[j]-ovrlapmgn < minfac[j] && minfac[j] < maxref[j]+ovrlapmgn) ||
       (minref[j]-ovrlapmgn < maxfac[j] && maxfac[j] < maxref[j]+ovrlapmgn) ||
       (minfac[j]-ovrlapmgn < minref[j] && minref[j] < maxfac[j]+ovrlapmgn) ||
       (minfac[j]-ovrlapmgn < maxref[j] && maxref[j] < maxfac[j]+ovrlapmgn))
	olap[j] = TRUE;
  }
  if(olap[0] == FALSE || olap[1] == FALSE) return(FALSE);
  else return(TRUE);
}

#if 1 == 0
/*
  returns TRUE if 1st face is deeper than second - uses Don's 3d X method
  - returns FALSE in don't care situations (1st need not be ordered before 2nd)
*/
int is1stFaceDeeper(first, second, view)
face *first, *second;
double *view;
{
  int fst, snd, ret = SAME;

  /* figure on which side of each face plane lies the other face */
  fst = whichSide(first, second); /* first relative to plane of second */
  snd = whichSide(second, first); /* second relative to plane of first */

  /* figure returned value */
  if((fst == NEG && snd == NEG) || (fst == POS && snd == POS) ||
      (fst == SAME && snd == SAME))
      ret = FALSE;		/* order arbitrary */
  else if((fst == NEG && snd == POS)) ret = TRUE;
  else if((fst == POS && snd == NEG)) ret = FALSE;
  else if(fst == NEG && snd == SPLIT) ret = TRUE;
  else if(fst == POS && snd == SPLIT) ret = FALSE;
  else if(fst == SPLIT && snd == POS) ret = TRUE;
  else if(fst == SPLIT && snd == NEG) ret = FALSE;
  /* next four added 22 Aug 91 w/o really thinking
     - catches false SAME's ? */
/*  else if(fst == NEG && snd == SAME) ret = TRUE;
  else if(fst == POS && snd == SAME) ret = FALSE;
  else if(fst == SAME && snd == POS) ret = TRUE;
  else if(fst == SAME && snd == NEG) ret = FALSE; */
  else if(snd == SPLIT && fst == SPLIT) {
    return(isThereProjOverlap(first, second, view));
#if 1 == 0
    /* if both are split, see if they interact at all */
    if(isThereProjOverlap(first, second, view) == FALSE) ret = FALSE;
    else {
      fprintf(stderr, 
	      "\nis1stFaceDeeper: unresolvable both faces split case\n");
      ret = TRUE;		/* don't order--try to order w/other check
				   sometimes works, sometimes not */
    }
#endif
  }
  if(ret == SAME) {
    fprintf(stderr,"\nis1stFaceDeeper: bad fst = %d and snd = %d\n", fst, snd);
    exit(0);
  }
  else if(ret == TRUE) {	/* make sure 2nd could overlap 1st */
    if(isThereProjOverlap(first, second, view) == FALSE) ret = FALSE;
  }
  return(ret);
}
/*
#else

int is1stFaceDeeper(first, second, view)
face *first, *second;
double *view;
{
  return(isThereProjOverlap(first, second, view));
}
*/
#endif

/*
  recursive guts of below
*/
int chkCycle(fac, ref, fp)
face *fac, *ref;
FILE *fp;
{
  int b, i;

  if(fac->mark == TRUE) return(FALSE);

  fac->mark = TRUE;

  if(fac->numbehind == 0) return(FALSE);
  else {
    for(b = 0; b < fac->numbehind; b++) {
      /*fprintf(fp, " %d (%d)", (fac->behind)[b]->depth,
	      (fac->behind)[b]->index);
      if(b % 5 == 0 && b != 0) fprintf(fp, "\n");*/
      if(fac->behind[b] == ref) return(TRUE);
      else if(chkCycle(fac->behind[b], ref, fp) == TRUE) return(TRUE);
    }
/*    if((i-1) % 5 != 0 || i == 1) fprintf(fp, "\n"); */
  }
  return(FALSE);
}

/*
  checks for cycles in the depth graph - BROKEN(?)
*/
void dumpCycles(faces, numfaces, file)
face **faces;
int numfaces;
FILE *file;
{
  int i, f, j, b, *cycled, *cyclei, numcycle, cycle = FALSE;
  face *fp;

  /* for each face, chase behind pointers until a leaf or same face is found */
  /*fprintf(file, "\nRecursive behind lists\n");*/
  for(f = 0; f < numfaces; f++) {
/*    fprintf(file, "%d (%d):", faces[f]->depth, faces[f]->index);*/
    for(j = 0; j < numfaces; j++) faces[j]->mark = FALSE;
    for(b = 0; b < faces[f]->numbehind; b++) {
      if(chkCycle(faces[f]->behind[b], faces[f], file) == TRUE) {
	cycle = TRUE;
	break;
      }
    }
    if(cycle == TRUE) break;
  }
  if(cycle == FALSE) fprintf(file, "Adjacency graph has no cycles\n");
  else fprintf(file, "Adjacency graph has cycles\n");
  for(j = 0; j < numfaces; j++) faces[j]->mark = FALSE;
}

    
/*
  recursively sets depths of faces
*/
void setDepth(fac)
face *fac;
{
  int i;

  /* mark so this face won't be renumbered */
  fac->mark = TRUE;

/*
  if(fac->index == 193) {
    i = i + 0;
  }
*/

  /* do adjacents if needed */
  for(i = 0; i < fac->numbehind; i++) {
    if((fac->behind[i])->mark == FALSE) setDepth(fac->behind[i]);
  }

  if(fac->index == 131 || fac->index == 193) {
    i = i + 0;
  }

  /* set depth, update counter */
  fac->depth = cnt--;
}

/*
  does a topological sorting of the faces using graph setup by getAdjGraph()
  - returns a new set of pointers with deepest (1st to print) face first
*/
face **depthSortFaces(faces, numfaces)
face **faces;
int numfaces;
{
  int f, i, facefound;
  face **rfaces;

  /* make sure all marks are cleared */
  for(f = 0; f < numfaces; f++) faces[f]->mark = FALSE;

  /* set depths recursively - zero is deepest (first to be rendered) */
  for(f = 0, cnt = numfaces-1; f < numfaces; f++) {
    if(faces[f]->mark == FALSE) setDepth(faces[f]);
  }

  /* make the new set of pointers */
  CALLOC(rfaces, numfaces, face *, ON, AMSC);
  for(f = 0; f < numfaces; f++) {
    for(i = 0, facefound = FALSE; i < numfaces; i++) {
      if(faces[i]->depth == f) {
	rfaces[f] = faces[i];
	facefound = TRUE;
	break;
      }
    }
    if(facefound == FALSE) {
      fprintf(stderr, "\ndepthSortFaces: can't find depth %d face\n", f);
      exit(0);
    }
  }
  return(rfaces);
}

/*
  sets up adjacency graph pointers in faces: pntr in i to j => face i behind j
*/
void getAdjGraph(faces, numfaces, view, rhs, normal)
face **faces;			/* array of face pntrs, faces[0] head of lst */
int numfaces;
double *view, rhs, *normal;
{
  int numbehind, f, i, check;
  face *fpcur, *fpchk, **behind;
  extern int k_;

  if(TRUE == TRUE) {		/* if memory can be sacrificed for speed */
    /* set up huge n^2 blocked face pointer arrays for each face */
    for(f = 0; f < numfaces; f++) {
      CALLOC(faces[f]->behind, numfaces, face *, ON, AMSC);
      faces[f]->numbehind = 0;
    }      

    /* for each face, check through all faces not previously checked */
    for(fpcur = faces[0], f = 0; fpcur != NULL; fpcur = fpcur->next, f++) {
      for(fpchk = fpcur->next, numbehind = i = 0; fpchk != NULL; 
	  fpchk = fpchk->next, i++) {
	if(fpchk == fpcur) continue;	/* a face can't be behind itself */
	if((check = is1stFaceDeeper(fpcur, fpchk, view, rhs, normal))==TRUE) { 
	  fpcur->behind[(fpcur->numbehind)++] = fpchk;
	}
	else if(check == REVERSE) fpchk->behind[(fpchk->numbehind)++] = fpcur;
      }
      if(f % 20 == 0 && f != 0) {
	fprintf(stdout, "%d ", f);
	fflush(stdout);
      }
      if(f % 200 == 0 && f != 0) fprintf(stdout, "\n");
    }
  }
  else {
    CALLOC(behind, numfaces, face *, ON, AMSC);

    /* for each face, check through all the other faces */
    for(fpcur = faces[0], f = 0; fpcur != NULL; fpcur = fpcur->next, f++) {
      for(fpchk = faces[0], numbehind = i = 0; fpchk != NULL; 
	  fpchk = fpchk->next, i++) {
	if(fpchk == fpcur) continue;	/* a face can't be behind itself */
	if(is1stFaceDeeper(fpcur, fpchk, view) == TRUE) { 
	  behind[numbehind++] = fpchk;
	}
      }
      /* transfer blocked face pointers to face struct */
      fpcur->numbehind = numbehind;
      if(numbehind > 0) CALLOC(fpcur->behind, numbehind, face *, ON, AMSC);
      for(; numbehind > 0; numbehind--)
          fpcur->behind[numbehind-1] = behind[numbehind-1];
      if(f % 20 == 0 && f != 0) {
	fprintf(stdout, "%d ", f);
	fflush(stdout);
      }
      if(f % 200 == 0 && f != 0) fprintf(stdout, "\n");
    }
  }
}
