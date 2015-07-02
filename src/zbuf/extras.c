#include "mulGlobal.h"
#include "zbufGlobal.h"

#define XI 0
#define YI 1
#define ZI 2
#define EQUIV_TOL 1.0e-9
/*#define PLANAR_TOL 1.0e-3*/
#define PLANAR_TOL 1.0e-5
#define MAX_REAL 1.0e+20;

#define Dot_Product(V1,V2) (V1[XI]*V2[XI]+V1[YI]*V2[YI]+V1[ZI]*V2[ZI])
#define DotP_Product(V1,R,S,T) ((V1[XI])*(R)+(V1[YI])*(S)+(V1[ZI])*(T))

#define ONE3 0.3333333333333

/*
  dumps a rows x cols matrix of doubles; assumes indices from zero 
*/
void dumpCorners(fp, mat, rows, cols)
int rows, cols;
double **mat;
FILE *fp;
{
  int i, j;
  for(i = 0; i < rows; i++) {
    fprintf(fp, "  corner%d ", i);
    for(j = 0; j < cols; j++) {
      if(mat[i][j] < 0.0) fprintf(fp, "%.5e ", mat[i][j]);
      else fprintf(fp, " %.5e ", mat[i][j]);
    }
    fprintf(fp, "\n");
  }
}

/* 
  dumps state of important compile flags
*/
void dumpConfig(fp, name)
char *name;
FILE *fp;
{
  int size = -1;		/* for '#define MAXITER size' case */

  fprintf(fp, "\n%s CONFIGURATION FLAGS:\n", name);

  fprintf(fp, " DISCRETIZATION CONFIGURATION\n");

  fprintf(fp, "   WRMETH");
  if(WRMETH == COLLOC)
      fprintf(fp, " == COLLOC (point collocation)\n");
  else if(WRMETH == SUBDOM)
      fprintf(fp, " == SUBDOM (not implemented - do collocation)\n");
  else if(WRMETH == GALKIN)
      fprintf(fp, " == GALKIN (not implemented - do collocation)\n");
  fprintf(fp, "   ELTYPE");
  if(ELTYPE == CONST)
      fprintf(fp, " == CONST (constant panel densities)\n");
  else if(ELTYPE == AFFINE)
      fprintf(fp, " == AFFINE (not implemented - use constant)\n");
  else if(ELTYPE == QUADRA)
      fprintf(fp, " == QUADRA (not implemented - use constant)\n");

  fprintf(fp, " MULTIPOLE CONFIGURATION\n");

  fprintf(fp, "   DNTYPE");
  if(DNTYPE == NOLOCL) 
      fprintf(fp, " == NOLOCL (no locals in dwnwd pass)\n");
  else if(DNTYPE == NOSHFT) 
      fprintf(fp, " == NOSHFT (no local2local shift dwnwd pass)\n");
  else if(DNTYPE == GRENGD) 
      fprintf(fp, " == GRENGD (full Greengard dwnwd pass)\n");
  fprintf(fp, "   MULTI");
  if(MULTI == ON) fprintf(fp, " == ON (include multipole part of P*q)\n");
  else fprintf(fp, " == OFF (don't use multipole part of P*q)\n");
  fprintf(fp, "   RADINTER");
  if(RADINTER == ON) 
      fprintf(fp," == ON (allow parent level interaction list entries)\n");
  else 
   fprintf(fp," == OFF (use only cube level interaction list entries)\n");
  fprintf(fp, "   NNBRS == %d (max distance to a nrst neighbor)\n", NNBRS);
  fprintf(fp, "   ADAPT");
  if(ADAPT == ON) 
      fprintf(fp, " == ON (adaptive - no expansions in exact cubes)\n");
  else fprintf(fp, " == OFF (not adaptive - expansions in all cubes)\n");
  fprintf(fp, "   OPCNT");
  if(OPCNT == ON) 
      fprintf(fp, " == ON (count P*q ops - exit after mat build)\n");
  else fprintf(fp, " == OFF (no P*q op count - iterate to convergence)\n");

  fprintf(fp, "   MAXDEP");
  fprintf(fp, 
	  " == %d (assume no more than %d partitioning levels are needed)\n",
	  MAXDEP, MAXDEP);

  fprintf(fp, "   NUMDPT");
  fprintf(fp, 
	  " == %d (do %d potential evaluations for each dielectric panel)\n",
	  NUMDPT, NUMDPT);

  fprintf(fp, " LINEAR SYSTEM SOLUTION CONFIGURATION\n");

  fprintf(fp, "   ITRTYP");
  if(ITRTYP == GCR)
      fprintf(fp, " == GCR (generalized conjugate residuals)\n");
  else if(ITRTYP == GMRES)
      fprintf(fp, " == GMRES (generalized minimum residuals)\n");
  else fprintf(fp, " == %d (not implemented - use GCR)\n", ITRTYP);

  fprintf(fp, "   PRECOND");
  if(PRECOND == BD) {
    fprintf(fp, 
	    " == BD (use block diagonal preconditioner)\n");
  }
  else if(PRECOND == OL) {
    fprintf(fp, 
	    " == OL (use overlap preconditioner)\n");
  }
  else if(PRECOND == NONE) {
    fprintf(fp, 
	    " == NONE (no preconditioner)\n");
  }
  else fprintf(fp, " == %d (not implemented - use BD, OL or NONE)\n", PRECOND);

  fprintf(fp, "   DIRSOL");
  if(DIRSOL == ON) 
      fprintf(fp, " == ON (do the whole calculation directly)\n");
  else fprintf(fp, " == OFF (do the calculation iteratively)\n");

  fprintf(fp, "   EXPGCR");
  if(EXPGCR == ON) 
      fprintf(fp, " == ON (do all P*q's explicitly w/full matrix)\n");
  else fprintf(fp, " == OFF (do P*q's with multipole)\n");

  fprintf(fp, "   MAXITER");
  if(MAXITER < 0) {
    fprintf(fp, " == size (for n panel system, do at most n iterations)\n");
  }
  else fprintf(fp, " == %d (stop after %d iterations if not converged)\n", 
	  MAXITER, MAXITER);

  fprintf(fp, "   EXRTSH");
  fprintf(fp, 
	  " == %g (exact/ttl cubes > %g on lowest level => stop refinement)\n",
	  EXRTSH, EXRTSH);
}

/*
  dump the contents of a face struct
*/
void dump_face(fp, fac)
face *fac;
FILE *fp;
{
  int i, j;
  face **behind = fac->behind;

  fprintf(fp, "Face %d, %d sides, depth %d, mark %d, greylev %d\n", 
	  fac->index, fac->numsides, fac->depth, fac->mark, fac->greylev);
  fprintf(fp, "  plane: n = (%g %g %g) rhs = %g\n",
	  fac->normal[0], fac->normal[1], fac->normal[2], fac->rhs);
  fprintf(fp, "  behind [depth(index)]:");
  for(i = 0; i < fac->numbehind; i++) {
    fprintf(fp, " %d(%d)", behind[i]->depth, behind[i]->index);
    if(i % 10 == 0 && i != 0) fprintf(fp, "\n");
  }
  i--;
  if(!(i % 10 && i != 0)) fprintf(fp, "\n");
  fprintf(fp, " Corners\n");
  dumpCorners(fp, fac->c, fac->numsides, 3);
}  

initcalcp(panel_list)
  charge *panel_list;
{
  charge *pq, *npq;
  double vtemp[3];
  double length, maxlength, minlength, length20, length31, sum, sum2, delta;
  double normalize();
  int i, j, next;

  for(i=0, pq = panel_list; pq != NULL; pq = pq->next) i++;

  pq = panel_list; 
  while (pq != NULL) {
    /* Calculate edge lengths. */
    maxlength = 0.0;
    minlength = MAX_REAL;
    for(i=0; i < pq->shape; i++) {    
      if(i == (pq->shape -1)) next = 0;
      else next = i + 1;
      for(sum= 0, j = 0; j < 3; j++) {
	delta = pq->corner[next][j] - pq->corner[i][j];
	sum += delta * delta;
      }
      pq->length[i] = length = sqrt(sum);
      maxlength = MAX(maxlength, length);
      minlength = MIN(minlength, length);
    }
    
    /* Get diags and lengths. */
    for(sum= 0, sum2 = 0, i = 0; i < 3; i++) {     
      pq->X[i] = delta = pq->corner[2][i] - pq->corner[0][i];
      sum += delta * delta;
      if(pq->shape == 3) pq->Y[i] = pq->corner[0][i] - pq->corner[1][i];      
      else {
	pq->Y[i] = delta = pq->corner[3][i] - pq->corner[1][i];      
	sum2 += delta * delta;
      }
    }
    length20 = sqrt(sum);
    length31 = sqrt(sum2);

    /* Check on lengths for quad. */
    if(pq->shape == 3) {
      pq->max_diag = maxlength;
      pq->min_diag = minlength;
    }
    else {
      length = MAX(length20, length31);
      pq->max_diag = length;
      pq->min_diag = MIN(length20, length31);
    }

    /* Z-axis is normal to two diags. */
    Cross_Product(pq->X, pq->Y, pq->Z);
/*#if 1 == 0*/
    if(flip_normal(pq)) {
      for(i = 0; i < 3; i++) {
	pq->Z[i] = -(pq->Z[i]);	/* flip the normal */
	pq->X[i] = -(pq->X[i]);	/* flip the x direction */
	/* interchange points 0 and 2 so that corner order will be
	   consistent with X flip (note this is OK for both quads and tris) */
	vtemp[0] = pq->corner[0][i];
	pq->corner[0][i] = pq->corner[2][i];
	pq->corner[2][i] = vtemp[0];
      }
      /* interchange length entries in length vector */
      vtemp[0] = pq->length[0];
      pq->length[0] = pq->length[1];
      pq->length[1] = vtemp[0];
      if(pq->shape == 4) {
	vtemp[0] = pq->length[2];
	pq->length[2] = pq->length[3];
	pq->length[3] = vtemp[0];
      }
    } 
/*#endif*/
    pq->area = 0.5 * normalize(pq->Z);
    normalize(pq->X);

    /* Real Y-axis is normal to X and Z (resulting system is left-handed). */
    Cross_Product(pq->X, pq->Z, pq->Y);

    /* Project the corner points into the plane defined by edge midpoints. */
    if(planarize(pq) == FALSE) {     
      /* Planarization too drastic, crack into two triangles. */
      CALLOC(npq, 1, charge, ON, AMSC);
      npq->next = pq->next;
      pq->next = npq;
      npq->shape = 3;
      pq->shape = 3;
      npq->cond = pq->cond;
      npq->surf = pq->surf;
      VCOPY(npq->corner[0], pq->corner[0]);
      VCOPY(npq->corner[1], pq->corner[2]);
      VCOPY(npq->corner[2], pq->corner[3]);
      continue;
    }

    /* Calculate the centroid. */
    centroid(pq, length20);      

    /* Put corners in the newly defined coord system. */
    for(i=0; i < pq->shape; i++) {
      pq->corner[i][XI] -= pq->x;
      pq->corner[i][YI] -= pq->y;
      pq->corner[i][ZI] -= pq->z;
    }
    for(i=0; i < pq->shape; i++) {
      vtemp[XI] = Dot_Product(pq->corner[i], pq->X);
      vtemp[YI] = Dot_Product(pq->corner[i], pq->Y);
      vtemp[ZI] = Dot_Product(pq->corner[i], pq->Z);
      pq->corner[i][XI] = vtemp[XI];
      pq->corner[i][YI] = vtemp[YI];
      pq->corner[i][ZI] = vtemp[ZI];
      if(fabs(pq->corner[i][ZI]) > (EQUIV_TOL * pq->min_diag)) {
	printf("FATAL PROGRAM ERROR: renormalized z=%g\n", pq->corner[i][ZI]);
	exit(0);
      }
      pq->corner[i][ZI] = 0.0;
    }

    /* Iterate for the next panel. */
    pq = pq->next;

  }

}

/* Calculates result_vector = vector1 X vector2. */
Cross_Product(vector1, vector2, result_vector)
  double vector1[], vector2[], result_vector[];
{
  result_vector[XI] = vector1[YI]*vector2[ZI] - vector1[ZI]*vector2[YI];
  result_vector[YI] = vector1[ZI]*vector2[XI] - vector1[XI]*vector2[ZI];
  result_vector[ZI] = vector1[XI]*vector2[YI] - vector1[YI]*vector2[XI];
}

double normalize(vector)
  double vector[3];
{
  double length;
  int i;

  length = sqrt( vector[0]*vector[0] 
		+ vector[1]*vector[1] 
		+ vector[2]*vector[2]);
    
  for (i=0; i<3; i++) vector[i] = vector[i] / length;

  return length;
}

/* 
Determines centroid of a panel (defined as the point which make the
first moments vanish.  Calculation begins by projection into the
coordinate system defined by the panel normal as the z-axis and
edge02 as the x-axis.
*/
centroid(pp, x2)
charge *pp;
double x2;
{
  double vertex1[3], vertex3[3];
  double sum, dl, x1, y1, x3, y3, xc, yc;
  int i;

  /* Use vertex 0 as the origin. */
  for(i=0; i< 3; i++) {
    vertex1[i] = pp->corner[1][i] - pp->corner[0][i];
    if(pp->shape == 4) vertex3[i] = pp->corner[3][i] - pp->corner[0][i];
    else vertex3[i] = pp->corner[2][i] - pp->corner[0][i];
  }

  /* Project into the panel axes. */
  y1 = Dot_Product(vertex1, pp->Y);
  y3 = Dot_Product(vertex3, pp->Y);
  x1 = Dot_Product(vertex1, pp->X);
  x3 = Dot_Product(vertex3, pp->X);

  yc = ONE3 * (y1 + y3);
  xc = ONE3 * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

  pp->x = pp->corner[0][XI] + xc * pp->X[XI] + yc * pp->Y[XI];
  pp->y = pp->corner[0][YI] + xc * pp->X[YI] + yc * pp->Y[YI];
  pp->z = pp->corner[0][ZI] + xc * pp->X[ZI] + yc * pp->Y[ZI];

}

/*
Changes the corner points so that they lie in the plane defined by the
panel diagonals and any midpoint of an edge.
*/
planarize(pq)
charge *pq;
{
  double origin[3], corner[3], delta[4][3], px, py, dx, dy, dz;
  int i, j, numcorners = pq->shape;
  double tolsq = PLANAR_TOL * pq->min_diag; 

  tolsq *= tolsq;

  /* Triangular panels are just fine already. */
  if(numcorners != 4) return(TRUE);

  /* Pick edge midpoint as origin. */
  for(i=0; i < 3; i++) origin[i] = 0.5 * (pq->corner[1][i] + pq->corner[0][i]);

  for(i=0; i < numcorners; i++) {
    for(j=0; j < 3; j++) corner[j] = pq->corner[i][j] - origin[j];
    px = Dot_Product(corner, pq->X);
    py = Dot_Product(corner, pq->Y);

    dx = px * pq->X[XI] + py * pq->Y[XI] + origin[XI] - pq->corner[i][XI];
    dy = px * pq->X[YI] + py * pq->Y[YI] + origin[YI] - pq->corner[i][YI];
    dz = px * pq->X[ZI] + py * pq->Y[ZI] + origin[ZI] - pq->corner[i][ZI];

    delta[i][XI] = dx;
    delta[i][YI] = dy;
    delta[i][ZI] = dz;
    
    /* Moved too much, crack the panel. */
    if((dx * dx + dy * dy + dz * dz) > tolsq) return(FALSE);
  }

  /* Still Here? Must be planarizing. */
  for(i=0; i < numcorners; i++) {
    for(j=0; j < 3; j++) {
      pq->corner[i][j] += delta[i][j];
    }
  }
  return(TRUE);
}

/*
  determine if normal needs to be flipped to get dielectric bdry cond right
  - this function uses 0.0 as a breakpoint when really machine precision
    weighted checks should be done (really not an issue if ref point far)
*/
flip_normal(panel)
charge *panel;
{
  int i;
  double x, y, z;
  double ctr_minus_n[3], ctr_plus_n[3], norm_minus, norm_plus, norm, norm_sq;
  surface *surf = panel->surf;
  int ref_inside = surf->ref_inside, flip_normal;
  double *ref = surf->ref, *normal, angle, norm_n;
  char *surf_name = surf->name;

  if(surf->type != DIELEC && surf->type != BOTH) return(FALSE);

  /* get panel corner (relative to reference point) and normal */
  x = panel->corner[0][0] - ref[0]; 
  y = panel->corner[0][1] - ref[1]; 
  z = panel->corner[0][2] - ref[2];
  norm_sq = x*x + y*y + z*z;
  norm = sqrt(norm_sq);
  normal = panel->Z;
  norm_n = 
      sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

  /* figure the dot product between the normal and the corner-ref point
     - ref_inside and angle <= 90 => flip
     - ref_inside and angle  > 90 => no flip
     - !ref_inside and angle <= 90 => no flip 
     - !ref_inside and angle > 90 => flip */
  /*    figure angle */
  angle = (x*normal[0] + y*normal[1] + z*normal[2])/(norm*norm_n);
  if(ref_inside && angle <= 0.0) flip_normal = TRUE;
  else if(ref_inside && angle > 0.0) flip_normal = FALSE;
  else if(!ref_inside && angle <= 0.0) flip_normal = FALSE;
  else if(!ref_inside && angle > 0.0) flip_normal = TRUE;
  else {
    fprintf(stderr, 
	    "flip_normal: inconclusive test for normal flipping\n");
    fprintf(stderr, "  Surface: %s\n", hack_path(surf_name));
    fprintf(stderr, "  Translation: (%g %g %g)\n", surf->trans[0],
	    surf->trans[1], surf->trans[2]);
    fprintf(stderr, "  Reference point: (%g %g %g)\n",
	    ref[0], ref[1], ref[2]);
    fprintf(stderr, "  Panel corner: (%g %g %g)\n",
	    panel->corner[0][0], panel->corner[0][1], panel->corner[0][2]);
    fprintf(stderr, "  Normal: (%g %g %g)\n",
	    normal[0], normal[1], normal[2]);
    exit(0);
  }
    
  return(flip_normal);
}
