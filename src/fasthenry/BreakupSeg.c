/* this breaks a segment that is too long, into many shorter segments */
#include "induct.h"

#define DIVFACT 1  /* a factor that is always 1 except for debugging */

DivideSegs(length, charges, indsys, is_initial)
double length;
charge *charges;
SYS *indsys;
int is_initial;    /* is this the initial refinement call */
{
  SEGMENT *seg;
  int warned, warn_refine;

  warned = warn_refine = 0;

  for(seg = indsys->segment; seg != NULL; seg = seg->next) {
    if (seg->length > length/DIVFACT) {
      if (indsys->opts->auto_refine == ON || is_initial == TRUE) {
	if (seg->type == NORMAL)
	  BreakupSeg(seg, length/DIVFACT, charges, indsys);
	else if (!warned) {
	  fprintf(stdout, "DivideSegs: Warning: tried to divide an indivisable segment.\n");
	  fprintf(stdout, "  Segment length: %lf,  maximum allowed length: %lf",
		  seg->length, length/DIVFACT);
	  fprintf(stdout, "  The segment is probably part of a ground plane.\n");
	  fprintf(stdout, "  If so, decrease the partitioning level by 1 or refine the ground plane\n");
	  warned = 1;
	}
      }
      else if (!warn_refine) {
	warn_refine = 1;
	if (indsys->opts->mat_vect_prod == MULTIPOLE)
	  fprintf(stdout, "Warning: couldn't refine segments as needed because auto_refine == OFF.\n This will decrease multipole accuracy.\n");
      }

    }
  }
}

BreakupSeg(seg, maxlength, charges, indsys)
SEGMENT *seg;
double maxlength;
charge *charges;
SYS *indsys;
{
  double oldlength, x, y, z, dx, dy, dz;
  int i, j, pieces;
  NODES *node0, *node1, *node;
  SEGMENT *newseg, *lastseg, *condseg, *origsegnext;
  charge *chg;
  NODES *nodelast, *newnode;
  NODES *makenode();
  static char name[80 + MAXDEP*3];
  charge *chgend;
  charge *oldnext;
  charge *assignFil();
  NPATH *apath;
  SPATH *condpath, *lastpath, *headpath, *condbeg, *condend;
  int backwards;

  oldlength = seg->length;
  pieces = seg->length/maxlength + 1;

  // Enrico, commented line because length
  // is re-calculated later, see bug fix
  //seg->length = seg->length/pieces;
  origsegnext = seg->next;

  node0 = seg->node[0];
  node1 = seg->node[1];

  dx = (node1->x - node0->x)/ pieces;
  dy = (node1->y - node0->y)/ pieces;
  dz = (node1->z - node0->z)/ pieces;

  nodelast = node0;
  sprintf(name, "%s_0",node0->name);
  x = nodelast->x + dx;
  y = nodelast->y + dy;
  z = nodelast->z + dz;

  // Enrico, bug fix
  // The following piece of code has a potential numerical problem.
  // The problem arises when 'node0', 'node1' end points coordinates
  // are big numbers, but the segment length is small.
  // In this case, breaking up the segment, the new length of the first
  // sub-segment was calculated
  // with 'seg->length = seg->length/pieces;' (see above), while the new
  // end points coodinates by calculating 'dx', 'dy', 'dz' (small) and
  // adding them to 'node0->x', 'node0->y', 'node0->z' (big).
  // A possible fix would be: x = node0->x + dx =
  // = node0->x + (node1->x - node0->x)/ pieces = (node0->x*(pieces-1) + node1->x) / pieces;
  // but since all new sub-segments nodes positions are calculated starting
  // from 'dx', 'dy', 'dz' it is more straightforward to re-calculate
  // the length, as it is done for all sub-segments but the first one
  // later on by 'makeseg()'.
  // Note that this bug caused the problems in 'joelself.c', see the bug
  // fix in 'exact_mutual()', since there the segment length was again
  // calculated (indirectly, as a projection) but starting from the nodes
  // positions and then compared with the length stored in 'seg->length'
  seg->length = sqrt( (node0->x - x)*(node0->x - x)
		            + (node0->y - y)*(node0->y - y)
		            + (node0->z - z)*(node0->z - z) );
  //
  // end of bug fix
  //

  if (nodelast->type != NORMAL) {
    printf("Internal bug.  nodelast->type != NORMAL\n");
    exit(1);
  }

  newnode = makenode(name, indsys->num_nodes++, x, y, z, nodelast->type, NULL);
  newnode->next = node0->next;
  node0->next = newnode;

  remove_from_connected_segs(seg->node[1], seg, NULL);
  seg->node[1] = newnode;
  add_to_connected_segs(newnode, seg, NULL);

  alterFils(seg, newnode, dx, dy, dz);   /* modify the segment's fils */

  lastseg = seg;
  nodelast = newnode;
  chgend = seg->filaments[seg->num_fils-1].pchg;
  oldnext = chgend->next;
  for(i = 1; i < pieces; i++) {

    if (i != pieces - 1) {
      x = nodelast->x + dx;
      y = nodelast->y + dy;
      z = nodelast->z + dz;
      sprintf(name, "%s_%d",node0->name, i);
      newnode = makenode(name, indsys->num_nodes++, x, y, z, nodelast->type,
			 NULL);
      newnode->next = nodelast->next;
      nodelast->next = newnode;

      if (nodelast->type != NORMAL) {
	printf("Internal bug.  nodelast->type != NORMAL\n");
	exit(1);
      }

    }
    else
      newnode = node1;

    sprintf(name, "%s_%d",seg->name,i);
    newseg = makeseg(name, nodelast, newnode, seg->height, seg->width,
		     seg->sigma, seg->hinc, seg->winc, seg->r_height,
		     seg->r_width, seg->widthdir,
		     indsys->num_segs++, NORMAL, NULL);

    newseg->next = lastseg->next;
    lastseg->next = newseg;

    chgend = assignFil(newseg, &(indsys->num_fils), chgend);

    lastseg = newseg;
    nodelast = newnode;

  } /* for(i = pieces..) */
  chgend->next = oldnext;

}

alterFils(seg, node, dx, dy, dz)
SEGMENT *seg;
NODES *node;
double dx, dy, dz;
{
  FILAMENT *fil;
  charge *chg;
  int i;

  for(i = 0; i < seg->num_fils; i++) {
    fil = &(seg->filaments[i]);
    fil->x[1] = fil->x[0] + dx;
    fil->y[1] = fil->y[0] + dy;
    fil->z[1] = fil->z[0] + dz;
    fil->length = seg->length;
    fil->lenvect[XX] = fil->x[1] - fil->x[0];
    fil->lenvect[YY] = fil->y[1] - fil->y[0];
    fil->lenvect[ZZ] = fil->z[1] - fil->z[0];
    chg = fil->pchg;
    chg->max_diag = chg->min_diag = fil->length;
    chg->x = (fil->x[0] + fil->x[1])/2.0;
    chg->y = (fil->y[0] + fil->y[1])/2.0;
    chg->z = (fil->z[0] + fil->z[1])/2.0;
  }
}







