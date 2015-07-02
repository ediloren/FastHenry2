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
 
*//* this breaks a segment that is too long, into many shorter segments */
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
  
  seg->length = seg->length/pieces;
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







