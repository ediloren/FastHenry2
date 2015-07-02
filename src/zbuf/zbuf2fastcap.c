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
 
*/

#include "mulGlobal.h"
#include "zbufGlobal.h"

/*#if CAPVEW == ON*/
/*
  main interface between old zbuf code and fastcap
  - replaces functionality of zbuf.c (main())
  - dumps a geometry in .ps file format
  - panels are shaded using the vector entries of q; 
    q = NULL => no shading; use_density = TRUE => divide q's by areas
  - file name used is <ps_file_base><iter>.ps; ps_file_base is either
    the list file base, the input file base or "stdin" (see get_ps_file_info())
*/
void dump_ps_geometry(chglist, q, cond, use_ttl_chg)
charge *chglist;
double *q;
int cond;
{
  int i, j, k, numlines, numfaces, use_density;
  face **faces, **sfaces, **fastcap2faces(), **depthSortFaces();
  double normal[3], rhs, temp, dot(), *getAvg();
  double *avg, pnt[3], radius, getSphere(), getNormal();
  charge *cur_chg;
  line **lines, **getLines();
  FILE *fp;
  char str[BUFSIZ];

  extern char *ps_file_base;
  extern double view[], moffset[], distance, rotation, scale;
  extern double azimuth, elevation;
  extern char **argvals;
  extern int argcnt, g_;
  extern char *line_file;

#if 1==0  
  /* set up use density flag---not too clean; saves changes in called funcs */
  if(use_ttl_chg) use_density = FALSE;
  else use_density = TRUE;
#endif

  use_density = FALSE;

  /* convert fastcap structs to zbuf structs - COULD ELIMINATE THIS */
  faces = fastcap2faces(&numfaces, chglist, q, use_density);
  
  /* get .fig format lines in file specified with -b option */
  lines = getLines(line_file, &numlines);

  /* figure the cntr of extremal (average) coordinates of all points */
  avg = getAvg(faces, numfaces, lines, numlines, OFF);

  /* get the radius of the smallest sphere enclosing all the lines */
  radius = getSphere(avg, faces, numfaces, lines, numlines);

  /* get normal to image plane, adjust view point to be (1+distance)*radius 
     away from object center point avg, view plane (1+distance/2)*radius away
     - view coordinates taken rel to obj center but changed to absolute */
  view[0] = azimuth; view[1] = elevation;
  rhs = getNormal(normal, radius, avg, view, distance);

#if 1 == 0
  fprintf(stderr, " %d faces read\n", numfaces);
  fprintf(stderr, " %d lines read\n", numlines);
  fprintf(stderr, " average obj point: (%g %g %g), radius = %g\n", 
	  avg[0], avg[1], avg[2], radius);
  fprintf(stderr, " absolute view point: (%g %g %g)\n",
	  view[0],view[1],view[2]);
#endif

  /* set up all the normals and rhs for the faces (needed for sort) */
  initFaces(faces, numfaces, view);

  /* set up ps file name */
#if 1 == 0
  /* no longer for FastHenry */
  if(q == NULL) sprintf(str, "%s.ps", ps_file_base);
  else sprintf(str, "%s%d.ps", ps_file_base, cond);
#endif
  sprintf(str, "%s.ps", ps_file_base);

  
  /* set up the adjacency graph for the depth sort */
  fprintf(stdout, "\nSorting %d faces for %s...", numfaces, str);
  fflush(stdout);
  getAdjGraph(faces, numfaces, view, rhs, normal);
  fprintf(stdout, "done.\n");

  /* depth sort the faces */
  /*fprintf(stderr, "Starting depth sort...");*/
  sfaces = depthSortFaces(faces, numfaces);
  /*fprintf(stderr, "done.\n");*/

  /* get the 2d figure and dump to ps file */
  image(sfaces, numfaces, lines, numlines, normal, rhs, view);
  flatten(sfaces, numfaces, lines, numlines, rhs, rotation, normal, view);
  makePos(sfaces, numfaces, lines, numlines);
  scale2d(sfaces, numfaces, lines, numlines, scale, moffset);
  if(g_ == TRUE) {
    dumpCycles(sfaces, numfaces, stdout); /* DANGER - doesnt work (?) */
    dumpFaceText(sfaces, numfaces, stdout);
  }
  else {
    if((fp = fopen(str, "w")) == NULL) {
      fprintf(stderr, "dump_ps_geometry: can't open\n `%s'\nto write\n", str);
      exit(0);
    }
    fprintf(stdout, "Writing %s...", str);
    dumpPs(sfaces, numfaces, lines, numlines, fp, argvals, argcnt,use_density);
    fprintf(stdout, "done.\n");
    fclose(fp);
  }
}
/*#endif*/
