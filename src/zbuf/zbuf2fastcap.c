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
