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
// Enrico
#include <string.h>

double black, white;		/* densities corresponding to shades */

/*
  loads axes' lines
*/
void setupLine(axi, index, x1, y1, z1, x2, y2, z2)
double x1, y1, z1, x2, y2, z2;
double ***axi;
int index;
{
  axi[index][0][0] = x1;
  axi[index][0][1] = y1;
  axi[index][0][2] = z1;
  axi[index][1][0] = x2;
  axi[index][1][1] = y2;
  axi[index][1][2] = z2;
}

#if 1 == 0
/*
  replaces Song's make_charges_all_patches when only patches are
  to be dispalyed
  - requires minor mods to main in patran.c
*/
charge *make_charges_from_patches()
{
  int i, j, k;
  PATCH *patch_pntr;
  extern PATCH *start_patch;
  charge *start_chg, *charge_pntr;
  GRID *current_grid;
  extern GRID *start_grid;

  start_chg = NULL;

  /* for each patch, */
  patch_pntr = start_patch;
  while(patch_pntr != NULL) {

    /* create a charge struct */
    if(start_chg == NULL) {
      CALLOC(charge_pntr, 1, charge *, ON, AMSC);
      start_chg = charge_pntr;
    }
    else {
      CALLOC(charge_pntr->next, 1, charge *, ON, AMSC);
      charge_pntr = charge_pntr->next;
    }

    /* for each patch grid, find and load the coordinates
       (assumes 4 sides, just like Song does throughout) */
    for(i = 0; i < 4; i++) {
      current_grid = start_grid;
      while(current_grid != NULL) {
	if(patch_pntr->corner[i] == current_grid->ID) {
	  if(i == 0) {
	    CALLOC(charge_pntr->corner0, 3, double *, ON, AMSC)
	    for(k = 0; k < 3; k++)
		charge_pntr->corner0[k] = current_grid->coord[k];
	  }
	  else if(i == 1) {
	    CALLOC(charge_pntr->corner1, 3, double *, ON, AMSC)
	    for(k = 0; k < 3; k++)
		charge_pntr->corner1[k] = current_grid->coord[k];
	  }
	  else if(i == 2) {
	    CALLOC(charge_pntr->corner2, 3, double *, ON, AMSC)
	    for(k = 0; k < 3; k++)
		charge_pntr->corner2[k] = current_grid->coord[k];
	  }
	  else if(i == 3) {
	    CALLOC(charge_pntr->corner3, 3, double *, ON, AMSC)
	    for(k = 0; k < 3; k++)
		charge_pntr->corner3[k] = current_grid->coord[k];
	  }
	  else {
	    fprintf(stderr, "make_charges_from_patches: panel has more than 4 sides\n");
	    exit(0);
	  }
	  break;
	}
	current_grid = current_grid->next;
      }
    }

    /* put in the shape - always rectangular */
    charge_pntr->shape = 4;
    patch_pntr = patch_pntr->next;
  }
  return(start_chg);
}
#endif

/*
  set the grey levels taking into account the rd_, rc_, and rb_ options
  - levels are scaled to match the range of densities present
  - sets the values of the extremal densities (`black' and `white')
*/
void figure_grey_levels(face_list, chgs, chglist, use_density)
double *chgs;
charge *chglist;
face **face_list;
int use_density;
{
  int first, i;
  double dif;
  charge *panel;
  extern int rc_, rb_, rd_;
  extern double black, white;
  extern ITER *kq_num_list;

  /* find the minimum and maximum charge density values */
  first = TRUE;
  for(panel = chglist, i = 0; panel != NULL; panel = panel->next, i++) {
    panel->index = i;  /* for fasthenry, since panel->index never gets assgnd*/
    if(panel->dummy) continue;

    /* skip if panel is on conductor in q picture kill list */
    if(panel->surf->type == CONDTR || panel->surf->type == BOTH) {
      if(want_this_iter(kq_num_list, panel->cond)) continue;
    }

    /* skip if removing DIELEC's and panel is on dielectric i/f
       - type both interfaces are considered conductors for removal purposes */
    if(panel->surf->type == DIELEC && rd_ == TRUE) continue;

    if(first == TRUE) {
      first = FALSE;
      if(use_density) black = white = chgs[panel->index]/panel->area;
      else black = white = chgs[panel->index];
    }
    else {
      if(chgs[panel->index] == 0.0) {
	dif = 0.0;
      }
      if(use_density) {
	black = MAX(black, chgs[panel->index]/panel->area);
	white = MIN(white, chgs[panel->index]/panel->area);
      }
      else {
	black = MAX(black, chgs[panel->index]);
	white = MIN(white, chgs[panel->index]);
      }
    }

  }

  /* assign the grey levels - 0.0 = white, 1.0 = black */
  dif = black - white;
  for(panel = chglist, i = 0; panel != NULL; panel = panel->next) {
    if(panel->dummy) continue;

    /* skip if panel is on conductor in q picture kill list */
    if(panel->surf->type == CONDTR || panel->surf->type == BOTH) {
      if(want_this_iter(kq_num_list, panel->cond)) continue;
    }

    /* skip if removing DIELEC's and panel is on dielectric i/f
       - type both interfaces are considered conductors for removal purposes */
    if(panel->surf->type == DIELEC && rd_ == TRUE) continue;

    if(use_density)
	face_list[i]->greylev = (chgs[panel->index]/panel->area - white)/dif;
    else face_list[i]->greylev = (chgs[panel->index] - white)/dif;
    i++;
  }
}

/*
  attempt to read charge density information to set panel grey levels
  file format:
     (....lines of garbage....)
     Panel charges, iteration <iter>
     <chg->index> q/A = <chargedensity> (...garbage...)
     <chg->index> q/A = <chargedensity> (...garbage...)
     (...)
     <chg->index> q/A = <chargedensity> (...garbage...)
     End panel charges
     (....lines of garbage....)
  CHARGE DENSITIES ARE NOW READ DIRECTLY FROM VECTORS IN SYSTEM STRUCTS
  THIS FUNCTION IS NOT USED
*/
void get_charge_densities(q, file, iter)
double *q;
char *file;
int iter;
{
  int index, linecnt, header_found, type;
  char str1[BUFSIZ], str2[BUFSIZ], str3[BUFSIZ], linein[BUFSIZ];
  double density;
  FILE *fp, *fopen();

  if((fp = fopen(file, "r")) == NULL) {
    fprintf(stderr,
	    "get_charge_densities: can't open charge file\n  `%s'\nto read\n",
	    file);
    exit(0);
  }

  linecnt = 0;
  header_found = FALSE;
  while(fgets(linein, sizeof(linein), fp) != NULL) {
    linecnt++;

    if(!header_found) {
      /* look for "Panel charges, iteration <iter>" */
      if(sscanf(linein, "%s %s %s %d", str1, str2, str3, &index) == 4) {
	if(!strcmp(str1, "Panel") && !strcmp(str2, "charges,")
	      && !strcmp(str3, "iteration") && index == iter) {
	  header_found = TRUE;
	}
      }
      continue;
    }

    /* check for q/A line */
    if(sscanf(linein, "%d %s %s %lf", &index, str1, str2, &density) == 4) {
      if(!strcmp(str1, "q/A") && !strcmp(str2, "=")) {
	q[index] = density;
      }
    }
    /* check for end of list line */
    else if(sscanf(linein, "%s %s %s", str1, str2, str3) == 3) {
      if(!strcmp(str1, "End") && !strcmp(str2, "panel")
	 && !strcmp(str3, "charges")) {
	break;
      }
    }
    else {
      fprintf(stderr,
	      "get_charge_densities: bad charge file format, line %d:\n%s\n",
	      linecnt, linein);
      exit(0);
    }
  }

  if(!header_found) {
    fprintf(stderr,
	    "get_charge_densities: can't find iteration %d data in\n `%s'\n",
	    iter, file);
    exit(0);
  }

}

/*
  figures the corner coordinates in absolute coordinates
*/
void getAbsCoord(vec, panel, num)
double *vec;
charge *panel;
int num;
{
  double *cor = panel->corner[num];
  double *x = panel->X, *y = panel->Y, *z = panel->Z;

  vec[0] = panel->x + cor[0]*x[0] + cor[1]*y[0] + cor[2]*z[0];
  vec[1] = panel->y + cor[0]*x[1] + cor[1]*y[1] + cor[2]*z[1];
  vec[2] = panel->z + cor[0]*x[2] + cor[1]*y[2] + cor[2]*z[2];

}


/*
  transfer fastcap panel info to face structs
*/
face **fastcap2faces(numfaces, chglist, q, use_density)
int *numfaces, use_density;	/* use_density = TRUE => use q/A not q */
charge *chglist;
double *q;
{
  int i, j, dummy;
  int autmom, autlev, numMom, numLev;
  char infile[BUFSIZ];
  double dot(), getPlane(), relperm;
  charge *chgp;
  face *head, *tail, **faces;
  extern int x_, k_, q_iter, q_, rc_, rd_, rb_;
  extern double ***axes;
  extern double axeslen;
  extern char *q_file;
  extern double black, white, linewd;
  surface *surf_list, *input_surfaces();
  ssystem *sys;
  int qindex=1, cindex=1;
  double tavg[3], lavg[3];
  extern ITER *kq_num_list;

  /* transfer info to face structs (a waste but saves wrtting new fnt end) */
  for(chgp = chglist, head = NULL, *numfaces = 0; chgp != NULL;
      chgp = chgp->next) {
    if(chgp->dummy) continue;

#if RMWEDGE == ON
    /* remove all dielectric panels in the first quadrant */
    tavg[0] = tavg[1] = tavg[2] = 0.0;
    for(i = 0; i < chgp->shape; i++) {
      getAbsCoord(lavg, chgp, i);
      for(j = 0; j < 3; j++) tavg[j] += lavg[j];
    }
    if(chgp->surf->type == DIELEC && tavg[0] > 1.3 && tavg[1] > 1.3
       && tavg[2] > 1.3)
	continue;
#endif

    /* skip if panel is on conductor in q picture kill list */
    if(chgp->surf->type == CONDTR || chgp->surf->type == BOTH) {
      if(want_this_iter(kq_num_list, chgp->cond)) continue;
    }

    /* skip if removing DIELEC's and panel is on dielectric i/f
       - type both interfaces are considered conductors for removal purposes */
    if(chgp->surf->type == DIELEC && rd_ == TRUE) continue;

    /* create and link in a new face */
    if(head == NULL) {
      CALLOC(head, 1, face, ON, AMSC);
      tail = head;
    }
    else {
      CALLOC(tail->next, 1, face, ON, AMSC);
      tail = tail->next;
    }
    tail->numsides = chgp->shape;
    /* allocate for corner coordinates */
    CALLOC(tail->c, tail->numsides, double *, ON, AMSC);
    for(i = 0; i < tail->numsides; i++)
	CALLOC(tail->c[i], 3, double, ON, AMSC);
    /* xfer corner coordinates */
    for(i = 0; i < tail->numsides; i++) {
      getAbsCoord(tail->c[i], chgp, i);
    }
    /* figure and store rhs and normal */
    tail->rhs = getPlane(tail->normal, tail->c[0], tail->c[1], tail->c[2]);
    /* load grey level and line width */
    tail->greylev = GREYLEV;
    tail->width = linewd;
    tail->index = *numfaces;
    (*numfaces)++;

  }

  /* extract an array of pointers to the faces */
  CALLOC(faces, *numfaces, face *, ON, AMSC);
  for(tail = head, i = 0; tail != NULL; tail = tail->next, i++)
      faces[i] = tail;

  if(q != NULL) {
    figure_grey_levels(faces, q, chglist, use_density);
  }

  /* set up axes lines (always done - needed for alignment) */
  setupLine(axes, 0, 0.0, 0.0, 0.0, axeslen, 0.0, 0.0); /* x axis */
  setupLine(axes, 1, 0.0, 0.0, 0.0, 0.0, axeslen, 0.0); /* y axis */
  setupLine(axes, 2, 0.0, 0.0, 0.0, 0.0, 0.0, axeslen); /* z axis */
  setupLine(axes, 3, 0.85*axeslen, -0.15*axeslen, 0.0,
	    1.15*axeslen, 0.15*axeslen, 0.0); /* x marker */
  setupLine(axes, 4, 1.15*axeslen, -0.15*axeslen, 0.0,
	    0.85*axeslen, 0.15*axeslen, 0.0); /* x marker */
  setupLine(axes, 5, 0.0, axeslen, 0.0,
	    -0.15*axeslen, 1.15*axeslen, 0.0); /* y marker */
  setupLine(axes, 6, 0.0, axeslen, 0.0,
	    0.15*axeslen, 1.15*axeslen, 0.0); /* y marker */

  return(faces);
}

/*
  allows recursive reads in a .fig file with line format
  read <filename>
  only lines (not faces or fills) may be read this way
  (ie encountering "FACES" or "fills" inside a recursive file read is an error,
   unless the ignore faces/fills option (f_) has been specified)
  **** above is only vaguely correct ****
  .fig file format for specifying lines, dots and arrows:
  # <comment>
  <fromx> <fromy> <fromz>           (line from point)
  <tox> <toy> <toz> [<line width>] [<arrow size>] [<dot size>]
  ...
  r <file name>                     (hierarchical read)
  e                                 (or \n)
  dot is always on from point; arrow is always on to point
  optional arguments are only partially so: if one arg is given,
    then all those to the left must also be given
  good sizes: width 0, 1 or 2 (goes straight to postscript)
              arrow size 3.0
	      dot size 2.0
  to get just a dot, put in tiny length line and arrow size
*/
void readLines(fp, head, tail, numlines)
line **head, **tail;
int *numlines;
FILE *fp;
{
  int flines = 0, falin = 0, fflag = 1, faflag = 1, getlines = 1, i, j;
  int f_;		/* f_ == 1 => ignore face and fill info */
  char linein[BUFSIZ], **chkp, *chk, *strtok(), *cp;
  char readfile[BUFSIZ], tempc[BUFSIZ];
  double arrowsize, dotsize;
  int temp, linewd;
  FILE *fpin, *fopen();

  f_ = 1;			/* hardwire to take fill/face info as
				   equivalent to end of file */

  chkp = &chk;                  /* pointers for error checking */
  /* input lines and add to linked list */
  while(fgets(linein, sizeof(linein), fp) != NULL) {
    if(linein[0] == 'e' || linein[0] == '\0') return;
    if(linein[0] == 'r') {	/* do a recursive read */
      if(sscanf(linein, "%s %s", tempc, readfile) != 2) {
	fprintf(stderr,
		"readLines: bad recursive read line format:\n%s\n", linein);
	exit(0);
      }
      if((fpin = fopen(readfile, "r")) == NULL) {
	fprintf(stderr,
		"readLines: can't open recursive read file\n `%s'\nto read\n",
		readfile);
	exit(0);
      }
      readLines(fpin, head, tail, numlines);
      fclose(fpin);
      continue;
    }
    if(linein[0] == 'F') {
      if(f_ == 0) {
	fprintf(stderr,
		"readLines: attempt to input faces with a recursive read\n");
	exit(0);
      }
      else {
	return;
      }
    }
    if(linein[0] == '#') continue;
    if(linein[0] == 'f') {
      if(f_ == 0) {
	fprintf(stderr,
		"readLines: attempt to input fills with a recursive read\n");
	exit(0);
      }
      else {
	return;
      }
    }

    /* input lines of line information */
    if(fflag == 1) {		/* allocate a line struct; input a from line */
      if(*numlines == 0) {
	CALLOC((*tail), 1, line, ON, AMSC);
	(*head) = (*tail);
	(*tail)->prev = NULL;
      }
      else {
	CALLOC((*tail)->next, 1, line, ON, AMSC); /* link forward */
	((*tail)->next)->prev = (*tail); /* link back */
	(*tail) = (*tail)->next;
      }
      if(sscanf(linein,"%lf %lf %lf",&((*tail)->from[0]), &((*tail)->from[1]),
		&((*tail)->from[2])) != 3) {
	fprintf(stderr,"readLines: from line %d bad, '%s'\n",flines+1,linein);
	exit(0);
      }
      (*tail)->index = *numlines;
      fflag = 0;
      flines++;
    }
    else if(fflag == 0) {		/* input a to line */
      /* if arrow heads are used, line width must be specified */
      if(sscanf(linein, "%lf %lf %lf %d %lf %lf",
		&((*tail)->to[0]), &((*tail)->to[1]),
		&((*tail)->to[2]), &linewd, &arrowsize, &dotsize) != 6) {
	if(sscanf(linein, "%lf %lf %lf %d %lf",&((*tail)->to[0]),
		  &((*tail)->to[1]), &((*tail)->to[2]),
		  &linewd, &arrowsize) != 5) {
	  if(sscanf(linein, "%lf %lf %lf %d", &((*tail)->to[0]),
		    &((*tail)->to[1]), &((*tail)->to[2]), &linewd) != 4) {
	    if(sscanf(linein, "%lf %lf %lf", &((*tail)->to[0]),
		      &((*tail)->to[1]), &((*tail)->to[2])) != 3) {
	      fprintf(stderr,
		      "readLines: to line %d bad, '%s'\n",flines+1, linein);
	      exit(0);
	    }
	    linewd = LINE;
	  }
	  arrowsize = 0.0;
	}
	dotsize = 0.0;
      }
      (*tail)->width = linewd;
      (*tail)->arrow = arrowsize;
      (*tail)->dot = dotsize;
      fflag = 1;
      flines++;
      (*numlines)++;
    }
  }
  if(fflag == 0) {
    fprintf(stderr, "readLines: file ended with unmatched from line\n");
    exit(0);
  }
}

/*
  opens a .fig file and reads only lines from it - closes if faces/fills found
*/
line **getLines(line_file, numlines)
int *numlines;
char *line_file;
{
  int i;
  FILE *fp;
  line *head, *tail, **linesout;

  *numlines = 0;

  if(line_file == NULL) return(NULL);

  if((fp = fopen(line_file, "r")) == NULL) {
    fprintf(stderr, "getLines: can't open .fig file\n `%s'\nto read\n",
	    line_file);
    exit(0);
  }

  readLines(fp, &head, &tail, numlines);
  fclose(fp);

  /* extract array of pointers to line structs */
  CALLOC(linesout, (*numlines), line *, ON, AMSC);
  for(i = 0; i < *numlines; i++) {
    linesout[i] = head;
    head = head->next;
  }
  return(linesout);
}

/*
  figure the bounding box and write ps file line
*/
void getBndingBox(faces, numfaces, lines, numlines, lowx, lowy, fp, axes)
face **faces;
line **lines;
double ***axes;
int numfaces, numlines, *lowx, *lowy;
FILE *fp;
{
  int upx, upy;
  double xmax, ymax, minx, miny;
  int i, j;
  extern int x_;

  /* find the smallest and largest x and y coordinates (assumed pos) */
  xmax = ymax = 0.0;
  for(i = 0; i < 7 && x_ == TRUE; i++) { /* check axes */
    for(j = 0; j < 2; j++) {
      if(i == 0 && j == 0) {
	minx = axes[i][j][0];
	miny = axes[i][j][1];
      }
      else {
	minx = MIN(minx, axes[i][j][0]);
	miny = MIN(miny, axes[i][j][1]);
      }
      xmax = MAX(xmax, axes[i][j][0]);
      ymax = MAX(ymax, axes[i][j][1]);
    }
  }
  for(i = 0; i < numfaces; i++) { /* check faces */
    for(j = 0; j < faces[i]->numsides; j++) {
      if(i == 0 && j == 0 && x_ == FALSE) {
	minx = faces[i]->c[j][0];
	miny = faces[i]->c[j][1];
      }
      else {
	minx = MIN(minx, faces[i]->c[j][0]);
	miny = MIN(miny, faces[i]->c[j][1]);
      }
      xmax = MAX(xmax, faces[i]->c[j][0]);
      ymax = MAX(ymax, faces[i]->c[j][1]);
    }
  }
  for(i = 0; i < numlines; i++) { /* check lines */
    if(i == 0 && x_ == FALSE && numfaces == 0) {
      minx = MIN(lines[i]->from[0], lines[i]->to[0]);
      miny = MIN(lines[i]->from[1], lines[i]->to[1]);
    }
    else {
      minx = MIN(minx, lines[i]->from[0]);
      miny = MIN(miny, lines[i]->from[1]);
      minx = MIN(minx, lines[i]->to[0]);
      miny = MIN(miny, lines[i]->to[1]);
    }
    xmax = MAX(xmax, lines[i]->to[0]);
    xmax = MAX(xmax, lines[i]->from[0]);
    ymax = MAX(ymax, lines[i]->to[1]);
    ymax = MAX(ymax, lines[i]->from[1]);
  }

  *lowx = minx-2;		/* note 2pnt offset and truncation */
  *lowy = miny-2;
  upx = xmax+2;
  upy = ymax+2;
  fprintf(fp, "%%%%BoundingBox: %d %d %d %d\n", *lowx, *lowy, upx, upy);
}

/*
  dump axes to ps file
*/
void dumpAxes(axi, fp)
double ***axi;
FILE *fp;
{
  int i;

  for(i = 0; i < 7; i++) {	/* loop on axes' lines (pointers too) */
    fprintf(fp, "%g %g moveto\n", axi[i][0][0], axi[i][0][1]);
    fprintf(fp, "%g %g lineto\n", axi[i][1][0], axi[i][1][1]);
    fprintf(fp, "%g setlinewidth %d setlinecap %d setlinejoin ",
	      AXEWID, LINCAP, LINJIN);
    fprintf(fp, " 0 setgray  stroke\n");
  }
  /*fprintf(stderr, "Axes inserted\n");*/
}


/*
  copy the body of the header to the output file
*/
void copyBody(fp)
FILE *fp;
{
  static char str[] = "\
%%%%DocumentProcSets: FreeHand_header 2 0 \n\
%%%%DocumentSuppliedProcSets: FreeHand_header 2 0 \n\
%%%%ColorUsage: Color \n\
%%%%CMYKProcessColor: 0 0 0 0.1  (10%% gray) \n\
%%%%+ 0 0 0 0.2  (20%% gray) \n\
%%%%+ 0 0 0 0.4  (40%% gray) \n\
%%%%+ 0 0 0 0.6  (60%% gray) \n\
%%%%+ 0 0 0 0.8  (80%% gray) \n\
%%%%EndComments \n\
%%%%BeginProcSet: FreeHand_header 2 0 \n\
/FreeHandDict 200 dict def \n\
FreeHandDict begin \n\
/currentpacking where{pop true setpacking}if \n\
/bdf{bind def}bind def \n\
/bdef{bind def}bdf \n\
/xdf{exch def}bdf \n\
/ndf{1 index where{pop pop pop}{dup xcheck{bind}if def}ifelse}bdf \n\
/min{2 copy gt{exch}if pop}bdf \n\
/max{2 copy lt{exch}if pop}bdf \n\
/dr{transform .25 sub round .25 add \n\
exch .25 sub round .25 add exch itransform}bdf \n\
/curveto{dr curveto}bdf \n\
/lineto{dr lineto}bdf \n\
/moveto{dr moveto}bdf \n\
/graystep 1 256 div def \n\
/bottom -0 def \n\
/delta -0 def \n\
/frac -0 def \n\
/left -0 def \n\
/numsteps -0 def \n\
/numsteps1 -0 def \n\
/radius -0 def \n\
/right -0 def \n\
/top -0 def \n\
/x -0 def \n\
/y -0 def \n\
/df currentflat def \n\
/tempstr 1 string def \n\
/clipflatness 3 def \n\
/inverted? \n\
0 currenttransfer exec .5 ge def \n\
/concatprocs{ \n\
/proc2 exch cvlit def/proc1 exch cvlit def \n\
/newproc proc1 length proc2 length add array def \n\
newproc 0 proc1 putinterval newproc proc1 length proc2 putinterval \n\
newproc cvx}bdf \n\
/storerect{/top xdf/right xdf/bottom xdf/left xdf}bdf \n\
/rectpath{newpath left bottom moveto left top lineto \n\
right top lineto right bottom lineto closepath}bdf \n\
/sf{dup 0 eq{pop df dup 3 mul}{dup} ifelse /clipflatness xdf setflat}bdf \n";

static char str2[] = "\
version cvr 38.0 le \n\
{/setrgbcolor{ \n\
currenttransfer exec 3 1 roll \n\
currenttransfer exec 3 1 roll \n\
currenttransfer exec 3 1 roll \n\
setrgbcolor}bdf}if \n\
/gettint{0 get}bdf \n\
/puttint{0 exch put}bdf \n\
/vms{/vmsv save def}bdf \n\
/vmr{vmsv restore}bdf \n\
/vmrs{vmr vms}bdf \n\
/CD{/NF exch def \n\
{exch dup/FID ne{exch NF 3 1 roll put} \n\
{pop pop}ifelse}forall NF}bdf \n\
/MN{1 index length/Len exch def \n\
dup length Len add string dup \n\
Len 4 -1 roll putinterval dup 0 4 -1 roll putinterval}bdf \n\
/RC{256 string cvs(|______)anchorsearch \n\
{1 index MN cvn/NewN exch def cvn \n\
findfont dup maxlength dict CD dup/FontName NewN put dup \n\
/Encoding MacVec put NewN exch definefont pop}{pop}ifelse}bdf \n\
/RF{dup FontDirectory exch known{pop}{RC}ifelse}bdf \n\
/FF{dup 256 string cvs(|______)exch MN cvn dup FontDirectory exch known \n\
{exch}if pop findfont}bdf \n\
userdict begin /BDFontDict 20 dict def end \n\
BDFontDict begin \n\
/bu{}def \n\
/bn{}def \n\
/setTxMode{pop}def \n\
/gm{moveto}def \n\
/show{pop}def \n\
/gr{pop}def \n\
/fnt{pop pop pop}def \n\
/fs{pop}def \n\
/fz{pop}def \n\
/lin{pop pop}def \n\
end \n\
/MacVec 256 array def \n\
MacVec 0 /Helvetica findfont \n\
/Encoding get 0 128 getinterval putinterval \n\
MacVec 127 /DEL put MacVec 16#27 /quotesingle put MacVec 16#60 /grave put \n\
/NUL/SOH/STX/ETX/EOT/ENQ/ACK/BEL/BS/HT/LF/VT/FF/CR/SO/SI \n\
/DLE/DC1/DC2/DC3/DC4/NAK/SYN/ETB/CAN/EM/SUB/ESC/FS/GS/RS/US \n\
MacVec 0 32 getinterval astore pop \n\
/Adieresis/Aring/Ccedilla/Eacute/Ntilde/Odieresis/Udieresis/aacute \n\
/agrave/acircumflex/adieresis/atilde/aring/ccedilla/eacute/egrave \n\
/ecircumflex/edieresis/iacute/igrave/icircumflex/idieresis/ntilde/oacute \n\
/ograve/ocircumflex/odieresis/otilde/uacute/ugrave/ucircumflex/udieresis \n\
/dagger/degree/cent/sterling/section/bullet/paragraph/germandbls \n\
/register/copyright/trademark/acute/dieresis/notequal/AE/Oslash \n\
/infinity/plusminus/lessequal/greaterequal/yen/mu/partialdiff/summation \n";

static char str3[] = "\
/product/pi/integral/ordfeminine/ordmasculine/Omega/ae/oslash \n\
/questiondown/exclamdown/logicalnot/radical/florin/approxequal/Delta/guillemotleft \n\
/guillemotright/ellipsis/nbspace/Agrave/Atilde/Otilde/OE/oe \n\
/endash/emdash/quotedblleft/quotedblright/quoteleft/quoteright/divide/lozenge \n\
/ydieresis/Ydieresis/fraction/currency/guilsinglleft/guilsinglright/fi/fl \n\
/daggerdbl/periodcentered/quotesinglbase/quotedblbase \n\
/perthousand/Acircumflex/Ecircumflex/Aacute \n\
/Edieresis/Egrave/Iacute/Icircumflex/Idieresis/Igrave/Oacute/Ocircumflex \n\
/apple/Ograve/Uacute/Ucircumflex/Ugrave/dotlessi/circumflex/tilde \n\
/macron/breve/dotaccent/ring/cedilla/hungarumlaut/ogonek/caron \n\
MacVec 128 128 getinterval astore pop \n\
/fps{currentflat exch dup 0 le{pop 1}if \n\
{dup setflat 3 index stopped \n\
{1.3 mul dup 3 index gt{pop setflat pop pop stop}if}{exit}ifelse \n\
}loop pop setflat pop pop \n\
}bdf \n\
/fp{100 currentflat fps}bdf \n\
/rfp{clipflatness currentflat fps}bdf \n\
/fcp{100 clipflatness fps}bdf \n\
/fclip{{clip}fcp}bdf \n\
/feoclip{{eoclip}fcp}bdf \n\
end %%. FreeHandDict \n\
%%%%EndProcSet \n\
%%%%BeginSetup \n\
FreeHandDict begin \n\
/ccmyk{dup 5 -1 roll sub 0 max exch}ndf \n\
/setcmykcolor{1 exch sub ccmyk ccmyk ccmyk pop setrgbcolor}ndf \n\
/setcmykcoloroverprint{4{dup -1 eq{pop 0}if 4 1 roll}repeat setcmykcolor}ndf \n\
/findcmykcustomcolor{5 /packedarray where{pop packedarray}{array astore readonly}ifelse}ndf \n\
/setcustomcolor{exch aload pop pop 4{4 index mul 4 1 roll}repeat setcmykcolor pop}ndf \n\
/setseparationgray{1 exch sub dup dup dup setcmykcolor}ndf \n\
/setoverprint{pop}ndf \n\
/currentoverprint false ndf \n\
/colorimage{pop pop \n\
[5 -1 roll/exec cvx 6 -1 roll/exec cvx 7 -1 roll/exec cvx 8 -1 roll/exec cvx \n\
/exch cvx/pop cvx/exch cvx/pop cvx/exch cvx/pop cvx/invbuf cvx]cvx image} \n\
%%. version 47.1 of Postscript defines colorimage incorrectly (rgb model only) \n\
version cvr 47.1 le{userdict begin bdf end}{ndf}ifelse \n\
/customcolorimage{pop image}ndf \n\
/separationimage{image}ndf \n\
/newcmykcustomcolor{6 /packedarray where{pop packedarray}{array astore readonly}ifelse}ndf \n\
/inkoverprint false ndf \n\
/setinkoverprint{pop}ndf \n\
/overprintprocess{pop}ndf \n\
/setspotcolor \n\
{spots exch get 0 5 getinterval exch setcustomcolor}ndf \n\
/currentcolortransfer{currenttransfer dup dup dup}ndf \n\
/setcolortransfer{systemdict begin settransfer end pop pop pop}ndf \n";

static char str4[] = "\
/setimagecmyk{dup length 4 eq \n\
{aload pop} \n\
{aload pop spots exch get 0 4 getinterval aload pop 4 \n\
{4 index mul 4 1 roll}repeat 5 -1 roll pop} ifelse \n\
systemdict /colorimage known{version cvr 47.1 gt}{false}ifelse \n\
not{pop 1 currentgray sub}if \n\
/ik xdf /iy xdf /im xdf /ic xdf \n\
}ndf \n\
/setcolor{dup length 4 eq \n\
{aload overprintprocess setcmykcolor} \n\
{aload 1 get spots exch get 5 get setinkoverprint setspotcolor} \n\
ifelse}ndf \n\
/bc2[0 0]def \n\
/bc4[0 0 0 0]def \n\
/c1[0 0 0 0]def \n\
/c2[0 0 0 0]def \n\
/absmax{2 copy abs exch abs gt{exch}if pop}bdf \n\
/calcstep \n\
{c1 length 4 eq \n\
{ \n\
0 1 3 \n\
{c1 1 index get \n\
c2 3 -1 roll get \n\
sub \n\
}for \n\
absmax absmax absmax \n\
} \n\
{ \n\
bc2 c1 1 get 1 exch put \n\
c1 gettint c2 gettint \n\
sub abs \n\
}ifelse \n\
graystep div abs round dup 0 eq{pop 1}if \n\
dup /numsteps xdf 1 sub dup 0 eq{pop 1}if /numsteps1 xdf \n\
}bdf \n\
/cblend{ \n\
c1 length 4 eq \n\
{ \n\
0 1 3 \n\
{bc4 exch \n\
c1 1 index get \n\
c2 2 index get \n\
1 index sub \n\
frac mul add put \n\
}for bc4 \n\
}{ \n\
bc2 \n\
c1 gettint \n\
c2 gettint \n\
1 index sub \n\
frac mul add \n\
puttint bc2 \n\
}ifelse \n\
setcolor \n\
}bdf \n";

static char str5[] = "\
/logtaper{/frac frac 9 mul 1 add log def}bdf \n\
/imbits 1 def \n\
/iminv false def \n\
/invbuf{0 1 2 index length 1 sub{dup 2 index exch get 255 exch sub 2 index 3 1 roll put}for}bdf \n\
/cyanrp{currentfile cyanbuf readhexstring pop iminv{invbuf}if}def \n\
/magentarp{cyanbuf magentabuf copy}bdf \n\
/yellowrp{cyanbuf yellowbuf copy}bdf \n\
/blackrp{cyanbuf blackbuf copy}bdf \n\
/fixtransfer{ \n\
dup{ic mul ic sub 1 add}concatprocs exch \n\
dup{im mul im sub 1 add}concatprocs exch \n\
dup{iy mul iy sub 1 add}concatprocs exch \n\
{ik mul ik sub 1 add}concatprocs \n\
currentcolortransfer \n\
5 -1 roll exch concatprocs 7 1 roll \n\
4 -1 roll exch concatprocs 6 1 roll \n\
3 -1 roll exch concatprocs 5 1 roll \n\
concatprocs 4 1 roll \n\
setcolortransfer \n\
}bdf \n";

fprintf(fp, "%s%s%s%s%s", str, str2, str3, str4, str5);

}

/*
  numbers the faces for checking
*/
void numberFaces(faces, numfaces, fp)
face **faces;
int numfaces;
FILE *fp;
{
  int i, j, mid[2];
  double cent[2];

  /* put face number at average point of each face */
  for(i = 0; i < numfaces; i++) {
    /* figure midpoint, truncate (truncation not really necessary) */
    for(j = 0, cent[0] = cent[1] = 0.0; j < faces[i]->numsides; j++) {
      cent[0] += faces[i]->c[j][0]; /* x coordinate sum */
      cent[1] += faces[i]->c[j][1]; /* y coordinate sum */
    }
    mid[0] = cent[0]/(((double) faces[i]->numsides)); /* average x */
    mid[1] = cent[1]/(((double) faces[i]->numsides)); /* average y */
    /* dump a label with associated garbage */
    fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
    fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
    fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
    fprintf(fp, "%d %d moveto\n", mid[0], mid[1]);
    fprintf(fp, "0 0 32 0 0 (F%d) ts\n}\n", i);
    fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");
  }
}

/*
  number a face for checking - used to cover up obscured faces' numbers
*/
void numberFace(fac, fp)
face *fac;
FILE *fp;
{
  int j, mid[2];
  double cent[2];

  /* figure midpoint, truncate (truncation not really necessary) */
  for(j = 0, cent[0] = cent[1] = 0.0; j < fac->numsides; j++) {
    cent[0] += fac->c[j][0]; /* x coordinate sum */
    cent[1] += fac->c[j][1]; /* y coordinate sum */
  }
  mid[0] = cent[0]/(((double) fac->numsides)); /* average x */
  mid[1] = cent[1]/(((double) fac->numsides)); /* average y */
  /* dump a label with associated garbage */
  fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
  fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
  fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
  fprintf(fp, "%d %d moveto\n", mid[0], mid[1]);
  fprintf(fp, "0 0 32 0 0 (%d) ts\n}\n", fac->index);
  fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");

}

/*
  dumps adjacency graph as a ps file - uses both input order and graph order
*/
void dumpAdjGraph(faces, numfaces, fp)
face **faces;
int numfaces;
FILE *fp;
{
  int f, i;
  double x, y;			/* current point in plot */
  double stepx, stepy;	       	/* step in x and y directions */
  double font;			/* font size */

  /* start the input numbered graph refered to lower left corner
     - row numbers on right because it's easier */
  /* set up the sizes - font never bigger than FONT; stepx, stepy <=1.25FONT */
  stepx = MIN(1.25*FONT, (IMAGEX-OFFSETX)/(double)numfaces);
  stepy = MIN(1.25*FONT, (IMAGEY-OFFSETY)/(double)numfaces);
  font = MIN(stepx, stepy)/1.25;
  x = OFFSETX + numfaces*stepx;
  y = OFFSETY + numfaces*stepy;

  /* number columns - mark those divisible by ten */
  for(f = 0; f < numfaces; f++) {
    if(f % 10 == 0 && f != 0) fprintf(fp, "%g %g dia\n", x-f*stepx, y+stepy);
  }

  /* number each row and fill it in - input ordering
  for(f = 0; f < numfaces; f++) { */
    /* dump a row label with associated garbage
    fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
    fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
    fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
    fprintf(fp, "%g %g moveto\n", x+stepx, y-faces[f]->index*stepy);
    fprintf(fp, "0 0 32 0 0 (%d) ts\n}\n", faces[f]->index);
    fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n"); */
    /* dump dot if an edge
    for(i = 0; i < faces[f]->numbehind; i++) {
      fprintf(fp, "%g %g dot\n",
	      x-(faces[f]->behind)[i]->index*stepx, y-faces[f]->index*stepy);
    }
  } */
  /* dump title
  fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
  fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
  fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
  fprintf(fp, "%g %g moveto\n", OFFSETX, y+FONT);
  fprintf(fp, "0 0 32 0 0 (Input Ordering) ts\n}\n");
  fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n"); */

  /* y += (numfaces*stepy + 3*FONT);	/* offset 2nd array */

  /* number each row and fill it in - graph ordering */
  for(f = 0; f < numfaces; f++) {
    fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
    fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
    fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
    fprintf(fp, "%g %g moveto\n", x+stepx, y-faces[f]->depth*stepy);
    fprintf(fp, "0 0 32 0 0 (%d) ts\n}\n", faces[f]->depth);
    fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");
    for(i = 0; i < faces[f]->numbehind; i++) {
      fprintf(fp, "%g %g dot\n",
	      x-(faces[f]->behind)[i]->depth*stepx, y-faces[f]->depth*stepy);
    }
  }
  fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
  fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
  fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", FONT, FONT);
  fprintf(fp, "%g %g moveto\n", OFFSETX, y+FONT);
  fprintf(fp, "0 0 32 0 0 (Graph Ordering) ts\n}\n");
  fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");
}

/*
  dump face graph as a text file
*/
void dumpFaceText(faces, numfaces, fp)
face **faces;
int numfaces;
FILE *fp;
{
  int f, i, first = 0, k;

  fprintf(fp, "depth order (input order) - lower numbers are deeper\n");
  for(f = 0; f < numfaces; f++) {
    fprintf(fp, "%d (%d):", faces[f]->depth, faces[f]->index);
    for(i = 0; i < faces[f]->numbehind && faces[f]->behind != NULL; i++) {
      fprintf(fp, " %d (%d)", (faces[f]->behind)[i]->depth,
	      (faces[f]->behind)[i]->index);
      if(i % 5 == 0 && i != 0) fprintf(fp, "\n");
    }
    if((i-1) % 5 != 0 || i == 1) fprintf(fp, "\n");
  }

  /* check to see that ordering is consistent with deeper lists */
  for(f = 0; f < numfaces; f++) {
    for(k = 0; k < faces[f]->numbehind; k++) {
      if(faces[f]->depth >= (faces[f]->behind)[k]->depth) {
	if(first == 0) {
	  first = 1;
	  fprintf(fp, "\nVertices whose depth lists are inconsistent\n");
	}
	fprintf(fp, "%d (%d):", faces[f]->depth, faces[f]->index);
	for(i = 0; i < faces[f]->numbehind && faces[f]->behind != NULL; i++) {
	  fprintf(fp, " %d (%d)", (faces[f]->behind)[i]->depth,
		  (faces[f]->behind)[i]->index);
	  if(i % 5 == 0 && i != 0) fprintf(fp, "\n");
	}
	if((i-1) % 5 != 0 || i == 1) fprintf(fp, "\n");
	break;
      }
    }
  }

}

/*
  dumps a line of chars in the Aldus FreeHand ps format
*/
void dump_line_as_ps(fp, psline, x_position, y_position, font_size)
char *psline;
double x_position, y_position, font_size;
FILE *fp;
{
  fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
  fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
  fprintf(fp, "f1 [%g 0 0 %g 0 0] makesetfont\n", font_size, font_size);
  fprintf(fp, "%g %g moveto\n", x_position, y_position);
  fprintf(fp, "0 0 32 0 0 (%s) ts\n}\n", psline);
  fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");
}

/*
  dump nblocks blocks with shades between white and black, labeled with density
*/
void dump_shading_key(fp, nblocks, precision, font_size, use_density)
int nblocks, precision, use_density;
double font_size;
FILE *fp;
{
  int i;
  double x_right, y_top, block_hgt, block_x, block_y, string_x, diddle_x;
  double grey_step, grey_lev, density, density_step, white_width;
  char linein[BUFSIZ], ctrl[BUFSIZ];
  extern double black, white, linewd;

  x_right = OFFSETX + IMAGEX;
  y_top = OFFSETY + IMAGEY;
  block_hgt = KEYHGT/(double)nblocks;
  block_x = x_right - KEYWID;
  block_y = y_top;
  /*string_x = block_x - font_size*(6.0/2.0 + (double)precision);*/
  string_x = block_x + font_size/2.0;

  /* writing raw ps so 1.0 = white, 0.0 = black */
  grey_lev = 0.0;
  grey_step = 1.0/(double)(nblocks-1);
  density = black;
  density_step = (black-white)/(double)(nblocks-1);
  /*white_width = 3.0 + (double)precision;*/
  white_width = KEYWID;

  /* write the key title */
  if(use_density) {
    strcpy(linein, "DENSITY, statC/m^2");
    diddle_x = font_size;
  }
  else {
    strcpy(linein, "CHARGE, statC");
    diddle_x = font_size/2.0;
  }
  dump_line_as_ps(fp, linein, string_x-diddle_x, y_top + font_size/2.0,
		  font_size);

  for(i = 0; i < nblocks; i++) {
    /* write a fill with border for the key block */
    /* dump the fill */
    fprintf(fp, "%g %g moveto\n", block_x, block_y);
    fprintf(fp, "%g %g lineto\n", block_x + KEYWID, block_y);
    fprintf(fp, "%g %g lineto\n", block_x + KEYWID, block_y - block_hgt);
    fprintf(fp, "%g %g lineto\n", block_x, block_y - block_hgt);
    fprintf(fp, "closepath\n");
    fprintf(fp, " %g setgray fill\n", grey_lev);

    /* dump the white out for the label */
    fprintf(fp, "%g %g moveto\n", block_x, block_y);
    fprintf(fp, "%g %g lineto\n", block_x + white_width, block_y);
    fprintf(fp, "%g %g lineto\n",
	    block_x + white_width, block_y - font_size - font_size/10.0);
    fprintf(fp, "%g %g lineto\n",
	    block_x, block_y - font_size - font_size/10.0);
    fprintf(fp, "closepath\n");
    fprintf(fp, " 1.0 setgray fill\n");

    /* dump the outline */
    fprintf(fp, "%g %g moveto\n", block_x, block_y);
    fprintf(fp, "%g %g lineto\n", block_x + KEYWID, block_y);
    fprintf(fp, "%g %g lineto\n", block_x + KEYWID, block_y - block_hgt);
    fprintf(fp, "%g %g lineto\n", block_x, block_y - block_hgt);
    fprintf(fp, "closepath\n");
    fprintf(fp, "%g setlinewidth %d setlinecap %d setlinejoin ",
	    linewd, LINCAP, LINJIN);
    fprintf(fp, " 0 setgray  stroke\n");

    /* dump the label */
    sprintf(ctrl, "%%.%dg", precision);
    sprintf(linein, ctrl, density);
    dump_line_as_ps(fp, linein, string_x, block_y - font_size, font_size);

    block_y -= block_hgt;
    density -= density_step;
    grey_lev += grey_step;
  }
}

/*
  numbers the lines for checking
*/
void numberLines(lines, numlines, fp)
line **lines;
int numlines;
FILE *fp;
{
  int i, mid[2];

  /* put line number on midpoint of each line */
  for(i = 0; i < numlines; i++) {
    /* figure midpoint, truncate (truncation not really necessary) */
    mid[0] = ((lines[i]->from)[0] + (lines[i]->to)[0])/2;
    mid[1] = ((lines[i]->from)[1] + (lines[i]->to)[1])/2;
    /* dump a label with associated garbage */
    fprintf(fp, "%%%%IncludeFont: Times-Roman\n");
    fprintf(fp, "/f1 /|______Times-Roman dup RF findfont def\n{\n");
    fprintf(fp, "f1 [%d 0 0 %d 0 0] makesetfont\n", FONT, FONT);
    fprintf(fp, "%d %d moveto\n", mid[0], mid[1]);
    fprintf(fp, "0 0 32 0 0 (%d) ts\n}\n", lines[i]->index);
    fprintf(fp, "[0 0 0 1]\nsts\nvmrs\n");
  }
}

/*
  lobotomized version of dumpPs in orthoPs.c - dumps lines/arrows
*/
void dumpLines(fp, lines, numlines)
line **lines;
int numlines;
FILE *fp;
{
  int i, j, w_;
  double temp[3], temp1[3], x, y;

  if(fp == NULL) {
    fprintf(stderr, "dumpLines: null ps file pointer\n");
    exit(0);
  }

  w_ = 0;			/* hardwire for no width override */

  /* dump the lines  */
  for(i = 0; i < numlines; i++) {
    fprintf(fp, "%%%% begin line %d\n", lines[i]->index);
    fprintf(fp, "0 sf\nnewpath\n");
    x = (lines[i]->from)[0];
    y = (lines[i]->from)[1];
    fprintf(fp, "%g %g moveto\n", x, y);
    x = (lines[i]->to)[0];
    y = (lines[i]->to)[1];
    fprintf(fp, "%g %g lineto\n", x, y);
    fprintf(fp, "gsave\n");
    if(lines[i]->width == DASHED) {
      if(w_ == 0)
	  fprintf(fp, "%d setlinewidth 1 setlinecap 0 setlinejoin 3.863693",
	      DASWTH);
      else fprintf(fp, "%d setlinewidth 1 setlinecap 0 setlinejoin 3.863693",
	      OVRWTH);
      fprintf(fp,
	      " setmiterlimit [0 0 0 1]setcolor [2 4] 0 setdash {stroke}fp\n");
    }
    else {
      if(w_ == 0)
	  fprintf(fp, "%d setlinewidth 1 setlinecap 0 setlinejoin 3.863693",
		  lines[i]->width);
      else
	  fprintf(fp, "%d setlinewidth 1 setlinecap 0 setlinejoin 3.863693",
		  OVRWTH);
      fprintf(fp, " setmiterlimit [0 0 0 1]setcolor  {stroke}fp\n");
    }
    fprintf(fp, "grestore\n");
    if(lines[i]->arrow > 0.0) {	/* put arrow head on to side if desired */
      /* figure unit vector from `to' point to `from' point */
      for(j = 0; j < 2; j++) temp[j] = lines[i]->from[j]-lines[i]->to[j];
      temp1[0] = sqrt(temp[0]*temp[0]+temp[1]*temp[1]);
      for(j = 0; j < 2; j++) temp[j] /= temp1[0];
      for(j = 0; j < 2; j++)	/* figure unit perpendicular */
	  temp1[j] =
	      1.0/(temp[j]*sqrt(1.0/(temp[0]*temp[0])+1.0/(temp[1]*temp[1])));
      temp1[0] = -temp1[0];
      /* draw the arrow */
      fprintf(fp, "%%%% Begin arrow head for line %d\n", i);
      fprintf(fp, "%g %g moveto\n", lines[i]->to[0], lines[i]->to[1]);
      fprintf(fp, "%g %g lineto\n",
	      lines[i]->to[0]+lines[i]->arrow*ALEN*temp[0]
	      +lines[i]->arrow*(AWID/2)*temp1[0],
	      lines[i]->to[1]+lines[i]->arrow*ALEN*temp[1]
	      +lines[i]->arrow*(AWID/2)*temp1[1]);
      fprintf(fp, "%g %g lineto\n",
	      lines[i]->to[0]+lines[i]->arrow*ALEN*temp[0]
	      -lines[i]->arrow*(AWID/2)*temp1[0],
	      lines[i]->to[1]+lines[i]->arrow*ALEN*temp[1]
	      -lines[i]->arrow*(AWID/2)*temp1[1]);
      fprintf(fp, "closepath\n");
      fprintf(fp, " 0 setgray fill\n");
      /* put dot on from end of line, if called for */
      if(lines[i]->dot > 0.0)
	  fprintf(fp, "%g %g %g 0 360 arc closepath fill\n",
		  lines[i]->from[0], lines[i]->from[1], lines[i]->dot*DOTSIZ);
    }
  }
}


/*
  dump faces in ps Aldus FreeHand format - assumes header body in afhpsheader
*/
void dumpPs(faces, numfaces, lines, numlines, fp, argv, argc, use_density)
face **faces;
line **lines;
char **argv;
int numfaces, numlines, argc, use_density;
FILE *fp;
{
  int i, j, f, lowx, lowy;
  extern int s_, n_, g_, c_, x_, q_, rk_, f_, m_; /* command line flags */
  extern double ***axes;
  extern double black, white;
  char linein[BUFSIZ];
  double dot(), len, temp[2], xc, yc;
  double x, y;

  /* print the lines before the bounding box */
  fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(fp, "%%%%Creator: FreeHand\n");
  fprintf(fp, "%%%%Title: test.ps\n");
  fprintf(fp, "%%%%CreationDate: 4/19/90 10:47 AM\n");

  getBndingBox(faces, numfaces, lines, numlines,
	       &lowx, &lowy, fp, axes); /* prnt bnding box */
  copyBody(fp);			/* copys the body of the header from
				   "afhpsheader" */

  /* dump the text header if needed */
  if(n_ == TRUE || g_ == TRUE || c_ == TRUE || q_ == TRUE) {
    fprintf(fp, "/textopf false def\n/curtextmtx{}def\n/otw .25 def\n");
    fprintf(fp, "/msf{dup/curtextmtx xdf makefont setfont}bdf\n");
    fprintf(fp, "/makesetfont/msf load def\n");
    fprintf(fp, "/curtextheight{.707104 .707104 curtextmtx dtransform\n");
    fprintf(fp, "dup mul exch dup mul add sqrt}bdf\n");
    fprintf(fp, "/ta{1 index\n{tempstr 0 2 index put tempstr 2 index\n");
    fprintf(fp, "gsave exec grestore\ntempstr stringwidth rmoveto\n");
    fprintf(fp, "5 index eq{6 index 6 index rmoveto}if\n");
    fprintf(fp, "3 index 3 index rmoveto\n");
    fprintf(fp, "}forall 7{pop}repeat}bdf\n");
    fprintf(fp,
	    "/sts{setcolor textopf setoverprint/ts{awidthshow}def exec}bdf\n");
    fprintf(fp, "/stol{setlinewidth setcolor textopf setoverprint newpath\n");
    fprintf(fp, "/ts{{false charpath stroke}ta}def exec}bdf\n");
  }

  /* print rest of header (starting with after /currentpacking where...) */
  fprintf(fp, "/currentpacking where{pop false setpacking}if\n");
  fprintf(fp, "%%%%EndSetup\n");
  fprintf(fp, "/spots[1 0 0 0 (Process Cyan) false newcmykcustomcolor\n");
  fprintf(fp, "0 1 0 0 (Process Magenta) false newcmykcustomcolor\n");
  fprintf(fp, "0 0 1 0 (Process Yellow) false newcmykcustomcolor\n");
  fprintf(fp,
	  "0 0 0 1 (Process Black) false newcmykcustomcolor\n]def\nvms\n");

  /* dump command line as a comment */
  fprintf(fp, "%%%% ");
  for(i = 0; i < argc; i++) fprintf(fp, " %s", argv[i]);
  fprintf(fp, "\n");

  if(x_ == TRUE) dumpAxes(axes, fp); /* dump axes if called for */

  /* for each face - dump fill, then outline - assumes depth ordering */
  for(f = 0; f < numfaces; f++) {
    if(!f_) {
      /* dump the fill */
      fprintf(fp, "%%%% Begin face %d\n", f);
      /* fprintf(fp, "0 sf\nnewpath\n"); */
      /* fprintf(fp, "newpath\n");*/
      fprintf(fp, "%g %g moveto\n", faces[f]->c[0][0], faces[f]->c[0][1]);
      for(i = 1; i < faces[f]->numsides; i++) {
	fprintf(fp, "%g %g lineto\n", faces[f]->c[i][0], faces[f]->c[i][1]);
      }
      fprintf(fp, "closepath\n");
      /* fprintf(fp, "gsave\n");*/
      /* fprintf(fp, "[0 0 0 %g]setcolor  {fill}fp\ngrestore\n", GREYLEV); */
      /* fprintf(fp, "[0 0 0 %g]setcolor  fill\n", GREYLEV);*/
      fprintf(fp, " %g setgray fill\n", 1.0-faces[f]->greylev);
      /* fprintf(fp, "grestore\n");*/
    }

    /* dump the outline */
    /* fprintf(fp, "0 sf\nnewpath\n"); */
    /* fprintf(fp, "newpath\n");*/
    fprintf(fp, "%g %g moveto\n", faces[f]->c[0][0], faces[f]->c[0][1]);
    for(i = 1; i < faces[f]->numsides; i++) {
      fprintf(fp, "%g %g lineto\n", faces[f]->c[i][0], faces[f]->c[i][1]);
    }
    fprintf(fp, "closepath\n");
    /* fprintf(fp, "gsave\n");*/
    if(faces[f]->width == DASHED) {
      fprintf(fp, "%d setlinewidth %d setlinecap %d setlinejoin 3.863693",
	      DASWTH, LINCAP, LINJIN);
      /* fprintf(fp,
	 " setmiterlimit [0 0 0 1]setcolor [2 4] 0 setdash {stroke}fp\n");*/
      fprintf(fp,
	      " setmiterlimit [0 0 0 1]setcolor [2 4] 0 setdash stroke\n");
    }
    else {
      /* fprintf(fp, "%g setlinewidth %d setlinecap %d setlinejoin 3.863693",
	 faces[f]->width, LINCAP, LINJIN); */
      fprintf(fp, "%g setlinewidth %d setlinecap %d setlinejoin ",
	      faces[f]->width, LINCAP, LINJIN);
      /* fprintf(fp, "[0 0 0 1]setcolor  {stroke}fp\ngrestore\n"); */
      /* fprintf(fp, "[0 0 0 1]setcolor  stroke\n");*/
      fprintf(fp, " 0 setgray  stroke\n");
      /* fprintf(fp, "grestore\n");*/
    }
    if(n_ == TRUE) numberFace(faces[f], fp);
  }

  dumpLines(fp, lines, numlines);

  /* if this is just to check placement, number the faces */
  if(n_ == TRUE) {
    /* numberFaces(faces, numfaces, fp); */
    numberLines(lines, numlines, fp);
    /*fprintf(stderr, "Faces and lines numbered\n");*/
  }

  /* if fills were not included, say so
  if(f_) fprintf(stderr, "Face fills not written to ps file\n"); */

  /* print shading key if not disabled and charge density info was inputed */
  if(q_ && !rk_ && !m_)
      dump_shading_key(fp, KEYBLKS, KEYPREC, KEYFONT, use_density);

  /* print footer */
  if(c_ == TRUE) {			/* print command line if asked for */
    for(f = 0, linein[0] = '\0'; f < argc; f++) {
      strcat(linein, argv[f]);
      strcat(linein, " ");
    }
    dump_line_as_ps(fp, linein, OFFSETX+2*CMDFONT, IMAGEY-2*CMDFONT, CMDFONT);
    /*fprintf(stderr, "Command line printed\n");*/
  }

  fprintf(fp, "vmr\nend  %% FreeHandDict\n");
  if(s_ == FALSE) {
    fprintf(fp, "showpage\n");
    /*fprintf(stderr, "Showpage inserted\n");*/
  }
  fprintf(fp, "%%%%EndDocument: _\n");
}

