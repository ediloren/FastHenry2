/* this will read the input file */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "induct.h"

#define MAXLINE 1000
#define XX 0
#define YY 1
#define ZZ 2

typedef struct Defaults {
  int isx,
      isy,
      isz,
      ish,
      isw,
      issigma,
      ishinc,
      iswinc,
      isrw,           /* width ratio */
      isrh;           /* height ratio */

  double x,y,z,h,w,sigma;
  int hinc, winc;
  double rw, rh;
} DEFAULTS;

char *getaline();

static DEFAULTS defaults = {0, 0 , 0, 0, 0, 1,     1 ,1,  1 ,  1,
			    0, 0 , 0, 0, 0, 5.8e7, 1, 1, 2.0, 2.0};

static double units = 1.0;

int plane_count;                     /* CMS 7/1/92 */
int sizenodes = 0;
int sizesegs = 0;

int readGeom(fp, indsys)
FILE *fp;
SYS *indsys;
{
  char *line;
  int end, error;
  NODES *node;
  SEGMENT *seg;
  GROUNDPLANE *plane;               /* CMS 7/2/92 */

  plane_count = indsys->num_extern = 0;  /* CMS 7/2/92 */
  indsys->num_segs = indsys->num_nodes = 0;

  /* read title */
  line = getaline(fp);
  printf("Title:\n%s\n", line);
  indsys->title = (char *)MattAlloc(strlen(line)+1,sizeof(char));
  strcpy(indsys->title,line);

  end = 0;
  error = 0;
  while (!end ){

    line = getaline(fp);
    tolowercase(line);
    switch (line[0]) {
    case '*':
      break;
    case '.':
      end = dodot(line + 1, indsys);
      break;
    case 'n':
      end = addnode(line,indsys, &node, NORMAL);
      break;
    case 'e':
      end = addseg(line, indsys, NORMAL, &seg);    /* CMS 8/21/92 normal seg */
      break;
    case 'g':                                      /* CMS 7/1/92 */
      end = addgroundplane(line,indsys, &plane);
      break;
    default:
      end = nothing(line);
      break;
    }
    if (end == 1) {
      error = 1;
      end = 0;
    }
  }

  indsys->num_planes = plane_count;           /* CMS 7/2/92 */

  if (end == 2) end = error;
  return end;
}

dodot(line,indsys)
char *line;
SYS *indsys;
{
  int end;

  if (strncasecmp("uni",line, 3) == 0)
    end = changeunits(line, indsys);
  else if (strncasecmp("ext",line, 3) == 0)
    end = addexternal(line, indsys);
  else if (strncasecmp("fre",line, 3) == 0)
    end = choosefreqs(line, indsys);
  else if (strncasecmp("equ",line, 3) == 0)
    end = equivnodes(line, indsys);
  else if (strncasecmp("def",line, 3) == 0)
    end = dodefault(line, indsys);
  else if (strncasecmp("end",line, 3) == 0)
    end = 2;
  else {
    printf("unrecognized . command\n.%s\n",line);
    return 1;
  }
  return end;
}

changeunits(line, indsys)
char *line;
SYS *indsys;
{
  char unitname[20];

  if (sscanf(line, "%*s %s", unitname) != 1) {
    printf("couldn't read units line\n.%s\n");
    return 1;
  }
  else if (strncasecmp("mil",unitname, 3) == 0)
    units = 2.54e-5;                     /* mils to meters */
  else if (strncasecmp("in",unitname, 2) == 0)
    units = 0.0254;                      /* inches to meters */

  else if (strncasecmp("um",unitname, 2) == 0)
    units = 1e-6;                        /* microns to meters */
  else if (strncasecmp("mm",unitname, 2) == 0)
    units = 1e-3;                      /* millimeters to meters*/
  else if (strncasecmp("cm",unitname, 2) == 0)
    units = 1e-2;                      /* centimeters to m */
  else if (strncasecmp("m",unitname, 1) == 0)
    units = 1.0;                      /* meters to meters */
  else if (strncasecmp("km",unitname, 1) == 0)
    units = 1e3;                      /* kilometers to meters */
  else {
    printf("unrecognizable unit name \n");
    return (1);
  }

  indsys->units = units;

  printf("all lengths multiplied by %lg to convert to meters\n",units);
  return 0;
}

addexternal(line, indsys)
char *line;
SYS *indsys;
{
  int skip, i;
  NODES *node[2];
  SEGMENT *seg;
  static char name[80];
  static char names[2][80];
  static char portname[80];
  SPATH *pathtail;
  NPATH *npathtail;
  PSEUDO_SEG *Vsource;
  int err;
  EXTERNAL *ext;

  /* new variables CMS 6/1/92 ---------------------------------------*/
  int nodes1, nodes2, inconnectread, outconnectread, j;
  int deco1, deco2, counter1, counter2;
  int beg0, beg1, end0, end1, done;
  double dumb, temp1, temp2, temp3, input[3], output[3];
  char *templine;
  GROUNDPLANE *plane;
  NODES *outnode, *innode, *finode, *pnode;
  SPATH *pathpointer;

  int equivnodes();

  void findrefnodes();                    /* functions from addgroundplane.c */
  void makenpath();
  NODES *find_nearest_node();
  /*--------------------------------------------------------------*/

  templine = (char *)MattAlloc(100,sizeof(char));           /* CMS 9/4/92 */

  err = 0;

  if(sscanf(line, "%*s%n", &skip) != 0){
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  /* read the nodes which will be external */
  for(i = 0; i < 2; i++) {
    if(sscanf(line, "%s%n",name,&skip) != 1){
      printf("no node name to read for .external %s\n",name);
      return 1;
    }
    line += skip;

    strcpy(names[i], name);

    /* search for the named node in the list */
    node[i] = get_node_from_name(name, indsys);

    if(node[i] == NULL){
      printf("addexternal: No node read in yet named %s\n",name);
      return 1;
    }

    node[i] = getrealnode(node[i]);
  }

  if(sscanf(line, "%s%n",portname,&skip) != 1){
    portname[0] = '\0';
  }
  else {
    line += skip;
  }

  if ((ext = get_external_from_portname(portname,indsys)) != NULL) {
    fprintf(stderr, "Error: Cannot name port between %s and %s as %s.\n \
The name is already used for port between %s and %s\n",
	    names[0],names[1],portname,ext->name1,ext->name2);
    err = 1;
  }

  if (node[0] == node[1]) {
    fprintf(stderr, "Error in .external: Nodes %s and %s are the same node. (.equiv'ed maybe?)\n", names[0], names[1]);
    if (is_gp_node(node[0]))
      fprintf(stderr,"  Nodes could be coincident because of coarse ground plane.\n");
    err = 1;
  }

  Vsource = make_pseudo_seg(node[0], node[1], EXTERNTYPE);
/*
  indsys->externals = add_to_external_list(make_external(Vsource,
							 indsys->num_extern),
					   indsys->externals);
*/
  indsys->externals
    = add_to_external_list(make_external(Vsource, indsys->num_extern, names[0],
					 names[1], portname),
			   indsys->externals);
  indsys->num_extern += 1;

  return err;
}

choosefreqs(line, indsys)
char *line;
SYS *indsys;
{
  int skip;
  double dumb;
  int fminread = 0, fmaxread = 0, logofstepread = 0;
  double fmin, fmax, logofstep;

  if (sscanf(line,"%*s%n",&skip) != 0) {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  while(notblankline(line)) {
    if (sscanf(line," fmin = %lf%n",&dumb,&skip) == 1) {
      fminread = 1;
      fmin = dumb;
      line += skip;
    }
    else if (sscanf(line," fmax = %lf%n",&dumb,&skip) == 1) {
      fmaxread = 1;
      fmax = dumb;
      line += skip;
    }
    else if (sscanf(line," ndec = %lf%n",&dumb,&skip) == 1) {
      if (dumb == 0) {
        dumb = 0.01;
        printf("Warning: ndec == 0.  Resetting to a very small number instead (%lg)\n",
               dumb);
      }
      logofstep = 1.0/dumb;
      logofstepread = 1;
      line += skip;
    }
    else {
      printf("don't know what piece this is in freq:%s\n",line);
      return 1;
    }
  }  /* end while(notblank... */

  if (fminread == 0) {
      printf("fmin not given for frequency range\n");
      return 1;
  }
  if (fmaxread == 0) {
      printf("fmax not given for frequency range\n");
      return 1;
  }
  if (logofstepread == 0) {
    logofstep = 1.0;
  }

 indsys->fmin = fmin;
 indsys->fmax = fmax;
 indsys->logofstep = logofstep;

 return 0;

}
old_equivnodes(line, indsys)
char *line;
SYS *indsys;
{

  int skip, i;
  char name1[80], name2[80];
  NODES *realnode, *equivnode;

  if (sscanf(line,"%*s%n",&skip) != 0) {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  if (sscanf(line,"%s%n",name1,&skip) != 1) {
    printf("No first node name on equiv line\n");
    return 1;
  }
  line += skip;

  if (sscanf(line,"%s%n",name2,&skip) != 1) {
    printf("No second node name on equiv line\n");
    return 1;
  }
  line += skip;

  if (sscanf(line,"%s%n",name1,&skip) == 1) {
    printf("Only allowed one equiv per line!\n");
    return 1;
  }

  realnode = indsys->nodes;
  while( realnode != NULL && strcmp(name1,realnode->name) != 0)
    realnode = realnode->next;

  if (realnode == NULL) {
    printf("No node read in yet named %s\n",name1);
    return 1;
  }

  equivnode = indsys->nodes;
  while( equivnode != NULL && strcmp(name2,equivnode->name) != 0)
    equivnode = equivnode->next;

  if (equivnode == NULL) {
    printf("No node read in yet named %s\n",name1);
    return 1;
  }
  else
    equivnode->equiv = realnode;

  return 0;
}

dodefault(line)
char *line;
{
  int skip;
  double dumb;
  int dumbi;

  if (sscanf(line,"%*s%n",&skip) != 0) {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  while(notblankline(line)) {
    if (sscanf(line," x = %lf%n",&dumb,&skip) == 1) {
      defaults.isx = 1;
      defaults.x = units*dumb;
      line += skip;
    }
    else if (sscanf(line," y = %lf%n",&dumb,&skip) == 1) {
      defaults.isy = 1;
      defaults.y = units*dumb;
      line += skip;
    }
    else if (sscanf(line," z = %lf%n",&dumb,&skip) == 1) {
      defaults.isz = 1;
      defaults.z = units*dumb;
      line += skip;
    }
    else if (sscanf(line," h = %lf%n",&dumb,&skip) == 1) {
      defaults.ish = 1;
      defaults.h = units*dumb;
      line += skip;
    }
    else if (sscanf(line," w = %lf%n",&dumb,&skip) == 1) {
      defaults.isw = 1;
      defaults.w = units*dumb;
      line += skip;
    }
    else if (sscanf(line," rho = %lf%n",&dumb,&skip) == 1) {
      defaults.issigma = 1;
      defaults.sigma = (1.0/dumb)/units;
      line += skip;
    }
    else if (sscanf(line," sigma = %lf%n",&dumb,&skip) == 1) {
      defaults.issigma = 1;
      defaults.sigma = dumb/units;
      line += skip;
    }
    else if (sscanf(line," nhinc = %d%n",&dumbi,&skip) == 1) {
      defaults.ishinc = 1;
      defaults.hinc = dumbi;
      line += skip;
    }
    else if (sscanf(line," nwinc = %d%n",&dumbi,&skip) == 1) {
      defaults.iswinc = 1;
      defaults.winc = dumbi;
      line += skip;
    }
    else if (sscanf(line," rw = %lf%n",&dumb,&skip) == 1) {
      defaults.isrw = 1;
      defaults.rw = dumb;
      line += skip;
    }
    else if (sscanf(line," rh = %lf%n",&dumb,&skip) == 1) {
      defaults.isrh = 1;
      defaults.rh = dumb;
      line += skip;
    }
    else {
      printf("don't know what piece this is in default:%s\n",line);
      return 1;
    }
  }  /* end while(notblank... */

  return 0;
}

addnode(line,indsys,retnode, type)
char *line;
SYS *indsys;
NODES **retnode;
int type;
{
  double dumb;
  int skip;
  int xread = 0, yread = 0, zread = 0;
  NODES *node;
  static char name[80];
  NODES *makenode();
  double nodex, nodey, nodez;

  /* read name */
  sscanf(line,"%s%n",name,&skip);
  line += skip;

    while(notblankline(line)) {
      if (sscanf(line," x = %lf%n",&dumb,&skip) == 1) {
	nodex = units*dumb;
	line += skip;
	xread = 1;
      }
      else if (sscanf(line," y = %lf%n",&dumb,&skip) == 1) {
	nodey = units*dumb;
	line += skip;
	yread = 1;
      }
      else if (sscanf(line," z = %lf%n",&dumb,&skip) == 1) {
	nodez = units*dumb;
	line += skip;
	zread = 1;
      }
      else {
	printf("don't know what piece this is in node %s: %s\n",name,line);
	return (1);
      }
    } /* while not blank */

  if (zread == 0) {
    if (defaults.isz == 1)
      nodez = defaults.z;
    else {
      printf("Z not given for node %s\n",name);
      return 1;
    }
  }
  if (yread == 0) {
    if (defaults.isy == 1)
      nodey = defaults.y;
    else {
      printf("Y not given for node %s\n",name);
      return 1;
    }
  }
  if (xread == 0) {
    if (defaults.isx == 1)
      nodex = defaults.x;
    else {
      printf("X not given for node %s\n",name);
      return 1;
    }
  }

  /* make sure node isn't already defined */
  if (get_node_from_name(name, indsys) != NULL) {
    fprintf(stderr, "Node %s already defined. Cannot redefine it\n", name);
    exit(1);
  }

  node = makenode(name, indsys->num_nodes, nodex, nodey, nodez, type, indsys);

#if 1==0
  /* add node to linked list of nodes */
  if (indsys->nodes == NULL) {
    indsys->nodes = node;
    indsys->endnode = node;
  }
  else {
    indsys->endnode->next = node;
    indsys->endnode = node;
  }
#endif
  *retnode = node;
  return 0;
}

NODES *makenode(name, number, x, y, z, type, indsys)
char *name;
int number;
double x, y, z;
int type;
SYS *indsys;
{
  NODES *node;

  node = (NODES *)MattAlloc(1, sizeof(NODES));
  node->name = (char *)MattAlloc( (strlen(name) + 1), sizeof(char));
  strcpy(node->name, name);
  node->number = number;
  node->index = -1;     /* to be assigned later if needed */
  node->x = x;
  node->y = y;
  node->z = z;
  node->equiv = node;
  node->to_end = NULL;
  node->connected_segs = NULL;
  node->num_to_end = 0;  /* note: no longer accurate until 'fillM' is called */
  node->gp = NULL;
  node->treeptr = NULL;
  node->level = -1;
  node->pred.segp = NULL;
  node->next = NULL;

  node->gp_node = NULL;  /* nonuni gp only */

  node->type = type;
  if (is_gp_node(node))
    node->examined = 1;
  else
    node->examined = 0;

  if (indsys != NULL) {
    /* add node to linked list of nodes */
    if (indsys->nodes == NULL) {
      indsys->nodes = node;
      indsys->endnode = node;
    }
    else {
      indsys->endnode->next = node;
      indsys->endnode = node;
    }
    indsys->num_nodes++;
  }

  return node;
}

addseg(line, indsys, type, retseg)
char *line;
SYS *indsys;
int type;                     /* CMS 8/21/92 -- type of thing the seg is in */
SEGMENT **retseg;
{
  double dumb, *tmp;
  int skip, i, j, dumbi;
  int hread, wread, hincread, wincread, sigmaread, wxread, wyread, wzread,
      rwread, rhread;
  SEGMENT *seg;
  static char name[80], n1[80], n2[80];
  NODES *anode;

   char *segname;
   double *widthdir;   /*if width is not || to x-y plane and perpendicular to*/
                       /* the length, then this is 3 element vector in       */
                       /* in the direction of width*/
   int number;         /* an arbitrary number for the segment */
   double length;
   double area;        /* area of cross section */
   double width, height;  /*width and height to cross section */
   int hinc, winc;             /* number of filament divisions in each dir */
   NODES *node[2];                /* nodes at the ends */
   double sigma;               /* conductivity */
   int num_fils;               /* hinc*winc */
   double rh,rw;           /* for assignFil(), ratio of adjacent fils */
   FILAMENT *filaments;        /* this segment's filaments */
   struct pathlist *loops;   /* loops in which this segment is a member */

  hread = wread = hincread = wincread = sigmaread = wxread
    = wyread = wzread = rwread = rhread = 0;

  widthdir = NULL;

  /* read name */
  sscanf(line,"%s%n",name,&skip);
  segname = (char *)MattAlloc( (strlen(name) + 1), sizeof(char));
  strcpy(segname, name);
  line += skip;

  /* read nodename one */
  for (j = 0; j < 2; j++) {
    sscanf(line,"%s%n",name,&skip);
    line += skip;

    anode = get_node_from_name(name, indsys);

    if (anode == NULL) {
      printf("No node read in yet named %s\n",name);
      return 1;
    }
    else {
      node[j] = anode;
    }
  }

    while(notblankline(line)) {
      if (sscanf(line," h = %lf%n",&dumb,&skip) == 1) {
	height = units*dumb;
	line += skip;
	hread = 1;
      }
      else if (sscanf(line," w = %lf%n",&dumb,&skip) == 1) {
	width = units*dumb;
	line += skip;
	wread = 1;
      }
      else if (sscanf(line," sigma = %lf%n",&dumb,&skip) == 1) {
	sigma = dumb/units;
	line += skip;
	sigmaread = 1;
      }
      else if (sscanf(line," rho = %lf%n",&dumb,&skip) == 1) {
        sigma = (1.0/dumb)/units;
	line += skip;
	sigmaread = 1;
      }
      else if (sscanf(line," nhinc = %d%n",&dumbi,&skip) == 1) {
	hinc = dumbi;
	line += skip;
	hincread = 1;
      }
      else if (sscanf(line," nwinc = %d%n",&dumbi,&skip) == 1) {
	winc = dumbi;
	line += skip;
	wincread = 1;
      }
      else if (sscanf(line," wx = %lf%n",&dumb,&skip) == 1) {
	if (widthdir == NULL)
	  widthdir = (double *)MattAlloc(3,sizeof(double));
	widthdir[XX] = units*dumb;
	line += skip;
	wxread = 1;
      }
      else if (sscanf(line," wy = %lf%n",&dumb,&skip) == 1) {
	if (widthdir == NULL)
	  widthdir = (double *)MattAlloc(3,sizeof(double));
	widthdir[YY] = units*dumb;
	line += skip;
	wyread = 1;
      }
      else if (sscanf(line," wz = %lf%n",&dumb,&skip) == 1) {
	if (widthdir == NULL)
	  widthdir = (double *)MattAlloc(3,sizeof(double));
	widthdir[ZZ] = units*dumb;
	line += skip;
	wzread = 1;
      }
      else if (sscanf(line," rw = %lf%n",&dumb,&skip) == 1) {
	rw = dumb;
	line += skip;
	rwread = 1;
      }
      else if (sscanf(line," rh = %lf%n",&dumb,&skip) == 1) {
	rh = dumb;
	line += skip;
	rhread = 1;
      }
      else {
       printf("don't know what piece this is in seg %s: %s\n",segname,line);
       return (1);
      }
    } /* while not blank */

  if (hread == 0) {
    if (defaults.ish == 1)
      height = defaults.h;
    else {
      printf("H not given for seg %s\n",segname);
      return 1;
    }
  }
  if (wread == 0) {
    if (defaults.isw == 1)
      width = defaults.w;
    else {
      printf("W not given for seg %s\n",segname);
      return 1;
    }
  }
  if (sigmaread == 0) {
    if (defaults.issigma == 1)
      sigma = defaults.sigma;
    else {
      printf("Sigma or rho not given for seg %s\n",segname);
      return 1;
    }
  }
  if (hincread == 0) {
    if (defaults.ishinc == 1)
      hinc = defaults.hinc;
    else
      hinc = 1;
  }
  if (wincread == 0) {
    if (defaults.iswinc == 1)
      winc = defaults.winc;
    else
      winc = 1;
  }
  if (rwread == 0) {
    if (defaults.isrw == 1)
      rw = defaults.rw;
    else if (winc != 1) {
      printf("ratio of adjacent filaments along width not given for seg: %s\n",
	     segname);
      return 1;
    }
  }

  if (rhread == 0) {
    if (defaults.isrh == 1)
      rh = defaults.rh;
    else if (hinc != 1) {
     printf("ratio of adjacent filaments along height not given for seg: %s\n",
	    segname);
     return 1;
    }
  }

  if (rw < 1.0) {
    printf("Error: ratio of adjacent fils less than 1.0 for seg: %s\n",
	   segname);
    return 1;
  }
  if (rh < 1.0) {
    printf("Error: ratio of adjacent fils less than 1.0 for seg: %s\n",
	   segname);
    return 1;
  }

  if ( widthdir != NULL) {
    if ( (wxread + wyread + wzread) != 3) {
      printf("not all of wx wy and wz specified for seg %s\n",segname);
      return 1;
    }
    tmp = widthdir;
    dumb = sqrt(tmp[XX]*tmp[XX] + tmp[YY]*tmp[YY] + tmp[ZZ]*tmp[ZZ]);
    tmp[XX] /= dumb;
    tmp[YY] /= dumb;
    tmp[ZZ] /= dumb;
  }

  if (is_normal(node[0]) != TRUE || is_normal(node[1]) != TRUE) {
    printf("Error: Segment %s has a ground plane node as one of its ends.\n",
	   segname);
    printf("   Define a new normal node as its end and then '.equiv' it to the ground plane\n");
    return 1;
  }

  seg = makeseg(segname, node[0], node[1], height, width, sigma, hinc, winc,
		rh, rw, widthdir, indsys->num_segs, type, indsys);

#if 1 == 0
  /* add segment to linked list */
  if (indsys->segment == NULL) {
    indsys->segment = seg;
    indsys->endseg = seg;
  }
  else {
    indsys->endseg->next = seg;
    indsys->endseg = seg;
  }
#endif

  *retseg = seg;

  return 0;
}

SEGMENT *makeseg(name, node0, node1, height, width, sigma, hinc, winc,
		r_height, r_width, widthdir, number, type, indsys)
   char *name;
   double *widthdir;  /*if width is not || to x-y plane and perpendicular to*/
                       /* the length, then this is 3 element vector in       */
                       /* in the direction of width*/
   int number;         /* an arbitrary number for the segment */
   int type;    /* CMS 8/21/92 -- type of structure the segment is in */
   double width, height;  /*width and height to cross section */
   int hinc, winc;             /* number of filament divisions in each dir */
   NODES *node0, *node1;                /* nodes at the ends */
   double sigma;               /* conductivity */
   double r_height, r_width; /* ratio of adjacent fils for assignFil() */
   SYS *indsys;         /* nonNULL if we are to add to system linked list */
{

  SEGMENT *seg;

  seg = (SEGMENT *)MattAlloc(1, sizeof(SEGMENT));

  seg->name = (char *)MattAlloc(strlen(name) + 1, sizeof(char));
  strcpy(seg->name, name);
  seg->widthdir = widthdir;
  seg->number = number ;
  seg->type = type ;
  seg->width = width ;
  seg->height = height ;
  seg->hinc = hinc;
  seg->winc = winc;
  seg->r_height = r_height;
  seg->r_width = r_width;
  seg->sigma = sigma;
  seg->node[0] = node0;
  seg->node[1] = node1;
  add_to_connected_segs(node0, seg, NULL);
  add_to_connected_segs(node1, seg, NULL);

  seg->length = sqrt( (node0->x - node1->x)*(node0->x - node1->x)
		     + (node0->y - node1->y)*(node0->y - node1->y)
		     + (node0->z - node1->z)*(node0->z - node1->z)
		     );
  seg->area = width*height;

  seg->filaments = NULL;
  seg->num_fils = 0;

  seg->loops = NULL;
  seg->is_deleted = 0;
  /*  seg->gp_node[0] = seg->gp_node[1] = NULL; */

  /*seg->table = NULL;*/

  seg->next = NULL;

  if (indsys != NULL) {
    /* add segment to linked list */
    if (indsys->segment == NULL) {
      indsys->segment = seg;
      indsys->endseg = seg;
    }
    else {
      indsys->endseg->next = seg;
      indsys->endseg = seg;
    }
    indsys->num_segs++;
  }

  return seg;
}

/*-------------------------------------------------------------------------------*/
/*                      CMS code addition to support groundplanes                */
/*-------------------------------------------------------------------------------*/
  /* functions called from addgroundplane.c */
  void fillgrids();
  void doincrment();
  void dounitvector();
  void make_nodelist();
  int checkmiddlepoint();
  int checkplaneformula();
  double findsegmentwidth();
  double lengthof2();
  double find_coordinate();
  SPATH *path_through_gp();


addgroundplane(line, indsys, retplane)
char *line;
SYS *indsys;
GROUNDPLANE **retplane;
{

  /* read in variables */
  static char name[80], filename[255];
  int xred = 0, yred = 0, zred = 0, segheiread = 0, rhread = 0;
  int segread1 = 0, segread2 = 0, inred = 0, outred = 0, filenameread = 0;
  int skip;
  double input[3], output[3];

  NODELIST *list_of_nodes;       /* linked list of read in nodes */
  NODELIST *listpointer;         /* pointer to a single element in the list */
  HoleList *list_of_holes, *holep; /* linked lists of holes to be made. */
  ContactList *list_of_contacts, *contactp; /* linked lists of contacts. */
  static char nodename[80];
  char *coordinate_string;
  double dumbx, dumby, dumbz;

  /* plane calculation variables */
  int nodes1, nodes2, checksum;
  double x, y, z, seghei, rh;
  double wx, wy, wz, dx1, dy1, dz1, dx2, dy2, dz2;
  double xinit, yinit, zinit;
  double segwid1 = -1;
  double segwid2 = -1;
  double segfull1, segfull2;

  /* layout variables */
  char *templine;
  int i, j, b0, b1, e0, e1;
  int signofelem;
  int o1 = 0, mid = 1, o2 = 2;
  int a = 0, b = 0;
  int err;

  /* temporary storage structures and variables */
  NODES *tempnode;
  SEGMENT *tempseg;
  PSEUDO_SEG *p_seg;
  GROUNDPLANE *grndp, *onegp;
  SPATH *tempath, *pathp;
  double dumb, dontcare = 2.0;
  int dummy1, dummy2;
  double xt[4], yt[4], zt[4];
  int dumbi, tseg1, tseg2;
  double get_perimeter();

  /* file for nonuniform plane hierarchy */
  FILE *nonuni_fp;

  double relx, rely, relz;  /* relative coords for plane.  MK 10/92 */
  double *widthdir;  /* reallocated for each seg */

  relx = rely = relz = 0;

  err = 0;

  /* allocate space for the groundplane */
  grndp = (GROUNDPLANE *)MattAlloc(1, sizeof(GROUNDPLANE));
  grndp->indsys = indsys;
  grndp->usernodes = NULL;
  grndp->fake_seg_list = NULL;
  grndp->sigma = defaults.sigma;
  grndp->rh = defaults.rh;
  grndp->hinc = 1;
  grndp->filename = NULL;
  grndp->nonuni = NULL;
  coordinate_string = (char *)MattAlloc(1000, sizeof(char));
  list_of_nodes = NULL;
  list_of_holes = holep = NULL;
  list_of_contacts = contactp = NULL;

  /* read in the name of the groundplane */
  if (sscanf(line, "%s%n",name,&skip) != 1) {
    fprintf(stderr, "addgroundplane: hey, no fair\n");
    exit(1);
  }

  grndp->name = (char *)MattAlloc(strlen(name) + 1, sizeof(char));
  strcpy(grndp->name, name);
  line += skip;

  /* read in groundplane specifications */
  while(notblankline(line)){
    if(sscanf (line, " x1 = %lf%n",&dumb, &skip) == 1){
      grndp->x[0] = xt[0] = dumb*units;
      line += skip;
      xred++;
    }
    else if(sscanf (line, " y1 = %lf%n",&dumb, &skip) == 1){
      grndp->y[0] = yt[0] = dumb*units;
      line += skip;
      yred++;
    }
    else if(sscanf (line, " z1 = %lf%n",&dumb, &skip) == 1){
      grndp->z[0] = zt[0] = dumb*units;
      line += skip;
      zred++;
    }
    else if(sscanf (line, " x2 = %lf%n",&dumb, &skip) == 1){
      grndp->x[1] = xt[1] = dumb*units;
      line += skip;
      xred++;
    }
    else if(sscanf (line, " y2 = %lf%n",&dumb, &skip) == 1){
      grndp->y[1] = yt[1] = dumb*units;
      line += skip;
      yred++;
    }
    else if(sscanf (line, " z2 = %lf%n",&dumb, &skip) == 1){
      grndp->z[1] = zt[1] = dumb*units;
      line += skip;
      zred++;
    }
    else if(sscanf (line, " x3 = %lf%n",&dumb, &skip) == 1){
      grndp->x[2] = xt[2] = dumb*units;
      line += skip;
      xred++;
    }
    else if(sscanf (line, " y3 = %lf%n",&dumb,&skip) == 1){
      grndp->y[2] = yt[2] = dumb*units;
      line += skip;
      yred++;
    }
    else if(sscanf (line, " z3 = %lf%n",&dumb,&skip) == 1){
      grndp->z[2] = zt[2] = dumb*units;
      line += skip;
      zred++;
    }
    else if(sscanf (line, " thick = %lf%n",&dumb,&skip) == 1){
      seghei = dumb*units;
      line += skip;
      segheiread = 1;
    }
    else if(sscanf (line, " segwid1 = %lf%n",&dumb,&skip) == 1){
      segwid1 = dumb*units;
      line += skip;
    }
    else if(sscanf (line, " segwid2 = %lf%n",&dumb,&skip) == 1){
      segwid2 = dumb*units;
      line += skip;
    }
    else if(sscanf (line, " relx = %lf%n",&dumb,&skip) == 1){
      relx = dumb*units;
      line += skip;
    }
    else if(sscanf (line, " rely = %lf%n",&dumb,&skip) == 1){
      rely = dumb*units;
      line += skip;
    }
    else if(sscanf (line, " relz = %lf%n",&dumb,&skip) == 1){
      relz = dumb*units;
      line += skip;
    }
    else if(sscanf (line, " seg1 = %d%n",&dumbi,&skip) == 1){
      grndp->seg1 = dumbi;
      line += skip;
      segread1 = 1;
    }
    else if(sscanf (line, " seg2 = %d%n",&dumbi,&skip) == 1){
      grndp->seg2 = dumbi;
      line += skip;
      segread2 = 1;
    }
    else if (sscanf(line," sigma = %lf%n",&dumb,&skip) == 1) {
      grndp->sigma = dumb/units;
      line += skip;
    }
    else if (sscanf(line," rho = %lf%n",&dumb,&skip) == 1) {
      grndp->sigma = (1.0/dumb)/units;
      line += skip;
    }
    else if (sscanf(line," nhinc = %d%n",&dumbi,&skip) == 1) {
      grndp->hinc = dumbi;
      line += skip;
    }
    else if (sscanf(line," rh = %lf%n",&dumb,&skip) == 1) {
      grndp->rh = dumb;
      line += skip;
    }
    else if (sscanf(line," file = %s%n",filename,&skip) == 1) {
      grndp->filename = (char *)MattAlloc(strlen(filename) + 1, sizeof(char));
      strcpy(grndp->filename, filename);
      filenameread = 1;
      line += skip;
    }
    else if (is_next_word("hole",line) == TRUE) {
      list_of_holes = make_holelist(list_of_holes, line,
				    units, relx, rely, relz, &skip);
      line += skip;
    }
    else if (is_next_word("contact",line) == TRUE) {
      list_of_contacts = make_contactlist(list_of_contacts, line,
				    units, relx, rely, relz, &skip);
      line += skip;
    }
    /* Read in a ground plane node reference */
    else if(sscanf (line, "%s %s%n", nodename, coordinate_string, &skip) == 2
	    && ((nodename[0] == 'n') && (coordinate_string[0] == '(') )){
      /* allocate space for the list_of_nodes */
      if(list_of_nodes == NULL){         /* first one in the list */
	list_of_nodes = (NODELIST *)Gmalloc(1* sizeof(NODELIST));
	list_of_nodes->name = (char *)Gmalloc( (strlen(nodename) + 1)* sizeof(char));
	listpointer = list_of_nodes;
	listpointer->next = NULL;
      } else {
	listpointer->next = (NODELIST *)Gmalloc(1* sizeof(NODELIST));
	listpointer = listpointer->next;
	listpointer->name = (char *)Gmalloc( (strlen(nodename) + 1)* sizeof(char));
	listpointer->next = NULL;
      }

      if(sscanf (coordinate_string, "( %lf , %lf , %lf )",&dumbx,
		 &dumby, &dumbz) == 3){
	/* all the coordinates are given */
      }
      else if(sscanf (coordinate_string, "( , %lf, %lf )", &dumby, &dumbz) == 2){
	/* missing the x-coordinate */
	dumbx = 0; /*find_coordinate(grndp, 0.0,  dumby, dumbz, 0);*/
      }
      else if(sscanf (coordinate_string, "( %lf ,  , %lf)", &dumbx, &dumbz) == 2){
	/* missing the y-coordinate */
	dumby = 0 ;  /* find_coordinate(grndp, dumbx, 0.0, dumbz, 1); */
      }
      else if(sscanf (coordinate_string, "( %lf , %lf , )", &dumbx, &dumby) == 2){
	/* missing the z-coordinate */
	dumbz = 0; /*find_coordinate(grndp, dumbx, dumby, 0.0, 2);*/
      } else {
	printf("Error: coordinates do not match correctly %s.\n",coordinate_string);
	printf("  Are there spaces in the string?\n");
	exit(1);
      }
      make_nodelist(listpointer, nodename, dumbx*units, dumby*units, dumbz*units);
      line += skip;
    }
    else {
      printf("don't know which piece this is in groundplane %s: %s\n",name,
	     line);
      return(1);
    }
  }/* while not a blank line */

  /* do checks to verify the line was read correctly */
  if(xred != 3) {
    printf(" x coordinate not given for plane\n");
    for (i = 0; i < 3; i++){
      printf("x-coordinate, point %d: %lg\n",i,xt[i]);
    }
    return 1;
  }

  if(yred != 3) {
    printf(" y coordinate not given for plane\n");
    for (i = 0; i < 3; i++){
      printf("y-coordinate, point %d: %lg\n",i,yt[i]);
    }
    return 1;
  }

  if(zred != 3) {
    printf(" z coordinate not given for plane\n");
    for (i = 0; i < 3; i++){
      printf("z-coordinate, point %d: %lg\n",i,zt[i]);
    }
    return 1;
  }

  if(segheiread == 0){
    printf("Thickness not given for plane %s\n",grndp->name);
    return 1;
  }

  /* Check discretization info */
  if (filenameread == 1 && (segread1 == 1 || segread2 == 1)) {
    printf("Can't specify both uniform plane discretization (seg1,seg2) and \
nonuniform discretization file %s\n",grndp->filename);
    return 1;
  }
  else if (filenameread == 1) {
    /* open file */
    if (strcmp(grndp->filename,"none") == 0)
      nonuni_fp = NULL;
    else
      if ( (nonuni_fp = fopen(grndp->filename,"r")) == NULL) {
	printf("Couldn't open nonuniform hierarchy file %s for plane %s\n",
	       grndp->filename, grndp->name);
	return 1;
      }
  }
  else {
    if(segread1 == 0){
      printf("Segments (seg1) not given for plane %s\n",grndp->name);
      return 1;
    }
    if(segread2 == 0){
      printf("Segments (seg2) not given for plane %s\n",grndp->name);
      return 1;
    }
  }

  if (grndp->rh < 1.0) {
    printf("Error: ratio of adjacent fils less than one for plane %s\n",
	   grndp->name);
    return 1;
  }

  /* end checks on reading line */

  /* checking point coordinates and their reference numbers */
  checksum = checkmiddlepoint(xt, yt, zt, o1, o2, mid);
  if(checksum != 1){
    printf("coordinates not set up correctly... not perpendicular\n");
    exit(1);
  }

  /* finding fourth corner point of the plane */
  grndp->x[3] = xt[3] = (xt[o2] + xt[o1] - xt[mid]);
  grndp->y[3] = yt[3] = (yt[o2] + yt[o1] - yt[mid]);
  grndp->z[3] = zt[3] = (zt[o2] + zt[o1] - zt[mid]);

  grndp->thick = seghei;

  /* length in "x" */
  grndp->length1 = mag(xt[o1] - xt[mid], yt[o1] - yt[mid], zt[o1] - zt[mid]);

  /* length in "y" */
  grndp->length2 = mag(xt[o2] - xt[mid], yt[o2] - yt[mid], zt[o2] - zt[mid]);

  grndp->list_of_contacts = list_of_contacts;

  if (filenameread == 0) {
    /***** it's a good old uniformly discretized plane ******/

    /*finding segment widths for both segment orientations */
    tseg1 = grndp->seg1;
    tseg2 = grndp->seg2;
    segfull1 = findsegmentwidth(xt, yt, zt, mid, o2, o1, tseg2);
    segfull2 = findsegmentwidth(xt, yt, zt, mid, o1, o2, tseg1);

    /* determine if segwid1 was specified and if it is ok */
    if (segwid1 < 0)
      segwid1 = segfull1;
    else
      if (segwid1 > segfull1) {
	printf("Warning: segwid1 greater than segment separation for plane %s\n",name);
	printf(" Using segwid1 = %lg\n", segfull1);
	segwid1 = segfull1;
      }

    /* determine if segwid2 was specified and if it is ok */
    if (segwid2 < 0)
      segwid2 = segfull2;
    else
      if (segwid2 > segfull2) {
	printf("Warning: segwid2 greater than segment separation for plane %s\n",name);
	printf(" Using segwid2 = %lg\n", segfull2);
	segwid2 = segfull2;
      }

    /* save some values */
    grndp->segwid1 = segwid1;
    grndp->segwid2 = segwid2;

    /* setup parameters for laying out nodes and segments */
    xinit = xt[mid];
    yinit = yt[mid];
    zinit = zt[mid];
    grndp->num_nodes1 = nodes1 = tseg1 + 1;
    grndp->num_nodes2 = nodes2 = tseg2 + 1;

    /* allocate the space for segments, nodes, and temporary variables */
    dumbi = 15 + ((tseg1 + tseg2) * 15);

    templine = (char *)MattAlloc(dumbi + 1, sizeof(char));
    grndp->pnodes = (NODES ***)MatrixAlloc(nodes1, nodes2, sizeof(NODES *));
    grndp->segs1 = (SEGMENT ***)MatrixAlloc(tseg1, nodes2, sizeof(SEGMENT *));
    grndp->segs2 = (SEGMENT ***)MatrixAlloc(nodes1, tseg2, sizeof(SEGMENT *));

    /* finding increments from (mid) to (o1) and from (mid) to (o2) */
    doincrement(xt[o1], yt[o1], zt[o1], xinit, yinit, zinit, tseg1, &dx1, &dy1, &dz1);
    doincrement(xt[o2], yt[o2], zt[o2], xinit, yinit, zinit, tseg2, &dx2, &dy2, &dz2);

    grndp->d1 = mag(dx1, dy1, dz1);
    grndp->d2 = mag(dx2, dy2, dz2);
    grndp->ux1 = dx1/grndp->d1;
    grndp->uy1 = dy1/grndp->d1;
    grndp->uz1 = dz1/grndp->d1;
    grndp->ux2 = dx2/grndp->d2;
    grndp->uy2 = dy2/grndp->d2;
    grndp->uz2 = dz2/grndp->d2;
    grndp->unitdiag = mag(dx1+dx2,dy1+dy2, dz1+dz2);

    /* lay out all the nodes and the segments from (mid) to (o1) */
    for(i = 0; i < nodes2; i++){
      x = xinit + (dx2 * i);
      y = yinit + (dy2 * i);
      z = zinit + (dz2 * i);

      /* verifying coordinates */
      checksum = checkplaneformula(xt, yt, zt, x, y, z, mid, o1, o2);

      sprintf(templine, "n%s0_%d",grndp->name, i);
      tempnode = makenode(templine, indsys->num_nodes, x, y, z, GPTYPE, indsys);
      tempnode->gp = grndp;	/* set plane pointer */
      tempnode->s1 = 0;
      tempnode->s2 = i;

      grndp->pnodes[0][i] = tempnode;

      for( j = 1; j < nodes1; j++){

	x = x + dx1;
	y = y + dy1;
	z = z + dz1;

	sprintf(templine, "n%s%d_%d",grndp->name, j, i);
	tempnode = makenode(templine, indsys->num_nodes, x, y, z,
			    GPTYPE, indsys);
	tempnode->gp = grndp;	/* set plane pointer */
	tempnode->s1 = j;
	tempnode->s2 = i;

	grndp->pnodes[j][i] = tempnode;

	/* layout the segments from (mid) to (o2) */

	/* verifying coordinates */
	checksum = checkplaneformula(xt, yt, zt, x, y, z, mid, o1, o2);

      }
    }

    /* Remove (mark) nodes that are part of holes */
    for(holep = list_of_holes; holep != NULL; holep = holep->next)
      make_holes(holep, grndp);

    /* finding width direction for segments along (mid -> o1) vector */
    /* [this direction is the unit vector from (mid) to (o2) */
    /* dounitvector(xt[o2] , yt[o2], zt[o2], xinit, yinit, zinit,
       &wx, &wy, &wz);*/
    /* I changed the order so that the height vector points in the same  */
    /*  direction for all segs in the plane.  Makes assignFil consistent */
    /*  and meshes tighter if nhinc > 1.   MK  8/93                      */
    dounitvector(xinit, yinit, zinit, xt[o2] , yt[o2], zt[o2], &wx, &wy, &wz);

    /* layout the segments from (mid) to (o1) */
    for(i = 0; i < nodes2; i++){
      for(j = 0; j < (nodes1 - 1); j++){

	/* Make sure nodes aren't part of a hole */
	if (!is_hole(grndp->pnodes[j][i]) && !is_hole(grndp->pnodes[j+1][i])) {
	  widthdir = (double *)MattAlloc(3, sizeof(double));
	  widthdir[XX] = wx; widthdir[YY] = wy; widthdir[ZZ] = wz;
	  sprintf(templine, "e1%s%d_%d",grndp->name, j, i);

	  tempseg = makeseg(templine, grndp->pnodes[j][i], grndp->pnodes[j+1][i],
			    seghei, segwid1, grndp->sigma, grndp->hinc, 1,
			    grndp->rh, dontcare, widthdir, indsys->num_segs,
			    GPTYPE, indsys);

	  grndp->segs1[j][i] = tempseg;
	}
	else
	  grndp->segs1[j][i] = NULL;
      }
    }

    /* finding the width direction of segments along the (mid - o2) vector */
    /* [this direction is the unit vector from (mid) to (o1) */
    dounitvector(xt[o1], yt[o1], zt[o1], xinit, yinit, zinit, &wx, &wy, &wz);

    /* layout the segments from (mid) to (o1) */
    for(i = 0; i < (nodes2 - 1); i++){
      for(j = 0; j < nodes1; j++){

	/* Make sure nodes aren't part of a hole */
	if (!is_hole(grndp->pnodes[j][i]) && !is_hole(grndp->pnodes[j][i+1])) {
	  widthdir = (double *)MattAlloc(3, sizeof(double));
	  widthdir[XX] = wx; widthdir[YY] = wy; widthdir[ZZ] = wz;
	  sprintf(templine, "e2%s%d_%d",grndp->name, j, i);
	  tempseg = makeseg(templine, grndp->pnodes[j][i], grndp->pnodes[j][i+1],
			    seghei, segwid2, grndp->sigma, grndp->hinc, 1,
			    grndp->rh, dontcare, widthdir, indsys->num_segs,
			    GPTYPE, indsys);

	  grndp->segs2[j][i] = tempseg;
	}
	else
	  grndp->segs2[j][i] = NULL;
      }
    }

    grndp->numesh = grndp->seg1 * grndp->seg2; /* number of meshes in plane */

  } /* if !is_nonuni_gp(), ie old gp code */
  else {
    /******** READ a NONUNIFORM plane ********/
    if (process_plane(grndp, nonuni_fp, indsys) != 0)
      return 1;

    templine = (char *)MattAlloc(100, sizeof(char));

  }

  /* equivalencing extra reference nodes to groundplane nodes */
  if(list_of_nodes != NULL){
    if (indsys->opts->debug == ON)
      printf("adding reference nodes...\n");

    for(i = 0, listpointer = list_of_nodes; listpointer != NULL;
	listpointer = listpointer->next, i++){

      if (!is_nonuni_gp(grndp)) {
	tempnode = find_nearest_gpnode(listpointer->x + relx,
				       listpointer->y + rely,
				       listpointer->z + relz,
				       grndp,&dummy1, &dummy2);
      }
      else {
	sprintf(templine, "n%s_pseudo_%d", grndp->name, i);
	tempnode = get_or_make_nearest_node(templine, indsys->num_nodes,
					    listpointer->x + relx,
			    listpointer->y + rely, listpointer->z + relz,
			    indsys, grndp->nonuni, grndp->usernodes);
      }

      if (is_hole(tempnode)) {
	fprintf(stderr,"Warning: ground plane node %s is in a hole.\n \
This could lead to an \"isolated section of conductor\" error if this \n \
node is referenced later.\n",listpointer->name);
      }
      /* add this name to the master list */
      append_pnlist(create_pn(listpointer->name, tempnode), indsys);
      if (is_node_in_list(tempnode, grndp->usernodes) == 1) {
	fprintf(stderr, "Warning: ground plane node %s coincident with another node\n", listpointer->name);
	if (is_nonuni_gp(grndp))
	  fprintf(stderr, "  at coordinates (%lg, %lg, %lg) (meters)\n",
		  tempnode->x, tempnode->y, tempnode->z);
      }
      else {
	if (is_nonuni_gp(grndp)) {
          if (indsys->opts->debug == ON)
            /* announce effective circumference of contact point */
            printf("Effective perimeter of contact point %s in meters: %lg\n",
                   listpointer->name, get_perimeter(tempnode->gp_node));
	}
	/* make a pseudo_seg between other ground plane usernodes */
	/* (to be used by graph searching algs) */
	grndp->fake_seg_list
	  = make_new_fake_segs(tempnode, grndp->usernodes, grndp->fake_seg_list);
	grndp->usernodes = add_node_to_list(tempnode, grndp->usernodes);
	tempnode->examined = 0;  /* examine these */
/*      sprintf(templine,".equiv %s %s \n",listpointer->name, tempnode->name);
        equivnodes(templine, indsys); */
      }
    }
  }

  /* removing the space allocated for the list_of_nodes */
  /* no, let's save it for regurgitation
  while(list_of_nodes != NULL){
    listpointer = list_of_nodes;
    list_of_nodes = list_of_nodes->next;
    free(listpointer->name);
    free(listpointer);
  }
  */
  grndp->usernode_coords = list_of_nodes;

  /* remove space allocated for list_of_holes */
  /* no, save it for regurgitation
  while(list_of_holes != NULL){
    holep = list_of_holes;
    list_of_holes = list_of_holes->next;
    free(holep->func);
    free(holep->vals);
    free(holep);
  }
  */
  grndp->list_of_holes = list_of_holes;


  plane_count++;
  grndp->next = NULL;

  /* check if another groundplane has the same name */
  onegp = indsys->planes;
  while(onegp != NULL && strcmp(onegp->name,grndp->name) != 0)
    onegp = onegp->next;

  if (onegp != NULL) {
    fprintf(stderr, "Error: Two ground planes with the name: %s\n",
	    grndp->name);
    err = 1;
  }

  /* adding plane to linked list of groundplanes in indsys */
  if(indsys->planes == NULL){
    indsys->planes = grndp;
    indsys->endplane = grndp;
  }
  else {
    indsys->endplane->next = grndp;
    indsys->endplane = grndp;
  }

  /* doing .external if the plane is a conductor */
  if(grndp->external == 1){
#if 1==0
    sprintf(templine, ".external %s ",grndp->innode->name);
    tempath = path_through_gp(grndp->innode, grndp->outnode, grndp);

    if(tempath == NULL){
      printf("tempath is NULL! \n");
      exit(1);
    }

    for(pathp = tempath; pathp != NULL; pathp = pathp->next){
      strcat(templine, " ");                 /* spacing between segments */
      strcat(templine, pathp->seg->name);
    }

    strcat(templine, " \n");             /* end of line marker */

    addexternal(templine, indsys);
#endif
  }
  /* end .external layout for the conductor */

  *retplane = grndp;
  return err;
  }
/*------------------------------------------------------------------------------*/
/*                          end of groundplane code                             */
/*------------------------------------------------------------------------------*/

nothing(line)
char *line;
{
  if (line[0] == '+')
    printf("Nothing to continue.\n%s\n",line);
  else
    printf("What is the following line? \n%s\n", line);
  return (1);
}

char *getaline(fp)
FILE *fp;
{
  static char *all_lines = NULL;
  static int length = 0;
  char *line;
  int newlength;
  char *getoneline(), *plusline();

  if (length == 0) {
    length = MAXLINE;
    all_lines = (char *)malloc(length*sizeof(char));
  }

  line = getoneline(fp);
  if (line == NULL) {
    fprintf(stderr, "Unexpected end of file\n");
    exit(1);
  }

  strcpy(all_lines, line);

  /* concatenate any lines beginning with a '+' to all_lines */
  while( (line = plusline(fp)) != NULL) {
    if ((newlength = strlen(all_lines) + strlen(line) + 1) > length) {
      if ( (all_lines = realloc(all_lines, MAX(newlength,length+MAXLINE)))
	  == NULL ) {
	fprintf(stderr,"couldn't get more space for a line. Needed %d chars\n",
		newlength);
	exit(1);
      }
      else
	length = MAX(newlength,length+MAXLINE);
    }
    strcat(all_lines, line);
  }

  return all_lines;
}



char *plusline(fp)
FILE *fp;
{
  char *tmpline, *getoneline();

  tmpline = getoneline(fp);

  while(tmpline != NULL && tmpline[0] == '*')
    tmpline = getoneline(fp);

  if (tmpline == NULL)
    return NULL;
  else if (tmpline[0] == '+') {
    tmpline++;
    return tmpline;
  }
  else {
    savealine(tmpline);
    return NULL;
  }
}

/* a variable just for the following functions */
static int keep = 0;

char *getoneline(fp)
FILE *fp;
{
  static char line[MAXLINE] = { '\0' };
  char *retchar;

  if (keep) {
    keep = 0;
    return line;
  }
  else
    do {
      retchar = fgets(line, MAXLINE, fp);
    } while(retchar != NULL && !notblankline(line));

  if (retchar != NULL && strlen(line) == MAXLINE - 1)
    fprintf(stderr,"Warning: line may be too long:\n%s\n",line);

  if (retchar == NULL)
    return NULL;
  else
    return line;
}

savealine(line)
char *line;
{
  if (keep != 0) {
    printf("already have one line stored\n");
    exit(1);
  }
  else
    keep = 1;
}

int notblankline(string)
char *string;
{
   while( *string!='\0' && isspace(*string))
     string++;

   if (*string == '\0') return 0;
     else return 1;
}

tolowercase(line)
char *line;
{
  while(*line != '\0') {
    *line = tolower(*line);
    line++;
  }
}

is_nonuni_gp(gp)
     GROUNDPLANE *gp;
{
  return (gp->nonuni != NULL);
}
