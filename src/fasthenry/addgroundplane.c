#include <stdio.h>
#include <math.h>
#include "induct.h"
// Enrico
#include <string.h>

#define PIOVER2 1.570796327

/*-------------------------MATH FUNCTIONS-------------------------------*/
/* calculates the square of x                                           */
double square(x)
double x;
{
  double val;

  val = x * x;
  return val;
}

/* function to find the length of a vector given the points that make up the */
/* vector.                                                                   */
double lengthof(p1, p2,x, y, z)
int p1, p2;
double x[4], y[4], z[4];
{
  double value;

  value = sqrt(((square((x[p2] - x[p1]))) +
		(square((y[p2] - y[p1]))) +
		(square((z[p2] - z[p1])))));

  return value;
}

/* function that calculates the length of a vector that connects two given */
/* points.                                                                 */
double lengthof2(a, x, y, z)
double a[3];
double x, y, z;
{
  double temp1, temp2, temp3, value;

  temp1 = a[0] - x;
  temp2 = a[1] - y;
  temp3 = a[2] - z;

  value = sqrt((square(temp1)) + (square(temp2)) + (square(temp3)));

  return value;
}

/* function that returns the innerproduct of two vectors given the points that */
/* make up the vectors.                                                        */
double innerproduct(p1,p2,p3,p4,x, y, z)
int p1,p2,p3,p4;
double x[4], y[4], z[4];
{
  double value;

  value = ((x[p2] - x[p1])*(x[p4] - x[p3])) +
          ((y[p2] - y[p1])*(y[p4] - y[p3])) +
	  ((z[p2] - z[p1])*(z[p4] - z[p3]));

  return value;
}

/* function that returns the angle between two vectors given the points that */
/* make up the vectors.                                                      */
double findangle(p1,p2,p3,p4,x, y, z)
int p1,p2,p3,p4;
double x[4], y[4], z[4];
{

  double temp1, temp2, temp3, value;

  temp1 = innerproduct(p1,p2,p3,p4,x, y, z);
  temp2 = lengthof(p1,p2,x, y, z);
  temp3 = lengthof(p3,p4,x, y, z);

  value = acos((temp1/(temp2 * temp3)));

  return value;
}

/*---------------------------SORTING ROUTINES (currently not used)----------*/
/* this is a primitive bubble sorting routing that takes is the order n*n... */
/* for the purpose of sorting a small array, it is effective.                */
void order(p, q)
int *p, *q;
{
  int temp;

  if(abs(*p) > abs(*q)){
    temp = *p;
    *p = *q;
    *q = temp;
  }
}

void bubblesort(array, n)
int *array, n;                 /* n is the size of the array of pointers */
{
  int i, j;

  for(i = 0; i < (n - 1); i++)
    for(j = (n - 1); i < j; --j)
      order(&array[j-1], &array[j]);
}

/*-------------------PLANE CALCULATION FUNCTIONS---------------------------*/
/* function to set the point numbering standard for a given ground plane.   */
/* x,y,z are matrices of the three input points, o1,mid,o2 are the returned */
/* point numbers.  (example: mid = 0 o1 = 1 o2 = 2 ... the middle is point 0*/
/* and the others are points 1 and 2.                                       */

int checkmiddlepoint(x, y, z, o1, o2, mid)
double x[4], y[4], z[4];
int o1, o2, mid;
{
  double temp1, temp2, temp3;
  int to1, to2, tmid;
  int value;

  temp1 = PIOVER2 - (findangle(o1,mid, o1, o2, x, y, z));
  temp2 = PIOVER2 - (findangle(o1, mid, mid, o2, x, y, z));
  temp3 = PIOVER2 - (findangle(o1, o2, o2, mid, x, y, z));

  if((fabs(temp1) < fabs(temp2)) && (fabs(temp1) < fabs(temp3))){
    tmid = 0; to1 = 1; to2 = 2;}
   else if((fabs(temp2) < fabs(temp3)) && (fabs(temp2) < fabs(temp1))){
    tmid = 1; to1 = 0; to2 = 2;}
  else {
    tmid = 2; to1 = 0; to2 = 1;}

  if((mid == tmid) && (o1 == to1) && (o2 == to2))
    return 1;
  else
    return 0;

}

/* function that checks if the given coordinates satisfy the equation of the */
/* given plane.                                                              */
int checkplaneformula(xcord, ycord, zcord,x,y,z,mid,o1,o2)
double xcord[4], ycord[4], zcord[4],x,y,z;
int mid,o1,o2;
{
  double xpart, ypart, zpart, value;
  double equation;

  xpart = ((ycord[o1] - ycord[mid])*(zcord[o2] - zcord[mid])) -
    ((zcord[o1] - zcord[mid])*(ycord[o2] - ycord[mid]));

  ypart = ((zcord[o1] - zcord[mid])*(xcord[o2] - xcord[mid])) -
    ((xcord[o1] - xcord[mid])*(zcord[o2] - zcord[mid]));

  zpart = ((xcord[o1] - xcord[mid])*(ycord[o2] - ycord[mid])) -
    ((ycord[o1] - ycord[mid])*(xcord[o2] - xcord[mid]));

  value = xpart*xcord[mid] + ypart*ycord[mid] + zpart*zcord[mid];

  equation = ( xpart * x) + (ypart * y) + (zpart * z);

  if (fabs(value - equation)/(xpart*xpart+ypart*ypart+zpart*zpart+1e-200)
      < EPS) {
    return 1;
  } else {
    printf("Possible internal error: something wrong with coordinates %lg %lg %lg\n",x,y,z);
    return 0;
  }
}

/* function that finds the width of a segment given its coordinates and the */
/* dimensions of a plane.                                                   */
double findsegmentwidth(x, y, z, mid, o1, o2, dim)
double x[4], y[4], z[4];
int mid, o1, o2, dim;
{
  double value, theta, temp1;

  theta = findangle(mid, o1, mid, o2, x, y, z);
  temp1 = lengthof(mid, o1, x, y, z);

  value = (temp1 * (sin(theta)))/dim;
  return value;
}

void doincrement(x, y, z, xi, yi, zi, dim, dx, dy, dz)
double x, y, z, xi, yi, zi;
int dim;
double *dx, *dy, *dz;
{
  *dx = (x - xi)/dim;
  *dy = (y - yi)/dim;
  *dz = (z - zi)/dim;
}

void dounitvector(x, y, z, xi, yi, zi, wx, wy, wz)
double x, y, z, xi, yi, zi;
double *wx, *wy, *wz;
{
  double temp1, w1, w2, w3;

  w1 = xi - x;
  w2 = yi - y;
  w3 = zi - z;
  temp1 = sqrt(w1*w1 + w2*w2 + w3*w3);

  *wx = w1/temp1;
  *wy = w2/temp1;
  *wz = w3/temp1;
}

/*---------------------------MISCELLANEOUS FUNCTIONS-------------------------*/
/* function that lays out and stores the grid matrix of a groundplane for use*/
/* in Matlab.(grids are used to visualize the current distribution in a plane*/
void fillgrids(plane)
GROUNDPLANE *plane;
{
  int rows, cols, i, j;
  int filnumber;
  SEGMENT *seg;

  cols = plane->seg2 + 1;
  rows = plane->seg1 + 1;
  plane->row[0] = rows;
  plane->col[0] = cols;

  /* allocate the space */
  plane->grid1 = (double **)MatrixAlloc(rows, cols, sizeof(double));

  for(i = 0; i < cols; i++){
    for(j = 1; j < rows; j++){
      seg = plane->segs1[j - 1][i];
      filnumber = seg->filaments[0].filnumber;

      plane->grid1[j][i] = (double)filnumber + 1.0;
    }
  }

  cols = plane->seg2 + 1;
  rows = plane->seg1 + 1;
  plane->row[1] = rows;
  plane->col[1] = cols;

  /* allocate the space for grid2 */
  plane->grid2 = (double **)MatrixAlloc(rows, cols, sizeof(double));

  for(i = 0; i < (cols - 1); i++){
    for(j = 0; j < rows; j++){
      seg = plane->segs2[j][i];
      filnumber = seg->filaments[0].filnumber;

      plane->grid2[j][i] = (double)filnumber + 1.0;
    }
  }
}

void make_nodelist(node, name, x, y, z)
NODELIST *node;
char *name;
double x, y, z;
{
  strcpy(node->name, name);
  node->x = x;
  node->y = y;
  node->z = z;

}

double find_coordinate(plane, x, y, z, flag)
GROUNDPLANE *plane;
double x, y, z;
int flag;
{
  double coordinate, xpart, ypart, zpart, value, equation, divpart;

  xpart = ((plane->y[0] - plane->y[1])*(plane->z[2] - plane->z[1])) -
    ((plane->z[0] - plane->z[1])*(plane->x[2] - plane->x[1]));
  ypart = ((plane->z[0] - plane->z[1])*(plane->x[2] - plane->x[1])) -
    ((plane->x[0] - plane->x[1])*(plane->z[2] - plane->z[1]));
  zpart = ((plane->x[0] - plane->x[1])*(plane->y[2] - plane->y[1])) -
    ((plane->y[0] - plane->y[1])*(plane->x[2] - plane->x[1]));

  value = xpart*plane->x[1] + ypart*plane->y[1] + zpart*plane->z[1];

  if(flag == 0){                           /* need to find the x_coordinate */
    equation = ypart*y + zpart*z;
    divpart = xpart;
  }
  else if(flag == 1){                      /* need to find the y_coordinate */
    equation = xpart*x + zpart*z;
    divpart = ypart;
  }
  else if(flag == 2){                      /* need to find the z_coordinate */
    equation = xpart*x + ypart*y;
    divpart = zpart;
  } else {
    printf("flag does not match coordinates \n");
    exit(1);
  }

  if(divpart == 0.0){
    return 0.0;
  } else {
    coordinate = (value - equation)/divpart;
    return coordinate;
  }

}

void findrefnodes(plane, begnode, endnode, b0, b1, e0, e1)
GROUNDPLANE *plane;
NODES *begnode, *endnode;
int *b0, *b1, *e0, *e1;
{
  int node1, node2, i, j;
  NODES *tempnode;

  node1 = plane->seg1 + 1;
  node2 = plane->seg2 + 1;

  for(i = 0; i < node1; i++){
    for(j = 0; j < node2; j++){
      tempnode = plane->pnodes[i][j];

      if(begnode == tempnode){
	*b0 = i;
	*b1 = j;
      }
      if(endnode == tempnode){
	*e0 = i;
	*e1 = j;
      }
    }
  }
}

SPATH *old_path_through_gp(nodein, nodeout, plane)
NODES *nodein, *nodeout;
GROUNDPLANE *plane;
{
  int deco1, deco2, counter1, counter2;
  int b0, b1, e0, e1, i, bc0, bc1, ec0,ec1;
  NODES *tempnode;
  SEGMENT *tempseg;
  SPATH *path;
  SPATH *tpath;
  SPATH *pathpointer;
  seg_ptr tempseg_ptr;

  NODES *find_nearest_node();

  tempseg_ptr.type = NORMAL;

  if((nodein->gp != plane) || (nodeout->gp != plane)){
    printf("nodes not in the groundplane !! \n");
    exit(1);
  }

  tempnode = find_nearest_gpnode(nodein->x, nodein->y, nodein->z, plane, &bc0, &bc1);
  tempnode = find_nearest_gpnode(nodeout->x, nodeout->y, nodeout->z, plane, &ec0, &ec1);

  findrefnodes(plane, nodein, nodeout, &b0, &b1, &e0, &e1);

  if (b0 == bc0 && b1 == bc1 && e0 == ec0 && e1 == ec1)
/*    printf("find_nearest_gpnode worked!\n"); */
    ;
  else
    printf("find_nearest_gpnode failed!!!\n");

  path = (SPATH *)MattAlloc(1, sizeof(SPATH));         /* allocate the space */
  tpath = path;
  tpath->next = NULL;

  deco1 = e0 - b0;
  deco2 = e1 - b1;
  counter1 = b0 + (deco1 > 0 ? 0.0 : -1.0);
  counter2 = b1 + (deco2 > 0 ? 0.0 : -1.0);

  for(i = 0; i < abs(deco1); i++){            /* segments from input to (o1) */

    pathpointer = tpath;                       /* pointer to previous path */

    tempseg = plane->segs1[counter1][b1];
    tempseg_ptr.segp = (void *)tempseg;
    tpath->seg = tempseg_ptr;                              /* assign the seg to the path */
    counter1 = counter1 + (deco1 > 0 ? 1.0 : -1.0);

    tpath->next = (SPATH *)MattAlloc(1, sizeof(SPATH));
    tpath = tpath->next;

  }

  if(deco1 > 0){
    if(counter1 != e0){
      printf("something wrong with .external stuff...\n");
      exit(1);
    }
  } else {
    if((counter1 + 1) != e0){
      printf("something wrong with .external stuff...\n");
      exit(1);
    }
  }

  for(i = 0; i < abs(deco2); i++){              /* segments from output to (o2) */

    pathpointer = tpath;

    tempseg = plane->segs2[e0][counter2];
    tempseg_ptr.segp = (void *)tempseg;
    tpath->seg = tempseg_ptr;                             /* assing the seg to the path */
    counter2 = counter2 + (deco2 > 0 ? 1.0 : -1.0);

    tpath->next = (SPATH *)MattAlloc(1,sizeof(SPATH));
    tpath = tpath->next;
  }

  tpath = pathpointer;
  tpath->next = NULL;

  return path;
}

NODES *findnode(node, plane, inout, var1, var2, samenode)
NODES *node, *samenode;
GROUNDPLANE *plane;
int inout[3];
double *var1, *var2;
{
  NODES *retnode;

  *var1 = lengthof2(inout, node->x, node->y, node->z);
  if(*var1 < *var2){
    *var2 = *var1;
    return node;
  }
  else {
    return samenode;
  }
}

#if 1==0
/*------------------------------ADDEXTERNAL FUNCTIONS------------------*/
/* makes the nodepath for a given segment in the function addextern()  */
void makenpath(seg, node)
SEGMENT *seg;
NODES *node;
{
  NPATH *npath;

  if(seg->conds == NULL){
    seg->conds = (NPATH *)MattAlloc(1, sizeof(NPATH));
    npath = seg->conds;
  }
  else{
    npath = seg->conds;
    while(npath->next != NULL)
      npath = npath->next;
    npath->next = (NPATH *)MattAlloc(1, sizeof(NPATH));
    npath = npath->next;
  }

  npath->node = node;
  npath->next = NULL;

}
#endif

/* makes the Mlist for the groundplane given a plane and parameters defining */
/* the current location of the overall Mlist.                                */
MELEMENT **old_makeMlist(plane, checksegs, M, mesh)
GROUNDPLANE *plane;
int *checksegs;
double **M;
int mesh;
{
  MELEMENT **pMlist, *melem, *melem2;
  FILAMENT **loopfils;
  SEGMENT *seg;
  NODES *node;
  int plusnode, tempmesh, number;
  int counter, i, j, k, n;
  int signofelem, *loop;

  counter = 0;
  tempmesh = mesh;

  pMlist = (MELEMENT **)MattAlloc(plane->numesh, sizeof(MELEMENT *));
  loop = (int *)MattAlloc(4, sizeof(int));
  loopfils = (FILAMENT **)MattAlloc(4, sizeof(FILAMENT *));

  for(i = 0; i < plane->seg2; i++){
    for(j = 0; j < plane->seg1; j++){

      node = plane->pnodes[j][i];
      plusnode = node->number;            /* initialize plusnode */

      for(k = 0; k < 4; k++){
	switch (k) {

	case 0:
	  seg = plane->segs1[j][i];
	  signofelem = -1.0;
	  break;
	case 1:
	  seg = plane->segs2[j + 1][i];
	  signofelem = -1.0;
	  break;
	case 2:
	  seg = plane->segs1[j][i + 1];
	  signofelem = 1.0;
	  break;
	case 3:
	  seg = plane->segs2[j][i];
	  signofelem = 1.0;
	  break;
	}

	number = seg->number;
	checksegs[number] = 1;             /* add to checksegs */
	                           /* MK 9/92: took out += 1 */

	for(n = 1; n < seg->num_fils; n++, tempmesh++){
	  melem = pMlist[tempmesh] = (MELEMENT *)MattAlloc(1, sizeof(MELEMENT));
	  melem->filindex = seg->filaments[n-1].filnumber;
	  melem->fil = &seg->filaments[n-1];
	  melem->sign = 1;
	  melem->mnext = (MELEMENT *)MattAlloc(1, sizeof(MELEMENT));

	  melem = melem->mnext;
	  melem->filindex = seg->filaments[n].filnumber;
	  melem->fil = &seg->filaments[n];
	  melem->sign = -1;
	  melem->mnext = NULL;
	}

	/* save the little mesh */
	loop[k] = signofelem * (seg->number + 1);
	loopfils[k] = &seg->filaments[0];
      }

      /* put out pMlist */
      for(k = 0; k < 4; k++){
/*	M[abs(loop[k]) - 1][tempmesh] += (loop[k] > 0 ? 1.0 : -1.0); */
	melem = (MELEMENT *)MattAlloc(1, sizeof(MELEMENT));
	melem->filindex = abs(loop[k]) - 1;
	melem->fil = loopfils[k];               /* assign filament pointer */
	melem->sign = (loop[k] > 0 ? 1.0 : -1.0);
	if(k == 0){
	  melem->mnext = NULL;
	  pMlist[counter] = melem;
	}
	else {
	  /* where to put melem */
	  melem2 = pMlist[counter];
	  if(melem2->filindex > melem->filindex){    /* put at the beginning */
	    melem->mnext = melem2;
	    pMlist[counter] = melem;
	  }
	  else{
	    while(melem2->mnext != NULL && melem2->mnext->filindex < melem->filindex)
	      melem2 = melem2->mnext;
	    /* insert in the middle */
	    melem->mnext = melem2->mnext;
	    melem2->mnext = melem;
	  }
	}
      }
      tempmesh++; counter++;                             /* mesh counter */
    }
  }

  if(counter != plane->numesh){
    printf("something wrong with meshes, numesh != counter \n");
    exit(1);
  }

  return pMlist;
}

SPATH *path_through_gp(nodein, nodeout, plane)
NODES *nodein, *nodeout;
GROUNDPLANE *plane;
{
  SPATH *path, *path_through_nonuni_gp();

#ifdef OLD_PATH
  return old_path_through_gp(nodein, nodeout, plane);
#endif

  if((nodein->gp != plane) || (nodeout->gp != plane)){
    printf("path_though_gp: nodes not in the groundplane !! \n");
    exit(1);
  }

  if (is_nonuni_gp(plane))
    path = path_through_nonuni_gp(nodein, nodeout, plane);
  else {
    /* we will use node->examined to mark if each node of a failed path */
    clear_plane_marks(plane);

    path = get_a_path(nodein, nodeout, plane, NULL, 0, 0);
  }

  if (path == NULL) {
    fprintf(stderr, "Error: Couldn't find a path from node %s to %s.\n",
	    nodein->name, nodeout->name);
    fprintf(stderr, "\
This may be due to \n\
  1) an isolated section of conductor in the ground plane formed\n\
     with multiple HOLE directives,\n\
  2) a .external statement which uses a node which is within a hole,\n\
  3) or a bug.\n");
    exit(1);
  }

  return path;
}

/* this returns a path from node to node_goal.  It takes one step and then
   calls itself recursively.  */
SPATH *get_a_path(node, node_goal, plane, nodes_so_far, s1_momentum,
		  s2_momentum)
NODES *node, *node_goal;
GROUNDPLANE *plane;
NPATH *nodes_so_far;
int s1_momentum, s2_momentum;
{
#define MAXchoices 4
  choice_list choices[MAXchoices];
  int nchoices, s1, s2, s1_goal, s2_goal, i;
  SEGMENT ***segs1 = plane->segs1;
  SEGMENT ***segs2 = plane->segs2;
  NODES ***pnodes = plane->pnodes;
  SPATH *path;
  seg_ptr sptr;
  int nodes1 = plane->num_nodes1;
  int nodes2 = plane->num_nodes2;

  nchoices = 0;
  s1 = node->s1;
  s2 = node->s2;
  s1_goal = node_goal->s1;
  s2_goal = node_goal->s2;

  sptr.type = NORMAL;

  /* gather all the choices of paths from this node and rank them */

  if (s1 > 0)
    nchoices += add_choice(&choices[nchoices], nodes_so_far, segs1[s1-1][s2],
			   pnodes[s1-1][s2],s1 - s1_goal, -1, 0, s1_momentum,
			   s2_momentum);

  if (s1 < nodes1 - 1)
    nchoices += add_choice(&choices[nchoices], nodes_so_far, segs1[s1][s2],
			   pnodes[s1+1][s2],s1_goal - s1, 1, 0, s1_momentum,
			   s2_momentum);

  if (s2 > 0)
    nchoices += add_choice(&choices[nchoices], nodes_so_far, segs2[s1][s2-1],
			   pnodes[s1][s2-1], s2 - s2_goal, 0, -1, s1_momentum,
			   s2_momentum);

  if (s2 < nodes2 - 1)
    nchoices += add_choice(&choices[nchoices], nodes_so_far, segs2[s1][s2],
			   pnodes[s1][s2+1],s2_goal - s2, 0, 1, s1_momentum,
			   s2_momentum);

  /* sort in increasing order by rank */
  sort_choices(choices, nchoices);

  /* add current node to list of passed through nodes */
  nodes_so_far = add_node_to_list(node, nodes_so_far);

  for(path = NULL, i = 0; path == NULL && i < nchoices; i++) {
    sptr.segp = (void *)choices[i].seg;
    if (node_goal == choices[i].node) {
      path = add_seg_to_list(sptr, path);
      increment_usage(choices[i].seg);
    }
    else {
      /* call myself with next node */
      path = get_a_path(choices[i].node, node_goal, plane, nodes_so_far,
			choices[i].s1_mom, choices[i].s2_mom);
      if (path != NULL) {
	/* there is a path, so lets add this seg to it */
	path = add_seg_to_list(sptr,path);
	increment_usage(choices[i].seg);
      }
    }
  }

  /* free this node */
  free(nodes_so_far);

  if (path == NULL)
    /* mark it is a node from which no path was found */
    node->examined++;

  return path;
}

/* sort choices by rank.  Not the most efficient */
sort_choices(choices, num)
choice_list *choices;
int num;
{
  int i, j;
  choice_list temp;

  for(i = 0; i < num - 1; i++)
    for(j = num - 1; j > i; j--)
      if (choices[j-1].rank > choices[j].rank) {
	temp = choices[j-1];
	choices[j-1] = choices[j];
	choices[j] = temp;
      }
}

clear_marks(indsys)
SYS *indsys;
{
  SEGMENT *seg;
  NODES *node;

  for(seg = indsys->segment; seg != NULL; seg = seg->next)
    seg->is_deleted = 0;

  for(node = indsys->nodes; node != NULL; node = node->next)
    node->examined = 0;
}

increment_usage(seg)
SEGMENT *seg;
{
  seg->is_deleted++;
}

dump_mesh_coords(indsys)
SYS *indsys;
{
  FILE *fp;
  int nmeshes, i, counter, j;
  double **line_list, *temp;
  int buffer_rows = 5000, buffer_cols = 4, nowrite;
  MELEMENT *melem;
  FILAMENT *fil;
  char outfname[80];

  concat4(outfname,"meshes",indsys->opts->suffix,".mat");
  if ( (fp = fopen(outfname,"w")) == NULL) {
    fprintf(stderr, "No open meshes.mat\n");
    exit(1);
  }

    j = 1000;
#ifdef DEC
    j = 0000;
#endif

  nmeshes = indsys->tree_meshes;

  temp = (double *)Gmalloc(buffer_rows*buffer_cols*sizeof(double));
  line_list = (double **)Make_C_array(temp, buffer_rows, buffer_cols,
				      sizeof(double));

  counter = 0;
  for(i = 0; i < nmeshes && counter < buffer_rows - 4; i++) {
    for(melem = indsys->Mlist[i]; melem != NULL && counter < buffer_rows - 4;
	melem = melem->mnext) {
      fil = melem->fil;
      line_list[counter][0] = fil->x[0];
      line_list[counter][1] = fil->y[0];
      line_list[counter][2] = fil->z[0];
      line_list[counter][3] = i + 1;
      line_list[counter+1][0] = fil->x[1];
      line_list[counter+1][1] = fil->y[1];
      line_list[counter+1][2] = fil->z[1];
      line_list[counter+1][3] = i + 1;
      counter += 2;
    }
    line_list[counter][3] = -1.0;
    counter += 2;
  }

  for(i = counter; i < buffer_rows; i++)
    line_list[i][3] = 0;

  savemat_mod(fp, j + 100, "meshes", buffer_rows, buffer_cols, 0,
	      line_list[0], (double *)NULL, 0,
	      buffer_rows*buffer_cols);

  fclose(fp);

}

dump_ascii_mesh_coords(indsys)
SYS *indsys;
{
  FILE *fp;
  int nmeshes, i, counter, j;
  double **line_list, *temp;
  int buffer_rows = 5000, buffer_cols = 4, nowrite;
  MELEMENT *melem;
  FILAMENT *fil;
  char outfname[80];

  concat4(outfname,"meshes",indsys->opts->suffix,".mat");
  if ( (fp = fopen(outfname,"w")) == NULL) {
    fprintf(stderr, "No open meshes.mat\n");
    exit(1);
  }

  nmeshes = indsys->tree_meshes;

  counter = 0;
  for(i = 0; i < nmeshes; i++) {
    for(melem = indsys->Mlist[i]; melem != NULL && counter < buffer_rows - 4;
	melem = melem->mnext) {
      fil = melem->fil;
      fprintf(fp, "%lg %lg %lg %d\n %lg %lg %lg %d\n",
	      fil->x[0], fil->y[0], fil->z[0], i + 1,
	      fil->x[1], fil->y[1], fil->z[1], i + 1);
      counter += 2;
    }
    fprintf(fp, "0.0 0.0 0.0 -1.0\n");
    fprintf(fp, "0.0 0.0 0.0 0.0\n");
    counter += 2;
  }

  fprintf(fp, "0.0 0.0 0.0 0.0\n");
  fclose(fp);

}

/* This makes a chunk of memory look like a doubly subscripted C array.
   It takes a continuous chunk of memory  of size rows*cols*size bytes,
   beginning at 'start',
   and returns an array of pointers to the beginning of each row.
*/
void **Make_C_array(start, rows, cols, size)
void *start;
int rows, cols, size;
{
  char **ptr;
  int i;

  ptr = (char **)Gmalloc(rows*sizeof(char *));

  if (sizeof(char) != 1) {
    fprintf(stderr,"oops");
    exit(1);
  }

  for(i = 0; i < rows; i++)
    ptr[i] = (char *)start + i*cols*size;

  return (void **)ptr;
}

/* fills the choice_list structure if this is a valid choice */
/* It returns 1 if it valid, 0 otherwise */
int add_choice(choice, nodes_so_far, seg, node, is_right_direction, new_s1_mom,
	       new_s2_mom, s1_momentum, s2_momentum)
choice_list *choice;
int is_right_direction;
NPATH *nodes_so_far;
SEGMENT *seg;
NODES *node;
int new_s1_mom, new_s2_mom;
int s1_momentum, s2_momentum;
{
  static int overlap = 0;  /* weight changed to 0 from 1 since new precond */
                           /* shouldn't care.   6/94  */
  static int wrong_way = 6;
  static int half_wrong_way = 2;
  static int against_momentum = 1;

  int with_momentum = new_s1_mom*s1_momentum + new_s2_mom*s2_momentum;

  if (seg != NULL && !is_orignode_in_list(node, nodes_so_far)
      && node->examined == 0) {
    choice->seg = seg;
    choice->node = node;
    choice->s1_mom = new_s1_mom;
    choice->s2_mom = new_s2_mom;

    choice->rank = 0;

    if (is_right_direction < 0) {
      /* This is the wrong way, avoid this direction if at all possible*/
      choice->rank += wrong_way;
      if (with_momentum != 1)
	/* this is against the momentum */
	choice->rank += against_momentum;
    }
    else {
      if (is_right_direction == 0)
	/* not completely the wrong way */
	choice->rank += half_wrong_way;

      choice->rank += (choice->seg->is_deleted)*overlap;
    }

    return 1;
  }
  else
    return 0;
}

clear_plane_marks(plane)
GROUNDPLANE *plane;
{
  int i,j;
  int nodes1 = plane->num_nodes1;
  int nodes2 = plane->num_nodes2;
  NODES ***pnodes = plane->pnodes;

  for(i = 0; i < nodes1; i++)
    for(j = 0; j < nodes2; j++)
      pnodes[i][j]->examined = 0;
}
