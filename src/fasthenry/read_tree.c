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
/* reads a tree from the standard input and fills the gcell structure */

#include <stdio.h>
#include "gp.h"
#include "induct.h"

#define GP_DEBUG FALSE

  
#define MAXLINE 1000

process_plane(grndp, fp, indsys)
     GROUNDPLANE *grndp;
     FILE *fp;
     SYS *indsys;
{
  Nonuni_gp *gp;
  FILE *fp2;
  ContactList *contactp;

  gp = (Nonuni_gp *)gp_malloc(sizeof(Nonuni_gp));

  /* have the two structures point to each other */
  grndp->nonuni = gp;
  gp->grndp = grndp;

  if (fp != NULL) {
    /* read tree from file */
    if (readTree(fp, gp) != 0)
      return 1;
  }
  else {
    gp->root_cell = new_Gcell(CLEAR);
    gp->num_cells = 1;
    gp->root_cell->index = 1;
    gp->root_cell->children = NULL;
    gp->root_cell->children_type = NONE;
    gp->root_cell->parent = NULL;
  }

  /*
  gp->thickness = grndp->thick;
  gp->n_hinc = grndp->hinc;
  */

  /* compute relative coordinate system for gp */
  set_gp_coord_system(grndp, gp);

  /* print leaf cells so we can see what this looks like */
  /* dump_leaf_cells_to_file(gp->root_cell, "hier.qui");*/

  gp->num_z_pts = 1;
  gp->z_pts = (double *)gp_malloc(1*sizeof(double));
  gp->z_pts[0] = 0.5;  /* doesn't matter */

  /* determine nodes and adjacencies */
  process_tree(gp);

  /* scan contacts for and initial grid to insure it's done first */
  for(contactp = gp->grndp->list_of_contacts; contactp != NULL; 
      contactp = contactp->next) {
    /* there should be only one, so if more than one gets called, we'll
       get an error */
    if (strncmp("initial",contactp->func,7) == 0) {
      make_contacts(contactp, gp);
      contactp->done = TRUE;
    }
  }

  for(contactp = gp->grndp->list_of_contacts; contactp != NULL; 
      contactp = contactp->next)
    /* do all the non-initial_grid  and non equiv */
    if (contactp->done == FALSE && strncmp("equiv_rect",contactp->func,10) != 0)
      {
        make_contacts(contactp, gp);
        contactp->done = TRUE;
      }

  for(contactp = gp->grndp->list_of_contacts; contactp != NULL; 
      contactp = contactp->next)
    /* do all the non-initial_grid  and non equiv */
    if (contactp->done == FALSE && strncmp("equiv_rect",contactp->func,10) == 0)
      {
        make_contacts(contactp, gp);
        contactp->done = TRUE;
      }

  if (indsys->opts->makeFastCapFile & HIERARCHY)
    /* print leaf cells after contacts */
    dump_leaf_cells_to_file(gp->root_cell, "hier.qui");

  /* from node adjacency, generate fasthenry segments */
  generate_segs(gp, indsys);

  gp->grndp->numesh = gp->num_leaves + (gp->num_z_pts - 1)*gp->num_seg_groups;

  return 0;

}

/* the nonuni gp wants an origin in one corner and x and y axes on
   the edges and z in the thickness.  Let's set up a relative coord
   to global */

set_gp_coord_system(grndp, gp)
     GROUNDPLANE *grndp;
     Nonuni_gp *gp;
{
  double dx,dy,dz;
  double magx, magy;
  double *x = grndp->x;
  double *y = grndp->y;
  double *z = grndp->z;
  double thickness = gp->grndp->thick;

  /* compute x-axis unit vector */
  
  dx = x[1] - x[0];
  dy = y[1] - y[0];
  dz = z[1] - z[0];
  magx = sqrt(dx*dx + dy*dy + dz*dz);

  gp->ux_x = dx/magx;
  gp->ux_y = dy/magx;
  gp->ux_z = dz/magx;

  dx = x[2] - x[1];
  dy = y[2] - y[1];
  dz = z[2] - z[1];
  magy = sqrt(dx*dx + dy*dy + dz*dz);

  gp->uy_x = dx/magy;
  gp->uy_y = dy/magy;
  gp->uy_z = dz/magy;

  if (magx != grndp->length1 || magy != grndp->length2)
    GP_PANIC("How can nonuni magx and length1 be different?");
  
  /* set the root cell coordinates */
  set_cell_coords(gp->root_cell, 0.0, 0.0, magx, magy);

  /* do z direction so we can figure out where the origin is */
  gp->uz_x = gp->ux_y * gp->uy_z - gp->ux_z * gp->uy_y;
  gp->uz_y = gp->ux_z * gp->uy_x - gp->ux_x * gp->uy_z;
  gp->uz_z = gp->ux_x * gp->uy_y - gp->ux_y * gp->uy_x;

  /* origin is (x[0],y[0],z[0]) - 0.5*thickness * uz; */
  gp->x0 = x[0] - 0.5*thickness* gp->uz_x;
  gp->y0 = y[0] - 0.5*thickness* gp->uz_y;
  gp->z0 = z[0] - 0.5*thickness* gp->uz_z;
}

/* convert global xyz to plane xyz */
get_nonuni_coords(x, y, z, gp, xr, yr, zr)
     double x,y,z,*xr,*yr,*zr;
     Nonuni_gp *gp;
{
  double rx,ry,rz;

  /*local(xr,yr,zr) = ((xyz - origin)*ux, (xyz - origin)*uy, (xyz-origin)*uz)*/

  /* xyz - origin */
  rx = x - gp->x0;
  ry = y - gp->y0;
  rz = z - gp->z0;

  *xr = rx*gp->ux_x + ry*gp->ux_y + rz*gp->ux_z;
  *yr = rx*gp->uy_x + ry*gp->uy_y + rz*gp->uy_z;
  *zr = rx*gp->uz_x + ry*gp->uz_y + rz*gp->uz_z;
}

/* convert plane xyz point to global xyz point */
get_global_coords(x, y, z, gp, xg, yg, zg)
     double x,y,z, *xg, *yg, *zg;
     Nonuni_gp *gp;
{
  double xv, yv, zv;

  /* get global vector direction */
  get_global_vec(x, y, z, gp, &xv, &yv, &zv);
  
  /* shift global vector from origin */
  *xg = gp->x0 + xv;
  *yg = gp->y0 + yv;
  *zg = gp->z0 + zv;
}

/* convert plane xyz vector to global xyz vector */
get_global_vec(x, y, z, gp, xg, yg, zg)
     double x,y,z, *xg, *yg, *zg;
     Nonuni_gp *gp;
{
  *xg = x * gp->ux_x + y * gp->uy_x + z * gp->uz_x;
  *yg = x * gp->ux_y + y * gp->uy_y + z * gp->uz_y;
  *zg = x * gp->ux_z + y * gp->uy_z + z * gp->uz_z;
}

int readTree(fp, gp)
FILE *fp;
Nonuni_gp *gp;
{
  static char line[MAXLINE];
  char *retchar;
  int num_cells, i;
  Gcell *cells, *cell;
  Gcell *root_cell;
  double x,y,thick;
  Bi *bi_kids;
  int dum, nchild1, nchild2, numread;
  char typename[100], bintype[100];

  retchar = fgets(line, MAXLINE, fp);

  if (retchar == NULL) {
    fprintf(stderr, "Couldn't read first line\n");
    return 1;
  }
 
/*
  if (sscanf(line, "%d %lg %lg %lg", &num_cells, &x, &y, &thick) != 4) {
    fprintf(stderr, "First line should be an integer with total number of cells xlength ylength thickness\n");
    return 1;
  }
*/

  if (sscanf(line, "%d", &num_cells) != 1) {
    fprintf(stderr, "First line should be an integer with total number of cells\n");
    return 1;
  }

  gp->num_cells = num_cells;
  cells = new_Gcells(num_cells, CLEAR);

  if (cells == NULL) {
    fprintf(stderr, "Couldn't get space for %d cells\n",num_cells);
    return 1;
  }

  for(i = 0; i < num_cells; i++) {
    retchar = fgets(line, MAXLINE, fp);

    if (retchar == NULL) {
      fprintf(stderr, "Unexpected end of file after %d lines\n", i+1);
      return 1;
    }

    /* assume binary tree */
    numread = sscanf(line,"%d %s %s %d %d", &dum, typename, bintype, &nchild1, 
		     &nchild2);

    if (numread != 5 && !(numread == 2 && strcmp(typename,"NONE") == 0)) {
      fprintf(stderr, "At line %d: Couldn't read binary tree data\n",i+2);
      return 1;
    }

    cells[dum - 1].index = dum;

    if (strcmp(typename,"NONE") == 0) {
      cells[dum - 1].children_type = NONE;
      cells[dum - 1].children = NULL;
    }
    else {
      if (strcmp(typename,"B") != 0) {
	fprintf(stderr, "Err at line %d: only binary trees allowed\n",i+2);
	return 1;
      }
      
      if (dum > num_cells || nchild1 > num_cells || nchild2 > num_cells) {
	fprintf(stderr, "Err at line %d: index out of range\n",i+2);
	return 1;
      }
      
      if (cells[dum - 1].children != NULL) {
	fprintf(stderr, "uh oh, this cell has been defined before: %d\n",dum);
      }
      
      if (cells[nchild1 - 1].parent != NULL)
	fprintf(stderr, "uh oh, child %d for cell %d already has a parent\n",
		nchild1, dum);
      if (cells[nchild2 - 1].parent != NULL)
	fprintf(stderr, "uh oh, child %d for cell %d already has a parent\n",
		nchild2, dum);
      
      cell = &cells[dum - 1];
      
      /* assume binary tree */
      cell->children_type = BI;
      
      cell->children = (void *)gp_malloc( sizeof(Bi) );
      bi_kids = (Bi *)cell->children;
      
      if (strcmp(bintype,"NS") == 0)
	bi_kids->type = NS;
      else if (strcmp(bintype,"EW") == 0)
	bi_kids->type = EW;
      else {
	fprintf(stderr, "Unknown bintype at line %d\n", i+2);
	return 1;
      }
      
      bi_kids->child1 = &cells[nchild1 - 1];
      bi_kids->child2 = &cells[nchild2 - 1];
      
      cells[nchild1 - 1].parent = &cells[dum - 1];
      cells[nchild2 - 1].parent = &cells[dum - 1];
      
    }
  }
  

  if (cells[0].parent != NULL) {
    fprintf(stderr, "uh oh, cells[0] has a parent but should be root cell!\n");
    return 1;
  }

  gp->root_cell = root_cell = &cells[0];

/*
  gp->x0 = gp->y0 = gp->z0 = 0;
  gp->ux_x = 1; gp->ux_y = 0; gp->ux_z = 0;
  gp->uy_x = 0; gp->uy_y = 1; gp->uy_z = 0;
  gp->uz_x = 0; gp->uz_y = 0; gp->uz_z = 1;
  
  root_cell->x0 = root_cell->y0 = 0;
  root_cell->x1 = x;
  root_cell->y1 = y;
*/
  
  return 0;
}

/* set coordinates of this cell and its children */
set_cell_coords(cell,x0,y0,x1,y1)
     Gcell *cell;
     double x0, y0, x1, y1;
{
  cell->x0 = x0;
  cell->y0 = y0;
  cell->x1 = x1;
  cell->y1 = y1;

  switch(get_children_type(cell)) {
  case NONE:
    break;
  case BI:
    set_bi_coords( (Bi *)cell->children, x0, y0, x1, y1 );
    break;
  case GRID_2D:
    set_grid_coords( (Grid_2d *)cell->children, x0, y0, x1, y1);
    break;
  default:
    GP_PANIC("Unknown child type in set_cell_coords")
    break;
  }

}

set_bi_coords( two_kids, x0, y0, x1, y1)
  Bi *two_kids;
  double x0,y0,x1,y1;
{
  if (two_kids->type == NS) {
    /* north child */
    set_cell_coords(two_kids->child1, x0, (y0+y1)/2.0, x1, y1);
    /* south */
    set_cell_coords(two_kids->child2, x0, y0, x1, (y0+y1)/2.0);
  }
  else if (two_kids->type == EW) {
    /* east child */
    set_cell_coords(two_kids->child1, (x0+x1)/2.0, y0, x1, y1);
    /* west */
    set_cell_coords(two_kids->child2, x0, y0, (x0+x1)/2.0, y1);
  }
  else 
    GP_PANIC("Unknown bi child type in set_bi_coords");
}

set_grid_coords(grid, x0, y0, x1, y1)
  Grid_2d *grid;
  double x0,y0,x1,y1;
{
  int i,j;
  double xsize, ysize;
  double x0n, y1n;

  xsize = (x1 - x0)/grid->x_cells;
  ysize = (y1 - y0)/grid->y_cells;

  /* [0][0] is the top left NW (not SW) */
  for(i = 0; i < grid->y_cells; i++)
    for(j = 0; j < grid->x_cells; j++) {
      x0n = x0 + j*xsize;
      y1n = y1 - i*ysize;
      set_cell_coords(grid->kids[i][j], x0n, y1n - ysize, x0n + xsize, y1n);
    }
}
    

process_tree(gp)
     Nonuni_gp *gp;
{
  Info info;

  info.gp = gp;
  
  gp->num_leaves = 0;
  gp->num_nodes = 0;
  resolve_nodes(gp->root_cell, &info);

  /* we've set the edge info correctly */
  gp->is_edge_corrupted = FALSE;

#if GP_DEBUG == TRUE
  fprintf(stdout, "\nBefore deletion:\n");
  fprint_node_list(gp->nodelist, stdout);
#endif

  /* most are already deleted, but if they were the head of the list, they
     were left in */
  delete_dead_nodes(gp);

#if GP_DEBUG == TRUE
  fprintf(stdout, "\nAfter deletion:\n");
  fprint_node_list(gp->nodelist, stdout);
#endif

  determine_adjaceny(gp->nodelist);

#if GP_DEBUG == TRUE
  fprintf(stdout, "\nAfter adjacency:\n");
  fprint_node_list(gp->nodelist, stdout);
#endif

  compute_z_fils(gp);

}


/* recursive function to process tree */
resolve_nodes(cell, info)
     Gcell *cell;
     Info *info;
{
  switch (get_children_type(cell)) {
  case NONE:
    /* we are a leaf */
    if (!is_leaf(cell)) GP_PANIC("Not a leaf!");
    make_nodes(cell, info);
    info->gp->num_leaves++;    /* update number of leaves in tree */
    break;
  case BI:
    /* recursively call for BI children */
    resolve_bi_children(cell, info);
    break;
  case GRID_2D:
    GP_PANIC("You haven't implemented this yet!");
    /*resolve_grid_2d_children(cell, info); */
    break;
  default:
    GP_PANIC("Unknown child type in resolve_nodes");
    break;
  }
}

/* make "nodes" for the four corners of a leaf */
/* the information we want to collect is the cells attached to the node (4 max)
   and the adjacent nodes (4 max).  We start at the bottom with the leaves */
make_nodes(cell, info)
     Gcell *cell;
     Info *info;
{
  G_nodes *node;
  Nonuni_gp *gp = info->gp;

  /* NE */
  node = make_new_node(cell->x1, cell->y1, SW, cell, ++gp->num_nodes);
  cell->bndry.nodes[NE] = node;
  gp->nodelist = add_to_gnodelist(node, gp->nodelist);

  node = make_new_node(cell->x1, cell->y0, NW, cell, ++gp->num_nodes);
  cell->bndry.nodes[SE] = node;
  gp->nodelist = add_to_gnodelist(node, gp->nodelist);

  node = make_new_node(cell->x0, cell->y0, NE, cell, ++gp->num_nodes);
  cell->bndry.nodes[SW] = node;
  gp->nodelist = add_to_gnodelist(node, gp->nodelist);

  node = make_new_node(cell->x0, cell->y1, SE, cell, ++gp->num_nodes);
  cell->bndry.nodes[NW] = node;
  gp->nodelist = add_to_gnodelist(node, gp->nodelist);

#if GP_DEBUG == TRUE
  fprintf(stdout, "\nDone making nodes for: ");
  print_cell_and_kids(cell);
#endif

  /* set adjacency.  We will need to know only who is in either direction */
  /* ack, can't do this now since combining info will lose some of this.
     Instead, do it after tree is created */
/* removed and done later
  cell->nodes[NE]->adjacent[W] = cell->nodes[NW];
  cell->nodes[NE]->adjacent[S] = cell->nodes[SE];

  cell->nodes[NW]->adjacent[E] = cell->nodes[NE];
  cell->nodes[NW]->adjacent[S] = cell->nodes[SW];

  cell->nodes[SW]->adjacent[N] = cell->nodes[NW];
  cell->nodes[SW]->adjacent[E] = cell->nodes[SE];

  cell->nodes[SE]->adjacent[N] = cell->nodes[NE];
  cell->nodes[SE]->adjacent[W] = cell->nodes[SW];
*/
}

G_nodes *add_to_gnodelist(node, nodelist)
     G_nodes *node, *nodelist;
{
  if (nodelist == NULL) {
    node->prev = node->next = NULL;
  }
  else {
    node->prev = NULL;
    node->next = nodelist;
    nodelist->prev = node;
  }

  return node;
}

resolve_bi_children(cell, info)
     Gcell *cell;
     Info *info;
{
  Bi *two_kids;
  G_edges **edges;

  two_kids = (Bi *)cell->children;

  resolve_nodes(two_kids->child1, info);
  resolve_nodes(two_kids->child2, info);

#if GP_DEBUG == TRUE
  printf("\nBefore combining binary children edges: ");
  print_cell_and_kids(cell);
#endif

  /* definitely not a leaf... */
  edges = cell->bndry.edges;
  
  /* combine edge information of children for parent */
  /* an edge with pieces from both children will always have the north or east
     child first, and south or west second */
  if (two_kids->type == NS) {
    edges[NORTH] = make_one_edge(cell, two_kids->child1, NORTH);
    edges[SOUTH] = make_one_edge(cell, two_kids->child2, SOUTH);
    edges[EAST] = make_two_edge(cell, two_kids->child1, two_kids->child2,
			       EAST); /* make sure children in this order */
    edges[WEST] = make_two_edge(cell, two_kids->child1, two_kids->child2,
				       WEST);

    Combine_edges(two_kids->child1, SOUTH, two_kids->child2, NORTH);
  }
  else if (two_kids->type == EW) {
    edges[EAST] = make_one_edge(cell, two_kids->child1, EAST);
    edges[WEST] = make_one_edge(cell, two_kids->child2, WEST);
    edges[NORTH] = make_two_edge(cell, two_kids->child1, 
				 two_kids->child2, NORTH);
    edges[SOUTH] = make_two_edge(cell, two_kids->child1, 
					two_kids->child2, SOUTH);
 
    Combine_edges(two_kids->child1, WEST, two_kids->child2, EAST);
  }
  else
    GP_PANIC("UNKNOWN child type in resolve_bi_children");

#if GP_DEBUG == TRUE
  printf("\nAfter combining binary children edges: ");
  print_cell_and_kids(cell);
#endif

}

G_edges *make_one_edge(owner_cell, child_cell, dir)
     Gcell *owner_cell, *child_cell;
     char dir;
{
  G_edges *edge;
  int i;

  edge = (G_edges *)gp_malloc(sizeof(G_edges));

  edge->cells[0] = owner_cell;
  edge->dirs[0] = dir;
  edge->num_children = 1;
  edge->children = (Gcell **)gp_malloc(1*sizeof(Gcell *));
  edge->children[0] = child_cell;

  for(i = 1; i < NUM_E_CELLS; i++)
    edge->cells[i] = NULL;

  return edge;
}

G_edges *make_two_edge(owner_cell, child1, child2, dir)
     Gcell *owner_cell, *child1, *child2;
     char dir;
{
  G_edges *edge;
  int i;

  edge = (G_edges *)gp_malloc(sizeof(G_edges));

  edge->cells[0] = owner_cell;
  edge->dirs[0] = dir;
  edge->num_children = 2;
  edge->children = (Gcell **)gp_malloc(edge->num_children*sizeof(Gcell *));
  edge->children[0] = child1;
  edge->children[1] = child2;

  for(i = 2; i < NUM_E_CELLS; i++)
    edge->cells[i] = NULL;

  return edge;
}


Gcell *new_Gcells(num, flag)
     int num;
     int flag;
{
  Gcell *new;
  int i;
  
  new = (Gcell *)gp_malloc(num*sizeof(Gcell));

  /* gp_malloc uses calloc, so don't do this */
  /*
  if (flag == CLEAR) 
    for(i = 0; i < num; i++)
      init_Gcell(new[i]);
  */

  return new;

}

Gcell *new_Gcell(flag)
     int flag;
{
  Gcell *new;

  new = (Gcell *)gp_malloc(sizeof(Gcell));
  
  /* gp_malloc uses calloc, so this is unnecessary */
  /* 
  if (flag == CLEAR)
    init_Gcell(new);
  */

  return new;
}

/* this is never called */
init_Gcell(cell)
     Gcell *cell;
{
  int i;

  cell->children = cell->parent = NULL;
  cell->children_type = NONE;

  for(i = 0; i < NUMNODES; i++) /* but if edges is bigger, this won't work...*/
    cell->bndry.nodes[i] = NULL;

  cell->ishole = FALSE;
}

G_nodes *new_Gnode(flag)
     int flag;
{
  int i;
  G_nodes *new;
    
  new = (G_nodes *)gp_malloc(sizeof(G_nodes));

  if (flag == CLEAR) {
    for(i = 0; i < NUMADJ; i++) 
      new->adjacent[i] = NULL;
    for(i = 0; i < NUM_N_CELLS; i++)
      new->cells[i] = 0;
  }

  return new;

}
  
G_nodes *make_new_node(x, y, cell_dir, cell, index)
     double x,y;
     int cell_dir;
     Gcell *cell;
     int index;
{
  G_nodes *node;

  node = new_Gnode(CLEAR);
  node->cells[cell_dir] = cell;  /* i am the cell in the opposite direction */
  node->x = x;
  node->y = y;
  node->flag = ALIVE;
  node->n_segs = NULL;
  node->e_segs = NULL;
  node->prev = node->next = NULL;
  node->index = index;

  return node;
}

void *gp_malloc(size)
     int size;
{
  void *ptr;

  ptr = (void *)calloc(1, size);

  if (ptr == NULL)
    GP_PANIC("Couldn't get space!");
    
  return ptr;
}

/* this is the function that does the work in combining the information for 
   the nodes.  It takes two edges and combines them to one.

   The edges are represented by n-ary (or more) trees where an edge with 
   n children is divided into n equal sections and the n children are those
   edges.

   right now, edges can only be combined if, at a given depth in the tree,
   both edges have either:
   1. the same number of children
   2. one is a leaf and the other has multiple children.
*/

Combine_edges(cell1, dir1, cell2, dir2)
     Gcell *cell1, *cell2;
     char dir1, dir2;
{
  G_edges *edge1, *edge2;
  int i;

  while(!is_leaf(cell1) && cell1->bndry.edges[dir1]->num_children == 1)
    cell1 = cell1->bndry.edges[dir1]->children[0];

  while(!is_leaf(cell2) && cell2->bndry.edges[dir2]->num_children == 1)
    cell2 = cell2->bndry.edges[dir2]->children[0];

  if (is_leaf(cell1) || is_leaf(cell2))
    /* one is a leaf so time to propogate info to the other */
    combine_node_info(cell1, dir1, cell2, dir2);
  else {
    edge1 = cell1->bndry.edges[dir1];
    edge2 = cell2->bndry.edges[dir2];
    if (edge1->num_children != edge2->num_children) {
      fprintf(stderr, "Bad discretization: number of children doesn't match %d != %d\n", 
	      edge1->num_children, edge2->num_children);
      exit(1);
    }
    else {
      /* recursively call for children */
      for(i = 0; i < edge1->num_children; i++)
	Combine_edges(edge1->children[i], dir1, edge2->children[i], dir2);
    }
  }

}

combine_node_info(cell1, dir1, cell2, dir2)
     Gcell *cell1, *cell2;
     char dir1, dir2;
{
  char isleaf1, isleaf2;
  Gcell *leafcell, *nonleafcell;
  G_nodes *fareast_north, *farwest_south;
  char leafdir, nonleafdir;

  isleaf1 = is_leaf(cell1);
  isleaf2 = is_leaf(cell2);

  if (isleaf1) {
    leafcell = cell1;
    leafdir = dir1;
    nonleafcell = cell2;
    nonleafdir = dir2;
  }
  else {
    leafcell = cell2;
    leafdir = dir2;
    nonleafcell = cell1;
    nonleafdir = dir1;
  }

  /* put leaf cell node info into all nonleaf children nodes */
  /* and give back pointer to nodes on the end of edge which 
     will be filled by info in leaf cell nodes */
  give_cell_adjaceny(leafcell, leafdir, nonleafcell, nonleafdir, 
		     &fareast_north, &farwest_south);

  /* combine leafcell's two nodes with the far ones */
  combine_nodes(leafcell, leafdir, fareast_north, farwest_south);

}

/* put leafcell pointer into nodes on the nonleafcell edge */
/* and return the nodes on the end of the nonleafcell edge */
give_cell_adjaceny(leaf, leafdir, nonleaf, nonleafdir, 
		   pfareast_north, pfarwest_south)
     Gcell *leaf, *nonleaf;
     G_nodes **pfareast_north, **pfarwest_south;
     char leafdir, nonleafdir;
{
  int i;
  int num_kids;
  G_nodes *fareast_north, *farwest_south;
  G_edges *edge;
  
  if (is_leaf(nonleaf)) {
    if (nonleafdir == NORTH) {   
      nonleaf->bndry.nodes[NE]->cells[NW] = leaf;
      nonleaf->bndry.nodes[NW]->cells[NE] = leaf;
      *pfareast_north = nonleaf->bndry.nodes[NE];  /* far east */
      *pfarwest_south = nonleaf->bndry.nodes[NW];  /* far west */
    }
    if (nonleafdir == SOUTH) {
      nonleaf->bndry.nodes[SE]->cells[SW] = leaf;
      nonleaf->bndry.nodes[SW]->cells[SE] = leaf;
      *pfareast_north = nonleaf->bndry.nodes[SE];  /* far east */
      *pfarwest_south = nonleaf->bndry.nodes[SW];  /* far west */
    }
    if (nonleafdir == EAST) {
      nonleaf->bndry.nodes[NE]->cells[SE] = leaf;
      nonleaf->bndry.nodes[SE]->cells[NE] = leaf;
      *pfareast_north = nonleaf->bndry.nodes[NE];  /* far north */
      *pfarwest_south = nonleaf->bndry.nodes[SE];  /* far south */
    }
    if (nonleafdir == WEST) {
      nonleaf->bndry.nodes[NW]->cells[SW] = leaf;
      nonleaf->bndry.nodes[SW]->cells[NW] = leaf;
      *pfareast_north = nonleaf->bndry.nodes[NW];  /* far north */
      *pfarwest_south = nonleaf->bndry.nodes[SW];  /* far south */
    }
  }
  else {
    /* call recursively for the children.  Assumed in north/east 
       to south/west order */
    edge = nonleaf->bndry.edges[nonleafdir];
    for(i = 0; i < edge->num_children; i++) {
      give_cell_adjaceny(leaf, leafdir, edge->children[i], nonleafdir, 
			 &fareast_north, &farwest_south);
      if (i == 0)
	*pfareast_north = fareast_north;

      if (i == edge->num_children-1)  /* note: could be also i == 0 */
	*pfarwest_south = farwest_south;
    }
  }
}

/* combine the nodes on leafcell's leafdir with the given nodes */	
combine_nodes(leafcell, leafdir, fareast_north, farwest_south)
     Gcell *leafcell;
     char leafdir;
     G_nodes *fareast_north, *farwest_south;
{

  if (leafdir == NORTH) {
    replace_node(leafcell->bndry.nodes[NW], farwest_south);
    replace_node(leafcell->bndry.nodes[NE], fareast_north);
  }
  else if (leafdir == SOUTH) {
    replace_node(leafcell->bndry.nodes[SW], farwest_south);
    replace_node(leafcell->bndry.nodes[SE], fareast_north);
  }
  else if (leafdir == EAST) {
    replace_node(leafcell->bndry.nodes[SE], farwest_south);
    replace_node(leafcell->bndry.nodes[NE], fareast_north);
  }
  else if (leafdir == WEST) {
    replace_node(leafcell->bndry.nodes[SW], farwest_south);
    replace_node(leafcell->bndry.nodes[NW], fareast_north);
  }
  else
    GP_PANIC("combine_nodes: Unknown direction");
}

/* replace the cell's node in the node_dir direction with new node
   */
replace_node(old_node, new_node)
     G_nodes *new_node, *old_node;
{
  int i;

  if (old_node == new_node)
    return;     /* these were already combined by my edge sibling probably */

  for(i = 0; i < NUM_N_CELLS; i++) {
    if (old_node->cells[i] != NULL) {
      if (new_node->cells[i] != NULL) {
	if (new_node->cells[i] != old_node->cells[i])
	  GP_PANIC("In replace_node: Can't combine conflicting cell info!");
	/* else do nothing */
      }
      else {
	/* give the new node everything we care about for the old node */
	new_node->cells[i] = old_node->cells[i];
      }

      /* have the other cells point to this new node instead */
      old_node->cells[i]->bndry.nodes[opposite_dir(i)] = new_node;
    }
  }

#if GP_DEBUG==TRUE
  if (old_node->x != new_node->x || old_node->y != new_node->y)
    printf("replace_nodes: nodes don't have the same coords!\n");
#endif
      
  /* mark the old node as dead */
  kill_node(old_node);
}

kill_node(node)
     G_nodes *node;
{
#if GP_DEBUG == TRUE
  printf("killing node:");
  dump_node(node, stdout);
#endif

  node->flag = DEAD;
  remove_and_free(node);
}

/* delete first node in nodelist if it is marked as dead. NOT USED */
delete_first_node(gp)
     Nonuni_gp *gp;
{
  G_nodes *node;

  if (gp->nodelist != NULL && gp->nodelist->flag == DEAD) {
    node = gp->nodelist;
    gp->nodelist = node->next;
    node->next->prev = NULL;
    free_g_node(node);
  }
}
  
/* remove nodes from linked list marked "DEAD".  */
delete_dead_nodes(gp)
     Nonuni_gp *gp;
{
  G_nodes *node, *free_me;

  node = gp->nodelist;
  while(node != NULL) {
    if (node->flag == DEAD) {
      free_me = node;
      if (node->prev != NULL)
	node->prev->next = node->next;
      else
	gp->nodelist = node->next;

      if (node->next != NULL)
	node->next->prev = node->prev;
      
      node = node->next;
      free_g_node(free_me);
    }
    else
      node = node->next;
  }
}

remove_and_free(node)
     G_nodes *node;
{

  if (node->prev != NULL)
    node->prev->next = node->next;
  else
    return;  /*This is the head node.  Must take care of it later */
  
  if (node->next != NULL)
    node->next->prev = node->prev;
  
  free_g_node(node);
}

free_g_node(node)
     G_nodes *node;
{
  free(node);
}

determine_adjaceny(nodelist)
     G_nodes *nodelist;
{
  G_nodes *node;
  int i;

  for(node = nodelist; node != NULL; node = node->next) {
    /* do north */
    if (node->cells[NE] != node->cells[NW]) 
      node->adjacent[N] = get_adjacent_node(node, node->cells[NE], SW, WEST,
					    node->cells[NW], SE, EAST);

    /* do east */
    if (node->cells[NE] != node->cells[SE])
      node->adjacent[E] = get_adjacent_node(node, node->cells[NE], SW, SOUTH,
					     node->cells[SE], NW, NORTH);

    /* do south */
    if (node->cells[SE] != node->cells[SW])
      node->adjacent[S] = get_adjacent_node(node, node->cells[SE], NW, WEST,
					     node->cells[SW], NE, EAST);
    /* do west */
    if (node->cells[NW] != node->cells[SW])
      node->adjacent[W] = get_adjacent_node(node, node->cells[NW], SE, SOUTH,
					     node->cells[SW], NE, NORTH);
  }
}

/* given two cells (cell1, cell2), which share a common node, find 
   the other nodes along their common edge and determine which is closer
   to the common node */
G_nodes *get_adjacent_node(node, cell1, node_dir1, edge_dir1, cell2, 
			   node_dir2, edge_dir2)
     G_nodes *node;             /* the shared node */
     Gcell *cell1, *cell2;
     char edge_dir1, edge_dir2; /*which edge of each cell is the shared edge*/
     char node_dir1, node_dir2; /*which node of each cell is the shared one */
{
  G_nodes *node1 = NULL, *node2 = NULL;

  if (cell1 != NULL)
    node1 = get_other_gnode(cell1, edge_dir1, node_dir1);

  if (cell2 != NULL)
    node2 = get_other_gnode(cell2, edge_dir2, node_dir2);
  
  if (node1 == NULL && node2 != NULL)
    return node2;
  else if (node2 == NULL && node1 != NULL)
    return node1;
  
  /* the direction the node is from the cell (node_dir) happens to be the 
     same direction we want to look for the other cell */
  if (cell2 == node1->cells[node_dir1])
    return node1;
  else if (cell1 == node2->cells[node_dir2])
    return node2;
  else
    GP_PANIC("Impossible! can't determine node adjacency!");
}

/* get the other node on the edge_dir edge of cell */
G_nodes *get_other_gnode(cell, edge_dir, node_dir)
     Gcell *cell;
     char edge_dir, node_dir;
{
  G_nodes **nodes;

  if (cell == NULL)
    GP_PANIC("get_other_gnode given a null cell!");
    
  if (!is_leaf(cell))
    GP_PANIC("get_other_gnode requested other node of a nonleaf node!");

  nodes = cell->bndry.nodes;
  if (edge_dir == NORTH) {
    if (node_dir == NW)
      return nodes[NE];
    else if (node_dir == NE)
      return nodes[NW];
  }
  else if (edge_dir == SOUTH) {
    if (node_dir == SW)
      return nodes[SE];
    else if (node_dir == SE)
      return nodes[SW];
  }
  else if (edge_dir == EAST) {
    if (node_dir == SE)
      return nodes[NE];
    else if (node_dir == NE)
      return nodes[SE];
  }
  else if (edge_dir == WEST) {
    if (node_dir == SW)
      return nodes[NW];
    else if (node_dir == NW)
      return nodes[SW];
  }
  else
    GP_PANIC("get_other_node: bad edge_dir, node_dir combination");
  
}

compute_z_fils(gp)
     Nonuni_gp *gp;
{
  double *z_c, *thick, *z_pts, thickness;
  int num_z_pts;
  int i;

  thick = gp->thick = (double *)gp_malloc(gp->num_z_pts * sizeof(double));
  z_c = gp->z_c = (double *)gp_malloc(gp->num_z_pts * sizeof(double));
  z_pts = gp->z_pts;
  thickness = gp->grndp->thick;
  num_z_pts = gp->num_z_pts;

  if (num_z_pts == 1) {
    thick[0] = thickness;
    z_c[0] = 0.5*thickness;
    return;
  }
  
  if (z_pts[0] != 0.0 || z_pts[num_z_pts - 1] != 1.0)
    GP_PANIC("first and last z_pts should be 0 and 1");

  i = 0;
  thick[i] = (z_pts[i+1] - z_pts[i])/2.0;
  z_c[i] = z_pts[i] + thick[i]/4.0;

  for(i = 1; i < num_z_pts - 1; i++) {
    thick[i] = (z_pts[i+1] - z_pts[i-1])/2.0;

    fprintf(stderr,"This isn't right, you didn't fix this\n");
    z_c[i] = z_pts[i - 1] + thick[i]/2.0;
  }

  /* last one, like the first */
  thick[i] = (z_pts[i] - z_pts[i - 1])/2.0;
  z_c[i] = z_pts[i] - thick[i]/4.0;

  for(i = 0; i < num_z_pts; i++) {
    z_c[i] *= thickness;
    thick[i] *= thickness;
  }
}
  

/* this generates the segments for FastH */
generate_segs(gp, indsys)
     Nonuni_gp *gp;
     SYS *indsys;
{
  G_nodes *node, *othernode;
  double thickness, leftwidth, rightwidth;
  static complain = 0;
  static complain2 = 0;
  double x_shift, x_width, y_shift, y_width;
  
  gp->num_seg_groups = 0; 

  /* create segments in the north and east directions, and down(z),if needed */
  for(node = gp->nodelist; node != NULL; node = node->next) {
    /* do north */
    if (node->adjacent[N] != NULL) {
      /* width of segment and shift of center point. X_DIR is width dir */
      get_width_and_shift(X_DIR, node, node->cells[NW], node->cells[NE], 
			   &x_width, &x_shift);
      /* if this segment isn't a hole on both sides */
      if (x_width != 0) {
        make_segs(N, node, node->adjacent[N], x_width, x_shift, 0.0, gp, indsys);
        gp->num_seg_groups++;
      }
    }
    if (node->adjacent[E] != NULL) {
      get_width_and_shift(Y_DIR, node, node->cells[SE], node->cells[NE],
			  &y_width, &y_shift);
      if (y_width != 0) {
        make_segs(E, node, node->adjacent[E], y_width, 0.0, y_shift, gp, indsys);
        gp->num_seg_groups++;
      }
    }

    /* use this for z-directed segs later */
    node->x_shift = x_shift;
    node->y_shift = y_shift;
  }

  /* now do z-directed since we have all the x-segs */
  if (gp->num_z_pts != 1) {
    for(node = gp->nodelist; node != NULL; node = node->next) {
      /* - get south and west widths and shifts */
      /* - pick smaller of n/s and e/w for width and height of z dir segs */
      GP_PANIC("Z-directed segs not supported!");
      /*
      if (complain2++ == 0)
	fprintf(stderr, "\nYou didn't do z-dir yet\n\n");
	*/
    }
  }
}

/* figure out width of segment and where its center point will be
   (center point is given as a shift from node point) */
get_width_and_shift(width_dir, node, leftcell, rightcell, ret_width, 
		     ret_shift)
     char width_dir;
     G_nodes *node;
     Gcell *leftcell, *rightcell;
     double *ret_width, *ret_shift;
{
  double x_min, x_max;  /* put the two cells in a box and these are the
			      extremal values of the enclosing box*/
  double center;
  double node_x;
    
  if (width_dir == X_DIR) {
    get_x_cell_vals(leftcell, node, rightcell, &x_min, &x_max);
    node_x = node->x;
  }
  else if (width_dir == Y_DIR) {
    /* really y_left and y_right */
    get_y_cell_vals(leftcell, node, rightcell, &x_min, &x_max);
    node_x = node->y;
  }
  else
    GP_PANIC("Bad width_dir in get_width_and_center()");

  /* half the width of the combined width of the adjacent cells */
  *ret_width = (x_max - x_min)/2.0;

  /* center of this half-width section */
  center = node_x / 2.0 + (x_max + x_min) / 4.0;

  /* shift from node */
  *ret_shift = center - node_x; 
}  

get_x_cell_vals(left, node, right, x_left, x_right)
     Gcell *left, *right;
     G_nodes *node;
     double *x_left, *x_right;
{
  if (left != NULL && !is_hole(left))
    *x_left = get_x0(left);
  else
    *x_left = node->x;
  
  if (right != NULL && !is_hole(right))
    *x_right = get_x1(right);
  else 
    *x_right = node->x;
  
#if GP_DEBUG == TRUE
  /* do some data checking */
  if (left != NULL && right != NULL)
    if (get_x1(left) != node->x || get_x0(right) != node->x)
      GP_PANIC("Node doesn't seem to be at the right cell point!");
#endif
  
}

/* left and right of direction */
get_y_cell_vals(left, node, right, y_min, y_max)
     Gcell *left, *right;
     G_nodes *node;
     double *y_min, *y_max;
{
  if (left != NULL && !is_hole(left))
    *y_min = get_y0(left);
  else
    *y_min = node->y;
  
  if (right != NULL && !is_hole(right))
    *y_max = get_y1(right);
  else 
    *y_max = node->y;
  
#if GP_DEBUG == TRUE
  /* do some data checking */
  if (left != NULL && right != NULL)
    if (get_y1(left) != node->y || get_y0(right) != node->y)
      GP_PANIC("Node doesn't seem to be at the right cell point!");
#endif
  
}

make_segs(direction, node, othernode, width, x_shift, y_shift, gp, indsys)
     char direction;
     G_nodes *node, *othernode;
     double width, x_shift, y_shift;
     Nonuni_gp *gp;
     SYS *indsys;
{
  int num_z_pts = gp->num_z_pts;
  double *z_c = gp->z_c;
  double *thick = gp->thick;
  SEGMENT **segs;
  NODES *node0, *node1;

  double hx,hy,hz,wx,wy,wz;
  int i;

  if (direction == N) {
    hx = hy = 0;
    hz = 1;
    wx = 1;
    wy = 0;
    wz = 0;
    node->n_segs = (SEGMENT **)gp_malloc(num_z_pts*sizeof(SEGMENT *));
    segs = node->n_segs;
  }
  else if (direction == E) {
    hx = hy = 0;
    hz = 1;
    wx = 0;
    wy = 1;
    wz = 0;
    node->e_segs = (SEGMENT **)gp_malloc(num_z_pts*sizeof(SEGMENT *));
    segs = node->e_segs;
  }
  else
    GP_PANIC("make_segs: bad direction");

  for(i = 0; i < num_z_pts; i++) {
    /*
    draw_one_seg(direction,
		 node->x + x_shift, node->y + y_shift, z_c[i], 
		 othernode->x + x_shift, othernode->y + y_shift, z_c[i], 
		 width, thick[i], hx, hy, hz, wx, wy, wz,
		 gp->n_hinc, gp);
    */
    segs[i] = 
      make_one_seg(node->x + x_shift, node->y + y_shift, z_c[i], 
		 othernode->x + x_shift, othernode->y + y_shift, z_c[i], 
		 width, thick[i], wx, wy, wz, gp, indsys, &node0, &node1);
    node0->gp_node = node;
    node1->gp_node = othernode;
  }
}

SEGMENT *make_one_seg(x0,y0,z0,x1,y1,z1, width, height, wx, wy, wz, gp, indsys,
		      pnode0, pnode1)
     double x0, y0, z0, x1, y1, z1, width, height, wx, wy, wz;
     Nonuni_gp *gp;
     SYS *indsys;  /* needed in makeseg at the very least */
     NODES **pnode0, **pnode1;
{
  NODES *node0, *node1;
  static double dontcare = 2.0;
  double *widthdir;
  char name[100];
  NODES *makenode();
  double xg,yg,zg;
  double wxg, wyg, wzg;
  
  get_global_coords(x0,y0,z0,gp, &xg, &yg, &zg);
  sprintf(name, "%s_%d", gp->grndp->name, indsys->num_nodes);
  node0 = makenode(name, indsys->num_nodes, xg, yg, zg, GPTYPE, indsys);
  *pnode0 = node0;

  get_global_coords(x1,y1,z1,gp, &xg, &yg, &zg);
  sprintf(name, "%s_%d", gp->grndp->name, indsys->num_nodes);
  node1 = makenode(name, indsys->num_nodes, xg, yg, zg, GPTYPE, indsys);
  *pnode1 = node1;

  node0->gp = gp->grndp;
  node1->gp = gp->grndp;

  get_global_vec(wx, wy, wz, gp, &wxg, &wyg, &wzg);

  widthdir = (double *)gp_malloc(3*sizeof(double));
  widthdir[0] = wxg;
  widthdir[1] = wyg;
  widthdir[2] = wzg;

  sprintf(name, "%s_%d_%d", gp->grndp->name, node0->number , node1->number);

  return makeseg(name, node0, node1, height, width, gp->grndp->sigma, 
		 gp->grndp->hinc, 1, gp->grndp->rh, dontcare, widthdir, 
		 indsys->num_segs, GPTYPE, indsys); 

}

draw_one_seg(direction,
	     x1, y1, z1, x2, y2, z2, width, height, hx, hy, hz, wx, wy, wz,
	     nhinc, gp)
     char direction;
     double x1, y1, z1, x2, y2, z2, width, height;
     double hx, hy, hz, wx, wy, wz;
     int nhinc;
     Nonuni_gp *gp;
{

  FILE *fp = stdout;
  int i;
  double x,y,z;

  i = direction;

      /* first node end */
    fprintf(fp, "Q %d ",i);
    x = x1 - 0.5 * wx * width + 0.5 * hx * height;
    y = y1 - 0.5 * wy * width + 0.5 * hy * height;
    z = z1 - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);

    x = x1 - 0.5 * wx * width - 0.5 * hx * height;
    y = y1 - 0.5 * wy * width - 0.5 * hy * height;
    z = z1 - 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    x = x1 + 0.5 * wx * width - 0.5 * hx * height;
    y = y1 + 0.5 * wy * width - 0.5 * hy * height;
    z = z1 + 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = x1 + 0.5 * wx * width + 0.5 * hx * height;
    y = y1 + 0.5 * wy * width + 0.5 * hy * height;
    z = z1 + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    fprintf(fp,"\n");

    /* end at other node */
    fprintf(fp, "Q %d ",i);
    x = x2 - 0.5 * wx * width + 0.5 * hx * height;
    y = y2 - 0.5 * wy * width + 0.5 * hy * height;
    z = z2 - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);
  
    x = x2 - 0.5 * wx * width - 0.5 * hx * height;
    y = y2 - 0.5 * wy * width - 0.5 * hy * height;
    z = z2 - 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    x = x2 + 0.5 * wx * width - 0.5 * hx * height;
    y = y2 + 0.5 * wy * width - 0.5 * hy * height;
    z = z2 + 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = x2 + 0.5 * wx * width + 0.5 * hx * height;
    y = y2 + 0.5 * wy * width + 0.5 * hy * height;
    z = z2 + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
    fprintf(fp,"\n");

      /* left side */
      fprintf(fp, "Q %d ",i);
      x = x1 - 0.5 * wx * width + 0.5 * hx * height;
      y = y1 - 0.5 * wy * width + 0.5 * hy * height;
      z = z1 - 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z);

      x = x1 - 0.5 * wx * width - 0.5 * hx * height;
      y = y1 - 0.5 * wy * width - 0.5 * hy * height;
      z = z1 - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 

      x = x2 - 0.5 * wx * width - 0.5 * hx * height;
      y = y2 - 0.5 * wy * width - 0.5 * hy * height;
      z = z2 - 0.5 * wz * width - 0.5 * hz * height;
  
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = x2 - 0.5 * wx * width + 0.5 * hx * height;
      y = y2 - 0.5 * wy * width + 0.5 * hy * height;
      z = z2 - 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z);

      fprintf(fp,"\n");

    /* right side */
      fprintf(fp, "Q %d ",i);
      x = x1 + 0.5 * wx * width - 0.5 * hx * height;
      y = y1 + 0.5 * wy * width - 0.5 * hy * height;
      z = z1 + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
      x = x1 + 0.5 * wx * width + 0.5 * hx * height;
      y = y1 + 0.5 * wy * width + 0.5 * hy * height;
      z = z1 + 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
      x = x2 + 0.5 * wx * width + 0.5 * hx * height;
      y = y2 + 0.5 * wy * width + 0.5 * hy * height;
      z = z2 + 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
      x = x2 + 0.5 * wx * width - 0.5 * hx * height;
      y = y2 + 0.5 * wy * width - 0.5 * hy * height;
      z = z2 + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      fprintf(fp, "\n");
  
    /* top */
    fprintf(fp, "Q %d ",i);
    x = x1 - 0.5 * wx * width + 0.5 * hx * height;
    y = y1 - 0.5 * wy * width + 0.5 * hy * height;
    z = z1 - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);
   
    x = x1 + 0.5 * wx * width + 0.5 * hx * height;
    y = y1 + 0.5 * wy * width + 0.5 * hy * height;
    z = z1 + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = x2 + 0.5 * wx * width + 0.5 * hx * height;
    y = y2 + 0.5 * wy * width + 0.5 * hy * height;
    z = z2 + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
    x = x2 - 0.5 * wx * width + 0.5 * hx * height;
    y = y2 - 0.5 * wy * width + 0.5 * hy * height;
    z = z2 - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);

    fprintf(fp,"\n");

      /* bottom */
      fprintf(fp, "Q %d ",i);

      x = x1 - 0.5 * wx * width - 0.5 * hx * height;
      y = y1 - 0.5 * wy * width - 0.5 * hy * height;
      z = z1 - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 

      x = x1 + 0.5 * wx * width - 0.5 * hx * height;
      y = y1 + 0.5 * wy * width - 0.5 * hy * height;
      z = z1 + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = x2 + 0.5 * wx * width - 0.5 * hx * height;
      y = y2 + 0.5 * wy * width - 0.5 * hy * height;
      z = z2 + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = x2 - 0.5 * wx * width - 0.5 * hx * height;
      y = y2 - 0.5 * wy * width - 0.5 * hy * height;
      z = z2 - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
    fprintf(fp,"\n");
}     

print_cell_and_kids(cell)
{
  fprint_cell_and_kids(cell, stdout);
}

fprint_cell_and_kids(cell, fp)
     Gcell *cell;
     FILE *fp;
{
  dump_cell(cell, fp);
  
  switch(get_children_type(cell)) {
  case NONE:
    break;
  case BI:
    fprint_bi_kids( (Bi *)cell->children, fp);
    break;
  case GRID_2D:
    GP_PANIC("print_cell_and_kids:You haven't implemented grid_2d yet")
    break;
  default:
    GP_PANIC("Unknown child type in print_cell_and_kids")
    break;
  }
}

fprint_bi_kids(two_kids, fp)
     Bi *two_kids;
     FILE *fp;
{
  fprint_cell_and_kids(two_kids->child1, fp);
  fprint_cell_and_kids(two_kids->child2, fp);
}

dump_cell(cell, fp)
     Gcell *cell;
     FILE *fp;
{
  int i;

  fprintf(fp, "\nSelf: %d\t parent %d\t (x0 y0)(x1 y1) (%lg,%lg)(%lg,%lg)\n",
	  cell->index, cell->parent->index, 
	  cell->x0,cell->y0,cell->x1,cell->y1);

  fprintf(fp, "  area %lg \tChildren type: %d\tChildren: ",
	  (cell->y1-cell->y0)*(cell->x1 - cell->x0), get_children_type(cell));

  switch(get_children_type(cell)) {
  case NONE:
    fprintf(fp, "NONE\n");
    break;
  case BI:
    print_bi_addresses( (Bi *)cell->children, fp);
    break;
  case GRID_2D:
    GP_PANIC("Not implemented GRID_2D in dump_cell");
    break;
  default:
    GP_PANIC("Unknown child type in dump_cell");
    break;
  }
  
  if (is_leaf(cell)) {
    fprintf(fp, "   Nodes ne,se,sw,nw: ");
    for (i = 0; i < NUMNODES; i++)
      DUMP_INDEX(cell->bndry.nodes[i]);
    fprintf(fp, "\n");
  }
  
  fflush(fp);
}

print_bi_addresses(two_kids, fp)
     Bi *two_kids;
     FILE *fp;
{
  fprintf(fp, "%d %d\n",two_kids->child1->index, two_kids->child2->index);
}

print_node_list(node)
{
  fprint_node_list(node, stdout);
}

fprint_node_list(node, fp)
     G_nodes *node;
     FILE *fp;
{
  while(node != NULL) {
    dump_node(node, fp);
    node = node->next;
  }
}

dump_node(node,fp)
     FILE *fp;
     G_nodes *node;
{
  int i;

  fprintf(fp, "\nself: %d \tx,y (%lg,%lg)\t flag %d\tnesw adj: ", 
	  node->index,node->x, node->y,node->flag);
  for (i = 0; i < NUMADJ; i++)
    DUMP_INDEX(node->adjacent[i]);

  fprintf(fp,"\n cells ne,se,sw,nw ");
  for (i = 0; i < NUM_N_CELLS; i++)
    DUMP_INDEX(node->cells[i]);

  fprintf(fp, "\tprev ");
  DUMP_INDEX(node->prev);
  fprintf(fp, "\tnext ");
  DUMP_INDEX(node->next);
  fprintf(fp, "\n");
  
  /*fprintf(fp,"\tprev %d\tnext %d\n",node->prev->index, node->next->index);*/
}


debug_func()
{
  int i;

  i = 0;
}

dump_leaf_cells_to_file(cell, fname)
     Gcell *cell;
     char *fname;
{
  FILE *fp;

  if ( (fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "dump_leaf_cells_to_file: couldn't open file\n");
    return;
  }
    
  fprintf(fp,"0 A nonuniform plane hierarchy: %s\n",fname);
  dump_leaf_cells(cell, fp);
  fclose(fp);
}

dump_leaf_cells(cell, fp)
     Gcell *cell;
     FILE *fp;
{
  switch (get_children_type(cell)) {
  case NONE:
    print_leaf_cell(cell, fp);
    break;
  case BI:
    dump_leaf_cells( ((Bi *)cell->children)->child1, fp);
    dump_leaf_cells( ((Bi *)cell->children)->child2, fp);
    break;
  case GRID_2D:
    dump_grid_leaf_cells( (Grid_2d *)cell->children, fp);
    /* GP_PANIC("dump_leaf_cells:You haven't implemented grid_2d yet");*/
    break;
  default:
    GP_PANIC("Unknown child type in dump_leaf_cells")
    break;
  }
}

dump_grid_leaf_cells( grid, fp )
     Grid_2d *grid;
     FILE *fp;
{
  int i, j;

  for(i = 0; i < grid->y_cells; i++)
    for(j = 0; j < grid->x_cells; j++)
      dump_leaf_cells(grid->kids[i][j], fp);
}

print_leaf_cell(cell, fp)
     Gcell *cell;
     FILE *fp;
{
  fprintf(fp, "Q %d  %lg %lg %lg  %lg %lg %lg  %lg %lg %lg  %lg %lg %lg\n",
	  cell->index, cell->x0, cell->y0, 0.0, cell->x1, cell->y0, 0.0,
	  cell->x1, cell->y1, 0.0, cell->x0, cell->y1, 0.0);
}

dump_nonuni_plane_currents(gp, Ib, fp)
     Nonuni_gp *gp;
     CX *Ib;
     FILE *fp;
{
  G_nodes *gnode;
  double x,y,z;

  for(gnode = gp->nodelist; gnode != NULL; gnode=gnode->next) {
    get_global_coords(gnode->x, gnode->y, gp->z0, gp, &x, &y, &z);
    if (gnode->e_segs != NULL && gnode->n_segs != NULL)
      fprintf(fp, "%lg %lg %lg   %lg +j %lg    %lg +j %lg\n",
	      x, y, z, 
	      Ib[gnode->e_segs[0]->filaments[0].filnumber].real,
	      Ib[gnode->e_segs[0]->filaments[0].filnumber].imag,
	      Ib[gnode->n_segs[0]->filaments[0].filnumber].real,
	      Ib[gnode->n_segs[0]->filaments[0].filnumber].imag
	      );
    else if (gnode->e_segs != NULL)
      fprintf(fp, "%lg %lg %lg   %lg +j %lg    %lg +j %lg\n",
	      x, y, z, 
	      Ib[gnode->e_segs[0]->filaments[0].filnumber].real,
	      Ib[gnode->e_segs[0]->filaments[0].filnumber].imag,
	      0.0, 0.0);
    else if (gnode->n_segs != NULL)
      fprintf(fp, "%lg %lg %lg   %lg +j %lg    %lg +j %lg\n",
	      x, y, z, 
	      0.0,
	      0.0,
	      Ib[gnode->n_segs[0]->filaments[0].filnumber].real,
	      Ib[gnode->n_segs[0]->filaments[0].filnumber].imag
	      );
  }
}

/* return effective perimeter of a node in the plane */
/* This is essentially the area that can be used to calculate 
   the current density from an injection point and should match 
   the real perimeter of the contact for accurate resistance calculation */
double get_perimeter(node)
     G_nodes *node;
{
  double sum = 0;

  if (node->n_segs != NULL)
    sum += node->n_segs[0]->width;

  if (node->e_segs != NULL)
    sum += node->e_segs[0]->width;

  if (node->adjacent[S] != NULL)
    sum += node->adjacent[S]->n_segs[0]->width;

  if (node->adjacent[W] != NULL)
    sum += node->adjacent[W]->e_segs[0]->width;

  return sum;
}
