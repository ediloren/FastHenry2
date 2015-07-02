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

#ifndef _GP_H
#define _GP_H

#include "induct.h"

#define TRUE 1
#define FALSE 0

#define CLEAR 0
#define NOCLEAR 1

#define NOT_CHECKED 0
#define CHECKED 1

#define ALIVE 0
#define DEAD 1

#define X_DIR 0
#define Y_DIR 1

/* number of cells a node can be a part of */
#define NUM_N_CELLS 4

/* number of cells an edge can be a part of */
#define NUM_E_CELLS 2

/* number of edges a cell has and their directions */
#define NUMEDGES 4
#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3

/* number of nodes to a cell and their positions */
#define NUMNODES 4
#define NE 0
#define SE 1
#define SW 2
#define NW 3

#define NUMADJ 4
#define N 0
#define E 1
#define S 2
#define W 3

typedef struct llist {
  void *ptr;
  struct llist *next;
} Llist;

/*   
   In FastHenry 2.0, a plane was a grid of segments in both x and y where
   the four segments making a square we call a 'cell'.
   Now let's assume our plane is made up of a grid of these cells, and
   now we allow a cell to be broken into subcells in a hierarchical fashion
   for nonuniform discretizations */ 

typedef struct nonuni_gp {
  struct gcell *root_cell;  /* root of the discretization tree */

  int num_nodes;
  int num_cells;
  int num_leaves; /* number of leaves == number of meshes if num_z_pts=1 */
  int num_seg_groups; /* number of segs if no z-directed segs */

  /* define relative coordinate system for plane so that top and bottom
     of plane are in x-y plane */

  /* coordinates of one corner of 3D rectangular solid  in real coords */
  double x0, y0, z0;   /* origin_p */

  /* The plane axis directions */
  double ux_x, ux_y, ux_z;            /* ux=unit vector along gp x direction */
  double uy_x, uy_y, uy_z;            /* uy=unit vector along gp y direction */
  /* So, given plane coords (xp,yp), real coord is:  
     origin_p + xp * ux + yp * uy    */
  
  double uz_x, uz_y, uz_z;            /* unit vector along gp z dir */
  /*double thickness; */
  
  /* describe discretization in z direction */
  int num_z_pts;   
  double *z_pts;   /*on the interval [0,1] describe where points should be*/

       /* so if num_z_pts = 3,  there will be 3 layers of x-y segs and two
	  zdirected segments connecting the layers.  Example:
          z_pts[0] = 0, z_pts[1] = 0.5, z_pts[2] = 1, then the z-directed
	  segments will go from 0 - 0.5, and 0.5 - 1. where 0 is 
	  the bottom.
       */

  /* also the center points for filaments and their thicknesses */
  double *thick;
  double *z_c;

  /* linked list of all the nodes */
  struct g_nodes *nodelist;

  /* fastH old ground plane */
  struct Groundplane *grndp;

  /* has the cell edge info been corrupted because i didn't think i'd need it*/
  char is_edge_corrupted;

} Nonuni_gp;

/* hierarchal cells to represent a nonuniformly discretized plane of 
   material. */

/* a node in the tree */
typedef struct gcell {
  
  int index;            /* an index for debugging */
  void *children;       /* a pointer to a structure representing the subcells*/
  char children_type;    /* what is the structure for the children */
  /* the children types: */
#ifdef NONE
#undef NONE
#endif
#define NONE 0x0
#define BI 0x1         /* cell is divided in half */
#define QUAD 0x2        /* cell is divided into four */
#define GRID_2D 0x3     /* a 2D grid of cells */

  char ishole;        /* does this cell represent a hole */

  struct gcell *parent; /* The parent of this cell */

  /* in the root cell coorinate frame, here are two corners of the cell */
  double x0, y0;    /* corner closest to origin (SW corner) */
  double x1, y1;    /* corner farthest from origin (NE corner) */

  union {
    struct g_edges *edges[NUMEDGES]; /* struct describing the edge of a cell */
    struct g_nodes *nodes[NUMNODES]; /* nodes of cell (only if leaf!) */
  } bndry;

} Gcell;


/* a edge of a cell holds info for the adjacent cell. This structure holds
  the info in a tree like fashion */

typedef struct g_edges {
  
  struct gcell *cells[NUM_E_CELLS]; /* cells for which this is a boundary */
  char dirs[NUM_E_CELLS];           /* which direction this edge is within 
				    corresponding cell */

  /* a pointer to a structure holding children making up this edge
     if it is not a leaf */

  int num_children;        /* for an edge, this must be a binary tree! */
                           /*  or could be unary */
  Gcell **children;          /* points to cells */

} G_edges;

typedef struct g_nodes {
  
  int index;      /* an index for debugging */
  double x, y;
  
  struct gcell *cells[NUM_N_CELLS];

  /* need all four because of z-directed segs */
  struct g_nodes *adjacent[NUMADJ];

  char flag;          /* flag to mark things */

  double x_shift, y_shift; /* shift for center of z-directed segments */
  
  /* the number of segs is determined by gp->num_z_pts */
  SEGMENT **e_segs;  /* array of segs in east direction */
  SEGMENT **n_segs;  /* in north direction */

  struct g_nodes *prev;
  struct g_nodes *next;

} G_nodes;

/* a binary tree */
typedef struct bi {

  struct gcell *child1;  /* either North or East child */
  struct gcell *child2;  /*    South or West */

  char type; /* divided North/South or East/West */
#define NS 1
#define EW 2

} Bi;

/* a quad tree */
typedef struct quad {

  struct gcell *kids[4];

} Quad;

/* a 2d array of children */
typedef struct grid_2d {
  
  struct gcell ***kids;
  int x_cells, y_cells;

} Grid_2d;

typedef struct {
  Nonuni_gp *gp;
} Info;

typedef struct _nonuni_choice_list {
  G_nodes *node;
  SEGMENT *seg;
  double rank;                    /* weight for this choice */
} nonuni_choice_list;

G_nodes *add_to_gnodelist();
G_edges *make_one_edge();
G_edges *make_two_edge();
G_nodes *get_other_gnode();
Gcell *new_Gcells();
Gcell *new_Gcell();
G_nodes *new_Gnode();
G_nodes *make_new_node();
void *gp_malloc();
G_nodes *get_adjacent_node();
SEGMENT *make_one_seg();
G_nodes *find_nearest_nonuni_node();
Gcell *get_containing_cell();
Gcell *get_containing_grid_cell();
Gcell *get_containing_bi_cell();
SPATH *get_a_nonuni_path();
Llist *add_ptr_to_list();
int add_nonuni_choice();
G_nodes *find_nearest_edge_node();
G_nodes *scan_edge();
double get_node_dist();
double get_dist();
MELEMENT *add_edge_segs_to_list();
double get_perimeter();
NODES *get_nonuni_node_from_list();
NODES *make_new_node_with_nonuni();
Llist *get_nodes_inside_rect();
Llist *bi_get_nodes_inside_rect();
Llist *grid_get_nodes_inside_rect();
Llist *which_nodes_inside();

/* contact.c */
G_nodes *find_mid_node();
Gcell *find_next_cell_along_line();
double edge_coord();
Gcell *cut_cell();
G_nodes *add_new_node();
NODES *make_fastH_node();
ContactList *make_contact_connection();
Gcell *pick_cell_based_on_vec();

#define get_x0(cell) (cell->x0)
#define get_y0(cell) (cell->y0)
#define get_x1(cell) (cell->x1)
#define get_y1(cell) (cell->y1)

/* could also check if children_type == NONE */
#define is_leaf(cell) (cell->children == NULL)
#define get_children_type(cell) (cell->children_type)
#define get_bi_type(two) (two->type)
#define is_hole(cell) (cell->ishole == TRUE)

#define GP_PANIC(str) { fprintf(stderr,"Internal error in nonuniform plane code: %s\n",str); debug_func(); exit(1); }

#define DUMP_INDEX(ptr) {if ((ptr) == NULL) \
			   fprintf(fp, "x "); \
			 else\
			   fprintf(fp, "%d ", (ptr)->index);}

/* #ifndef _GP_GLOBAL 
   #define _GP_GLOBAL
   int opposite[4] = {SW, NW, NE, SE}
   #else
   extern opposite[4];
   #endif
*/

/* return opposite direction to dir */
#define opposite_dir(dir) ( (dir + 2)%4 )

#endif /* _GP_H */
