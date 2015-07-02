#include "gp.h"
#include "induct.h"

SPATH *path_through_nonuni_gp(nodein, nodeout, plane)
     NODES *nodein, *nodeout;
     GROUNDPLANE *plane;
{
  G_nodes *in, *out;
  Nonuni_gp *gp = plane->nonuni;
  SPATH *path, *get_a_nonuni_path();

  if (nodein->gp_node != NULL)
    in = nodein->gp_node;
  else
    GP_PANIC("path_through_nonuni_gp: nodein->gp_node == NULL!");

  if (nodeout->gp_node != NULL)
    out = nodeout->gp_node;
  else 
    GP_PANIC("path_through_nonuni_gp: nodeout->gp_node == NULL!");

  if (in == out) {
    GP_PANIC("path_through_nonuni_gp: Path starts and ends at same node!");
    /*
    fprintf(stderr, "Warning: user gp node at (%lg, %lg, %lg) (meters) and node at (%lg, %lg, %lg)\n", nodein->x, nodein->y, nodein->z, nodeout->x, 
	    nodeout->y, nodeout->z);
    fprintf(stderr, "in Nonuniform gp %s correspond to the same node\n",
	    gp->grndp->name);
    return NULL;  
    */
  }

  clear_nonuni_marks(gp->nodelist);
  path = get_a_nonuni_path(in, out, gp, NULL);

  /* if z-directed segs, add them to list */
  if (gp->num_z_pts != 1)
    GP_PANIC("find_non_uni_gp: You haven't added z-segs to path");

  return path;
}

/* clear all the flags as NOT_CHECKED */
clear_nonuni_marks(node)
     G_nodes *node;
{
  while(node != NULL) {
    node->flag = NOT_CHECKED;
    node = node->next;
  }
}

/*find the nearest node in the leaf cell that contains this global xyz point */
G_nodes *find_nearest_nonuni_node(xg, yg, zg, gp)
     double xg, yg, zg;
     Nonuni_gp *gp;
{
  Gcell *nearest;
  G_nodes *node;

  double x,y,z;

  get_nonuni_coords(xg, yg, zg, gp, &x, &y, &z);

  /* find the cell this point is in */
  nearest = get_containing_cell(x, y, gp->root_cell);

  if (!is_in_cell(x,y,nearest))
    return NULL;

  /* now find nearest node along the edges of the cell */
  return find_nearest_edge_node(x, y, nearest);

}
     
SPATH *get_a_nonuni_path(node, node_goal, plane, nodes_so_far)
     G_nodes *node, *node_goal;
     Nonuni_gp *plane;
     Llist *nodes_so_far;
{
#define MAXchoices 4
  nonuni_choice_list choices[MAXchoices];

  int nchoices;
  SPATH *path;
  seg_ptr sptr;
  double distance_to_goal;
  double new_dist;
  G_nodes *neighbor;
  int i;

  distance_to_goal = get_node_dist(node_goal->x, node_goal->y, node);

  nchoices = 0;

  sptr.type = NORMAL;

  if (node->adjacent[N] != NULL && node->n_segs != NULL) {
    neighbor = node->adjacent[N];
    new_dist = get_node_dist(node_goal->x, node_goal->y, neighbor);
    nchoices += add_nonuni_choice(&choices[nchoices], nodes_so_far, 
		      node->n_segs, neighbor,
		      distance_to_goal - new_dist); 
  }
  if (node->adjacent[E] != NULL && node->e_segs != NULL) {
    neighbor = node->adjacent[E];
    new_dist = get_node_dist(node_goal->x, node_goal->y, neighbor);
    nchoices += add_nonuni_choice(&choices[nchoices], nodes_so_far, 
		      node->e_segs, neighbor,
		      distance_to_goal - new_dist); 
  }
  if (node->adjacent[S] != NULL && node->adjacent[S]->n_segs != NULL) {
    neighbor = node->adjacent[S];
    new_dist = get_node_dist(node_goal->x, node_goal->y, neighbor);
    nchoices += add_nonuni_choice(&choices[nchoices], nodes_so_far, 
		      node->adjacent[S]->n_segs, neighbor,
		      distance_to_goal - new_dist); 
  }
  if (node->adjacent[W] != NULL && node->adjacent[W]->e_segs != NULL) {
    neighbor = node->adjacent[W];
    new_dist = get_node_dist(node_goal->x, node_goal->y, neighbor);
    nchoices += add_nonuni_choice(&choices[nchoices], nodes_so_far, 
		      node->adjacent[W]->e_segs, neighbor,
		      distance_to_goal - new_dist); 
  }

  /* sort in descreasing order by rank */
  sort_nonuni_choices(choices, nchoices);
  
  /* add current node to list of passed through nodes */
  nodes_so_far = add_ptr_to_list( (void *)node, nodes_so_far);
  
  for(path = NULL, i = 0; path == NULL && i < nchoices; i++) {
    sptr.segp = (void *)choices[i].seg;
    if (node_goal == choices[i].node) {
      path = add_seg_to_list(sptr, path);
      /* increment_usage(choices[i].seg); not used anymore (for overlap) */
    }
    else {
      /* call myself with next node */
      path = get_a_nonuni_path(choices[i].node, node_goal, plane, 
			       nodes_so_far);
      if (path != NULL) {
	/* there is a path, so lets add this seg to it */
	path = add_seg_to_list(sptr,path);
	/* increment_usage(choices[i].seg); not used anymore */
      }
    }
  }

  /* free this node */
  free(nodes_so_far);

  if (path == NULL)
    /* mark it is a node from which no path was found */
    node->flag = CHECKED;

  return path;

  
}

/* sort choices DECREASING by rank.  Not the most efficient */
sort_nonuni_choices(choices, num)
     nonuni_choice_list *choices;
     int num;
{
  int i, j;
  nonuni_choice_list temp;

  for(i = 0; i < num - 1; i++)
    for(j = num - 1; j > i; j--)
      if (choices[j-1].rank < choices[j].rank) {
	temp = choices[j-1];
	choices[j-1] = choices[j];
	choices[j] = temp;
      }
}
	

Llist *add_ptr_to_list(ptr, list)
     void *ptr;
     Llist *list;
{
  Llist *new_llist;

  new_llist = (Llist *)gp_malloc(sizeof(Llist));

  new_llist->ptr = ptr;
  new_llist->next = list;

  return new_llist;
}

/* add or not a possible direction to search and rank it */
int add_nonuni_choice(choice, nodes_so_far, segs, node, delta)
     nonuni_choice_list *choice;
     Llist *nodes_so_far;
     SEGMENT **segs;
     G_nodes *node;
     double delta;
{
  if (segs != NULL && !is_ptr_in_list((void *)node, nodes_so_far)
      && node->flag == NOT_CHECKED) {
    choice->seg = segs[0];
    choice->node = node;
    choice->rank = delta;
    
    return 1;
  }
  else
    return 0;
}

is_ptr_in_list(ptr, list)
     void *ptr;
     Llist *list;
{
  while (list != NULL) {
    if (list->ptr == ptr)
      return TRUE;
    list = list->next;
  }
  
  return FALSE;
}
      
/* find the cell that contains this point */
Gcell *get_containing_cell(x, y, cell)
     double x,y;
     Gcell *cell;
{

  switch(get_children_type(cell)) {
  case NONE:
    /* we are a leaf and the point is in us */
    return cell;
  case BI:
    return get_containing_bi_cell( x, y, cell);
    break;
  case GRID_2D:
    return get_containing_grid_cell( x, y, cell );
    break;
  default:
    GP_PANIC("Unknown child type in set_cell_coords")
    break;
  }
}

/* figure out which of the grid of children contains (x,y) */
Gcell *get_containing_grid_cell( x, y, cell )
     double x, y;
     Gcell *cell;
{
  double xn, yn;  /* normalized coords */
  Grid_2d *grid;
  int x_i, y_i;

  grid = (Grid_2d *)cell->children;

  /* you've written a function to do this: get_grid_indices() */
  
  xn = (x - cell->x0)/(cell->x1 - cell->x0);
  yn = (y - cell->y0)/(cell->y1 - cell->y0);

  if (xn >= 1) 
    x_i = grid->x_cells - 1;
  else if (xn < 0)
    x_i = 0;
  else
    /* round to the nearest cell */
    x_i = ((int) (xn * grid->x_cells));

  if (yn >= 1)
    y_i = grid->y_cells - 1;
  else if (yn < 0)
    y_i = 0;
  else
    /* round to the nearest cell */
    y_i = ((int) (yn * grid->y_cells));

  return get_containing_cell(x,y,grid->kids[grid->y_cells - y_i - 1][x_i]);
}


/* figure out which of the two children of cell contains (x,y) */
Gcell *get_containing_bi_cell( x, y, cell)
     double x, y;
     Gcell *cell;
{
  double xn, yn;  /* normalized coords */
  Bi *two_kids;

  two_kids = (Bi *)cell->children;

  xn = (x - cell->x0)/(cell->x1 - cell->x0);
  yn = (y - cell->y0)/(cell->y1 - cell->y0);

  if (get_bi_type(two_kids) == NS) {
    if (yn > 0.5)
      return get_containing_cell(x,y,two_kids->child1);
    else
      return get_containing_cell(x,y,two_kids->child2);
  }
  else if (get_bi_type(two_kids) == EW) {
    if (xn > 0.5)
      return get_containing_cell(x,y,two_kids->child1);
    else
      return get_containing_cell(x,y,two_kids->child2);
  }
  else 
    GP_PANIC("get_containing bi_kids: Unknown bi_type!");
}

is_in_cell(x, y, cell)
     double x,y;
     Gcell *cell;
{
  if (x < cell->x0 || x > cell->x1 || y < cell->y0 || y > cell->y1)
    return FALSE;
  else
    return TRUE;
}

/* looks for the nearest node to x,y along the edges of the cell */
/* Called by find_nearest_nonuni_node() and make_equiv_rect() */
G_nodes *find_nearest_edge_node(x, y, cell)
     double x,y;
     Gcell *cell;
{
  double xn, yn;
  double x0 = get_x0(cell);
  double y0 = get_y0(cell);
  double x1 = get_x1(cell);
  double y1 = get_y1(cell);
  int is_y_bigger;

  G_nodes **nodes;
  int quad;  /* quadrant */
  char direction;

  if (!is_leaf(cell))
    GP_PANIC("find_nearest_edge_node: This isn't a leaf cell!");

  nodes = &(cell->bndry.nodes[0]);

  /* put in normalized coords: 1x1 square with origin at center */

  xn = (x - x0) / (x1 - x0)  - 0.5;
  yn = (y - y0) / (y1 - y0)  - 0.5;

  /* which node is this closest to */

  if (xn >= 0 && yn >= 0)
    quad = NE;
  else if (xn <= 0 && yn >= 0)
    quad = NW;
  else if (xn <= 0 && yn <= 0)
    quad = SW;
  else if (xn >= 0 && yn <= 0)
    quad = SE;

  /* which direction points to the next nearest node */
  
  if (fabs(xn) < fabs(yn))
    is_y_bigger = TRUE;
  else
    is_y_bigger = FALSE;

  /* scan along the given direction to find the nearest node */
  if (quad == NE) {
    if (is_y_bigger == TRUE)
      direction = W;
    else
      direction = S;
  }
  else if (quad == NW) {
    if (is_y_bigger == TRUE)
      direction = E;
    else
      direction = S;
  }
  else if (quad == SW) {
    if (is_y_bigger == TRUE)
      direction = E;
    else
      direction = N;
  }
  else if (quad == SE) {
    if (is_y_bigger == TRUE)
      direction = W;
    else
      direction = N;
  }
    
  return scan_edge(x, y, nodes[quad], direction);
}

G_nodes *scan_edge(x, y, node, dir)
     double x,y;
     G_nodes *node;
     char dir;
{

  /* scan adjacency direction from node on until we minimize distance */

  double dist, tmp;

  dist = get_node_dist(x,y,node);

  /* there better always be an adjacent node since the adjacent
     node to this one should be farther than the original node! */
  if (node->adjacent[dir] == NULL)
    GP_PANIC("scan_edge: adjacent[dir] == NULL!");

  while( (tmp = get_node_dist(x,y,node->adjacent[dir])) < dist) {
    dist = tmp;
    node = node->adjacent[dir];
  }

  return node;
}

double get_node_dist(x,y,node)
     double x,y;
     G_nodes *node;
{
  return get_dist(x - node->x, y - node->y);
}

double get_dist(x,y)
     double x,y;
{
  return sqrt(x*x + y*y);
}

int make_nonuni_Mlist(plane, pMlist)
     GROUNDPLANE *plane;
     MELEMENT **pMlist;
{
  int counter = 0;

  /* make meshes of the kids */
  make_children_meshes(plane->nonuni->root_cell, pMlist, &counter);

  if (plane->nonuni->num_z_pts != 1) {
    GP_PANIC("make_nonuni_Mlist: z_mesh function not written!");
    /* make_z_meshes(plane->nonuni, pMlist, &counter); */
  }

  if (counter > plane->numesh)
    /* can be less (from edge holes), just not greater */
    GP_PANIC("make_nonuni_Mlist: number of meshes too big!");

  return counter;
}

make_children_meshes(cell, pMlist, pcount)
     Gcell *cell;
     MELEMENT **pMlist;
     int *pcount;
{
  switch (get_children_type(cell)) {
  case NONE:
    *pcount += make_leaf_mesh(cell, &(pMlist[*pcount]));
    break;
  case BI:
    make_children_meshes( ( (Bi *)cell->children )->child1, pMlist, pcount);
    make_children_meshes( ( (Bi *)cell->children )->child2, pMlist, pcount);
    break;
  case GRID_2D:
    make_grid_children_meshes( (Grid_2d *)cell->children, pMlist, pcount);
    break;
  default:
    GP_PANIC("make_children_meshes: Unknown child type!");
    break;
  }
}

make_grid_children_meshes( grid, pMlist, pcount)
     Grid_2d *grid;
     MELEMENT **pMlist;
     int *pcount;
{
  int i, j;

  for(i = 0; i < grid->y_cells; i++)
    for(j = 0; j < grid->x_cells; j++)
      make_children_meshes( grid->kids[i][j], pMlist, pcount);
}
    
make_leaf_mesh(cell, pMlist)
     Gcell *cell;
     MELEMENT **pMlist;
{
  G_nodes **nodes;
  int bad;

  *pMlist = NULL;

  if (!is_leaf(cell))
    GP_PANIC("make_leaf_mesh: not a leaf");

  nodes = cell->bndry.nodes;

  if (nodes[SW]->n_segs == NULL || nodes[SW]->e_segs == NULL 
      || nodes[NW]->e_segs == NULL || nodes[SE]->n_segs == NULL) {
    /* this is a hole on an edge.  If it's not on an edge, we
       must report an error since we haven't implemented a hole consisting
       of multiple cells (only one mesh for the lot) */
    if (!is_hole(cell))
      GP_PANIC("make_leaf_mesh: NULL segments bordering non-hole cell!");

    bad = FALSE;
    if (nodes[SW]->n_segs == NULL && nodes[SW]->cells[NE] != NULL 
        && nodes[SW]->cells[NW] != NULL)
      bad = TRUE;
    if (nodes[SW]->e_segs == NULL && nodes[SW]->cells[NE] != NULL 
        && nodes[SW]->cells[SE] != NULL)
      bad = TRUE;
    if (nodes[NW]->e_segs == NULL && nodes[NW]->cells[NE] != NULL 
        && nodes[NW]->cells[SE] != NULL)
      bad = TRUE;
    if (nodes[SE]->n_segs == NULL && nodes[SE]->cells[NE] != NULL 
        && nodes[SE]->cells[NW] != NULL)
      bad = TRUE;

    if (bad == TRUE)
      GP_PANIC("make_leaf_mesh: Adjacent hole cells not allowed!");

    /* if not bad, then do nothing */
    return 0;
  }
  else {
    /* make a mesh as required */
    *pMlist = add_edge_segs_to_list(nodes[SW], nodes[NW], N, 1, *pMlist);
    *pMlist = add_edge_segs_to_list(nodes[NW], nodes[NE], E, 1, *pMlist);
    *pMlist = add_edge_segs_to_list(nodes[SE], nodes[NE], N, -1, *pMlist);
    *pMlist = add_edge_segs_to_list(nodes[SW], nodes[SE], E, -1, *pMlist);

    return 1;
  }

}

MELEMENT *add_edge_segs_to_list(node0, node1, dir, signofelem, mfirst)
     G_nodes *node0, *node1;
     char dir;
     int signofelem;
     MELEMENT *mfirst;
{
  G_nodes *node;
  MELEMENT *plist = mfirst;
  FILAMENT *fil;
  MELEMENT *melem;

  for(node = node0; node != node1; node = node->adjacent[dir]) {
    if (dir == N)
      fil = &(node->n_segs[0]->filaments[0]);
    else if (dir == E)
      fil = &(node->e_segs[0]->filaments[0]);
    else
      GP_PANIC("add_edge_segs_to_list: bad dir");
      
    melem = make_melement(fil->filnumber, fil, signofelem);
    plist = insert_in_list(melem, plist);
  }

  return plist;
  
}

NODES *get_or_make_nearest_node(name, index, x, y, z, indsys, gp, node_list)
     char *name;
     int index;
     double x,y,z;
     SYS *indsys;
     Nonuni_gp *gp;
     NPATH *node_list;
{
  G_nodes *nonuni_node;
  NODES *node;
  NPATH *nodeL;
  NODES *makenode();

  nonuni_node = find_nearest_nonuni_node(x, y, z, gp);

  if (nonuni_node == NULL) {
    fprintf(stderr,"Error: Node for plane %s requested at %lg %lg %lg doesn't appear to be within the plane\n",
	    gp->grndp->name, x, y, z);
    exit(1);
  }

  /* See if G_nodes node has already been made into a NODES node.
     If not, creaet it */
  
  /* look for nonuni_node already in tempnode list */
  node = get_nonuni_node_from_list(nonuni_node, node_list);

  if (node != NULL)
    /* this non-uni node is the same as another read so far */
    return node;
  else {
    /* make a new one */
    return make_new_node_with_nonuni(nonuni_node,name, index, 
                                     x, y, z, indsys, gp);
  }
}

/* make a new NORMAL node that is derived from a nonuni_node */
NODES *make_new_node_with_nonuni(nonuni_node, name, index, x, y, z, indsys, gp)
     G_nodes *nonuni_node;
     char *name;
     int index;
     double x,y,z;
     SYS *indsys;
     Nonuni_gp *gp;
{
  NODES *node;
  NODES *makenode();

  /* make a new one */
  node = makenode(name, index, x, y, z, NORMAL, indsys);
  node->gp = gp->grndp;
  node->gp_node = nonuni_node;
  return node;
}

/* find a NODES node in list nodeL whose gp_node is nonuni_node */
/* Note: you used to compare getrealnode(nodeL->node)->gp_node != nonuni_node
   but you took out getrealnode for equiv_rect. Hopefully that won't 
   break anything */
NODES *get_nonuni_node_from_list(nonuni_node, nodeL)
    G_nodes *nonuni_node;
    NPATH *nodeL;    
{
  
  while(nodeL != NULL && (nodeL->node)->gp_node != nonuni_node)
    nodeL = nodeL->next;

  if (nodeL != NULL)
    return nodeL->node;
  else
    return NULL;
}

/* return a linked list of all the nodes of the children of 
   this cell that are inside the rectangle with corners (x0,y0),(x1,y1) */
Llist *get_nodes_inside_rect(x0, y0, x1, y1, cell, endoflist)
     double x0,y0,x1,y1;
     Gcell *cell;
     Llist **endoflist;
{
  Llist *nodelist;

  /* check if there is any overlap for this child */
  if (intersection(x0,y0,x1,y1,cell->x0,cell->y0,cell->x1,cell->y1) == FALSE)
    return NULL;

  switch(get_children_type(cell)) {
  case NONE:
    /* we are a leaf. figure which nodes are inside */
    return which_nodes_inside(x0,y0,x1,y1,cell, endoflist);
  case BI:
    return bi_get_nodes_inside_rect( x0, y0, x1, y1, cell, endoflist);
    break;
  case GRID_2D:
    return grid_get_nodes_inside_rect( x0, y0, x1, y1, cell, endoflist);
    break;
  default:
    GP_PANIC("Unknown child type in set_cell_coords")
    break;
  }
}

Llist *bi_get_nodes_inside_rect( x0, y0, x1, y1, cell, endoflist)
    double x0, y0, x1, y1;
    Gcell *cell;
    Llist **endoflist;
{
  Bi *bi = (Bi *)(cell->children);
  Llist *one;
  Llist *two;
  Llist *oneend;

  one = get_nodes_inside_rect(x0, y0, x1, y1, bi->child1, &oneend);
  two = get_nodes_inside_rect(x0, y0, x1, y1, bi->child2, endoflist);

  
  /* concatenate two lists */
  if (one == NULL)
    one = two;
  else {
    if (oneend->next != NULL)
      GP_PANIC("bi_get_nodes_inside_rect: oneend is not end of list!");

    oneend->next = two;

    if (two == NULL)
      *endoflist = oneend;
  }

  /* endoflist assigned in two's call */
  return one;
  
}

Llist *grid_get_nodes_inside_rect( x0, y0, x1, y1, cell, endoflist)
    double x0, y0, x1, y1;
    Gcell *cell;
    Llist **endoflist;
{
  Grid_2d *grid = (Grid_2d *)cell->children;
  Llist *current;
  Llist *current_end;
  int row_start, row_end, col_start, col_end;
  int i, j;
  Llist *entire;
  
  *endoflist = NULL;

  /* to save a little cpu, let's figure out which cells intersect this
     rectangle */

  /* get indices that correspond to top left and bottom right corner */

  get_grid_indices(cell, x0, y1, &row_start, &col_start);  /* top left */
  get_grid_indices(cell, x1, y0, &row_end, &col_end);  /* bot right */

  entire = NULL;
  for(i = row_start; i <= row_end; i++)
    for(j = col_start; j <= col_end; j++) {
      current = get_nodes_inside_rect( x0, y0, x1, y1, grid->kids[i][j], 
                                       &current_end);
      if (current != NULL) {
        if (*endoflist == NULL)
          *endoflist = current_end;
        current_end->next = entire;
        entire = current;
      }

    }

  return entire;
}

/* returns the cell indices that contain (x,y) */
get_grid_indices(cell, x, y, pi, pj)
    Gcell *cell;
    double x,y;
    int *pi, *pj;
{
  Grid_2d *grid = (Grid_2d *)cell->children;
  double xn, yn;  /* normalized coords */
  int y_i;

  xn = (x - cell->x0)/(cell->x1 - cell->x0);
  yn = (y - cell->y0)/(cell->y1 - cell->y0);

  if (xn == 1) 
    *pj = grid->x_cells - 1;
  else
    /* round to the nearest cell */
    *pj = ((int) (xn * grid->x_cells));

  if (yn == 1)
    y_i = grid->y_cells - 1;
  else
    /* round to the nearest cell */
    y_i = ((int) (yn * grid->y_cells));

  *pi = grid->y_cells - y_i - 1;

}

/* returns TRUE if the two rectangles overlap.  
   lower left=(x0,y0) upper right = (x1,y1) */
intersection(x0,y0,x1,y1,cx0,cy0,cx1,cy1)
     double x0,y0,x1,y1,cx0,cy0,cx1,cy1;
{
  if ((x0 < cx0 && x1 < cx0) || (x0 > cx1 && x1 > cx1)
      && (y0 < cy0 && y1 < cy0) || (y0 > cy1 && y1 > cy1)) 
    return FALSE;
  else
    return TRUE;

}

/* returns a linked list of the nodes of this cell contained in the rectangle*/
Llist *which_nodes_inside(x0,y0,x1,y1,cell, endoflist)
    double x0, y0, x1, y1;
    Gcell *cell;
    Llist **endoflist;
{
  Llist *list = NULL;

  if (!is_leaf(cell))
    GP_PANIC("which_nodes_inside: non-leaf cell!");

  /* is SW inside */
  if (cell->x0 >= x0 && cell->x0 <= x1 && cell->y0 >= y0 && cell->y0 <= y1) {
    list = add_ptr_to_list((void *)cell->bndry.nodes[SW], list);
    if (list->next == NULL)
      *endoflist = list;
  }

  /* is SE inside */
  if (cell->x1 >= x0 && cell->x1 <= x1 && cell->y0 >= y0 && cell->y0 <= y1) {
    list = add_ptr_to_list((void *)cell->bndry.nodes[SE], list);
    if (list->next == NULL)
      *endoflist = list;
  }

  /* is NE inside */
  if (cell->x1 >= x0 && cell->x1 <= x1 && cell->y1 >= y0 && cell->y1 <= y1) {
    list = add_ptr_to_list((void *)cell->bndry.nodes[NE], list);
    if (list->next == NULL)
      *endoflist = list;
  }

  /* is NW inside */
  if (cell->x0 >= x0 && cell->x0 <= x1 && cell->y1 >= y0 && cell->y1 <= y1) {
    list = add_ptr_to_list((void *)cell->bndry.nodes[NW], list);
    if (list->next == NULL)
      *endoflist = list;
  }

  return list;
}

free_Llist(list)
     Llist *list;
{
  Llist *node;

  while(list != NULL) {
    node = list;
    list = list->next;
    free(node);
  }
}
