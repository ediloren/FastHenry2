
#include "induct.h"
#include "gp.h"

/* The following are support functions for contacts in a ground plane. */
/*  Similar to hole.c code  and uses many functions defined there*/
/*  To write new contact functions, add the name of your function call to
    make_contacts() and then include the function itself in here.
    */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "induct.h"

/* to compare floating point numbers */
#define nearzero(x) (fabs(x) < EPS)

/* SRW */
ContactList *make_contactlist(ContactList*, char*, double, double, double,
    double, int*);
void contact_error(char*, char*, ContactList*);
void contact_error2(char*, char*, char*);
void contact_warning(char*, ContactList*);
void regurg_contact(FILE*, ContactList*);
void make_contacts(ContactList*, Nonuni_gp*);
void contact_point(ContactList*, Nonuni_gp*, double, double, double, double);
void contact_line(ContactList*, Nonuni_gp*, double, double, double, double);
void walk_along_line(double, double, double, double, double, double,
    double, double, Nonuni_gp*);
Gcell *find_next_cell_along_line(double, double, double, double, Gcell*,
    double*, double*);
void get_new_x_y(double, double, double, double, Gcell*, double*, double*,
    G_nodes**, char*);
double edge_coord(double, double, double, double, double);
Gcell *cut_cell(double, double, double, double, Gcell*, Nonuni_gp*);
void break_cell(Gcell*, char, Nonuni_gp*);
void update_bi_nodes(Gcell*, Nonuni_gp*);
void clear_edge_ptrs(Gcell*);
void fix_node_cell_ptrs(Gcell*);
void set_edge_nodes(G_nodes*, G_nodes*, char, char, char, Gcell*);
void find_or_make_node(Gcell*, Gcell*, Gcell*, char, char, double, double,
    Nonuni_gp*);
void fix_adjacency(G_nodes*, char, G_nodes*, G_nodes*);
G_nodes *find_mid_node(G_nodes*, char, G_nodes*);
void make_two_kids(Gcell*, char, Nonuni_gp*);
void contact_rect(ContactList*, Nonuni_gp*, double, double, double, double);
void cut_inside_rect(double, double, double, double, double, double, double,
    Nonuni_gp*);
void contact_decay_rect(ContactList*, Nonuni_gp*, double, double, double,
    double);
void do_decay_rect(double, double, double, double, double, double, double,
    double, double, Nonuni_gp*);
void limit_box(double, double, double, double, double, double, double,
    double, double*, double*, double*, double*);
void compute_new_widths(double*, double*, double*, double*);
void contact_equiv_rect(ContactList*, Nonuni_gp*, double, double, double,
    double);
void make_equiv_rect(double, double, double, double, double, Nonuni_gp*,
    char*); 
void walk_along_edge(double, double, double, double, Nonuni_gp*, char,
    NODES*, char*, double);
void equiv_nodes_on_edge(G_nodes*, char, G_nodes*, NODES*, Nonuni_gp*,
    char*, double);
NODES *make_fastH_node(G_nodes*, GROUNDPLANE*, Nonuni_gp*, double, double,
    double, char*);
void contact_initial_grid(ContactList*, Nonuni_gp*, double, double, double,
    double);
void make_initial_grid(Nonuni_gp*, int, int);
void update_grid_nodes(Gcell*, Nonuni_gp*);
void set_node_and_cell_info(G_nodes*, int, Gcell*);
void set_cell_node_adjacency(Gcell*);
void point_at_each_other(G_nodes*, int, G_nodes*);
G_nodes *add_new_node(double, double, int, Gcell*, int, Nonuni_gp*);
void make_grid_kids(Gcell*, int, int, Nonuni_gp*);
void contact_initial_mesh_grid(ContactList*, Nonuni_gp*, double, double,
    double, double);
void poke_holes(Nonuni_gp*);
ContactList *make_contact_connection(ContactList*, char*, double, double,
    double, double, int*);
void contact_trace(ContactList*, Nonuni_gp*, double, double, double, double);
void do_trace(double, double, double, double, double, double, double, double,
    Nonuni_gp*);
Gcell *pick_cell_based_on_vec(G_nodes*, double, double);


#define MAXNAME 20

/* this turns the info in line into a ContactList element */
ContactList *make_contactlist(ContactList *head, char *line, double units,
    double relx, double rely, double relz, int *skip)
{
  char name[MAXNAME];
  char nname[80];
  int skip1, skip2, skip3, skip4, skip25;
  ContactList *contactp;
  char *linep;
  double *vals;
  
  name[MAXNAME-1] = '\0';
  sscanf(line, "%19s%n",name,&skip1);
  line+=skip1;

  if (strcmp("contact",name) != 0) {
    fprintf(stderr, "Internal Error:  %s is not the word 'contact'\n",name);
    exit(1);
  }

  name[MAXNAME-1] = '\0';
  sscanf(line, "%19s%n",name,&skip2);
  line+=skip2;

  if (name[MAXNAME-1] != '\0')
    printf("Warning: contact function '%s' truncated to 19 chars\n",name);

  *skip = skip1 + skip2;

  if (strcmp(name,"connection") == 0)
    /* a connection is made of multiple contact keywords, so call it 
       to make them and then call this function again */
    /* note: adds to skip */
    return make_contact_connection(head, line, units, relx, rely, relz, skip);

  contactp = (ContactList *)Gmalloc(sizeof(ContactList));
  
  contactp->func = (char *)Gmalloc((strlen(name) + 1)*sizeof(char));
  strcpy(contactp->func, name);
  contactp->units = units;
  contactp->relx = relx;
  contactp->rely = rely;
  contactp->relz = relz;

  /* some contacts have a name before the values */
  skip25 = 0;
  if (strcmp(name,"equiv_rect") == 0) {
    sscanf(line,"%s%n",nname,&skip25);
    line += skip25;
    contactp->name = (char *)Gmalloc((strlen(nname) + 1)*sizeof(char));
    strcpy(contactp->name, nname);
  }
    
  skip3 = 0;
  while(isspace(*line) && *line != '\0') {
    line++;
    skip3++;
  }

  if (*line != '(')
    contact_error("Values for contact must start with '('",line,contactp); 
  
  /* let's count how many values we have first */
  linep = line;
  contactp->numvals = 0;

  while(*linep != ')' && !eos(*linep)) {
    linep++;
    while(isspace(*linep) || is_one_of(*linep, "1234567890.e+-"))
      linep++;
    contactp->numvals++;
  }
  if (*linep != ')')
    contact_error("Contact values did not end with ')'", line, contactp);

  contactp->vals = (double *)Gmalloc(contactp->numvals*sizeof(double));

  skip3 += linep - line + 1;

  /* now let's read them in */
  linep = line + 1;  /* skip the ( */
  vals = contactp->vals;

  while(*linep != ')') {
    if (sscanf(linep, "%lf%n",vals,&skip4) != 1)
      contact_error("Couldn't read value starting at:",linep,contactp);
    else {
      linep += skip4;
      vals++;
      linep += skipspace(linep);
      if (*linep == ',')
	linep++;
    }
  }

  *skip = skip1 + skip2 + skip25 + skip3;

  contactp->next = head;
  contactp->done = FALSE;

  return contactp;

}


void contact_error(char *errstr, char *line, ContactList *contactp)
{
  fprintf(stderr, "While reading in a %s contact:\n",contactp->func);
  fprintf(stderr, "%s\n%s",errstr,line);
  fprintf(stderr, "The contact format is 'contact <type> (val1,val2,...)'\n");
  fprintf(stderr, "Where type is a string, and valn is a floating point value\n");
  exit(1);
}

void contact_error2(char *errstr, char *line, char *nametype)
{
  fprintf(stderr, "While reading in a %s contact:\n",nametype);
  fprintf(stderr, "%s\n%s",errstr,line);
  fprintf(stderr, "The contact format is 'contact <type> <nodename> (val1,val2,...)'\n");
  fprintf(stderr, "Where type is a string, and valn is a floating point value\n");
  exit(1);
}
  
void contact_warning(char *errstr, ContactList *contactp)
{
  fprintf(stderr, "Warning: While reading in a %s contact:\n",contactp->func);
  fprintf(stderr, "%s\n",errstr);
  regurg_contact(stderr, contactp);
}
  
void regurg_contact(FILE *fp, ContactList *contactp)
{
  int i;

  fprintf(fp,"contact %s (",contactp->func);
  for(i = 0; i < contactp->numvals; i++) {
    if (i != 0) fprintf(fp,",");
    fprintf(fp,"%lg",contactp->vals[i]);
  }
  fprintf(fp,")\n");
}

/* This calls the functions which make the discretization around the contacts*/
void make_contacts(ContactList *contactp, Nonuni_gp *gp)
{
  double relx = contactp->relx;
  double rely = contactp->rely;
  double relz = contactp->relz;
  double units = contactp->units;

  if (strncmp("rect",contactp->func,4) == 0)
    contact_rect(contactp, gp, relx, rely, relz, units);
  else if (strncmp("line",contactp->func,4) == 0)
    contact_line(contactp, gp, relx, rely, relz, units);
  else if (strncmp("point",contactp->func,5) == 0)
    contact_point(contactp, gp, relx, rely, relz, units);
  else if (strncmp("decay_rect",contactp->func,10) == 0)
    contact_decay_rect(contactp, gp, relx, rely, relz, units);
  else if (strncmp("equiv_rect",contactp->func,10) == 0)
    contact_equiv_rect(contactp, gp, relx, rely, relz, units);
  else if (strncmp("initial_grid",contactp->func,10) == 0)
    contact_initial_grid(contactp, gp, relx, rely, relz, units);
  else if (strncmp("initial_mesh_grid",contactp->func,10) == 0)
    contact_initial_mesh_grid(contactp, gp, relx, rely, relz, units);
  else if (strncmp("trace",contactp->func,5) == 0)
    contact_trace(contactp, gp, relx, rely, relz, units);
  else 
    contact_error("Unknown type of contact","",contactp);
}

/* cut up cell for this point until it's contained in a cell no
   bigger than the 4th argument */
void contact_point(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double xc, yc, zc;
  Gcell *contain;

  if (contactp->numvals != 5) 
    contact_error("Exactly 5 values required for a point contact.","",
		  contactp);

  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &xc, &yc, &zc);
  
  contain = get_containing_cell(xc, yc, gp->root_cell);
  
  if (!is_in_cell(xc,yc,contain)) 
    contact_error("Contact point not in plane!","",contactp);

  if (vals[3] <= 0.0 || vals[4] <= 0.0)
    contact_error("Widths in x-dir and y-dir must be greater than 0.","",
		  contactp);

  cut_cell(xc, yc, vals[3]*units, vals[4]*units, contain, gp);

}

/* follow a line, cutting every cell you come to */
/* vals[0]-[2] first point of line, vals[3]-[5] other end point, 
   vals[6] = max width in plane's x-dir,  vals[7] = max width in y-dir */
void contact_line(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double x0, y0, z0, x1, y1, z1;
  Gcell *contain;

  if (contactp->numvals != 8) 
    contact_error("Exactly 8 values required for a line contact.","",
		  contactp);

  /* beginning point */
  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &x0, &y0, &z0);

  /* end point */
  get_nonuni_coords(vals[3]*units + relx, vals[4]*units + rely, 
		    vals[5]*units + relz, gp, &x1, &y1, &z1);
  
  contain = get_containing_cell(x0, y0, gp->root_cell);
  
  if (!is_in_cell(x0,y0,contain)) 
    contact_error("First point of contact line not in plane!","",contactp);

  contain = get_containing_cell(x1, y1, gp->root_cell);
  
  if (!is_in_cell(x1,y1,contain)) 
    contact_error("Second point of contact line not in plane!","",contactp);

  if (vals[6] <= 0.0 || vals[7] <= 0.0)
    contact_error("Widths in x-dir and y-dir must be greater than 0.","",
		  contactp);

  /* walk along line and cut up cells */
  walk_along_line(x0, y0, z0, x1, y1, z1, vals[6]*units, vals[7]*units, gp);

}

/*walk along line from (x0,y0) to (x1,y1) and cut cells to size (max_x,max_y)*/
void walk_along_line(double x0, double y0, double z0, double x1, double y1,
    double z1, double max_x, double max_y, Nonuni_gp *gp)
{
  Gcell *contain, *smallcell;
  double xv, yv;  /* unit vector along line */
  double mag;
  double x,y, x_new, y_new;

  mag = sqrt( SQUARE(x1-x0) + SQUARE(y1 - y0));
  xv = (x1 - x0)/mag;
  yv = (y1 - y0)/mag;

  x = x0;
  y = y0;

  contain = get_containing_cell(x, y, gp->root_cell);

  while( (x - x0)*xv + (y - y0)*yv < mag*(1 - EPS) ) {
    if (contain == NULL)
      GP_PANIC("walk_along_line: contain==NULL. line out of plane, probably");

    /* cut down to this point and update contain */
    smallcell = cut_cell(x, y, max_x, max_y, contain, gp);

    /* find next cell and coordinates */
    contain = find_next_cell_along_line(x, y, xv, yv, smallcell, 
					&x_new, &y_new);

    x = x_new;
    y = y_new;
  }
}

Gcell *find_next_cell_along_line(double x, double y, double xv, double yv,
    Gcell *cell, double *ret_x, double *ret_y)
{
  char edge_dir;
  G_nodes *node;
  char dir;
  Gcell *newcell;

  /* we get new x,y and a node and dir to search for nearest node */
  get_new_x_y(x, y, xv, yv, cell, ret_x, ret_y, &node, &dir);

  /* we assume dir == N or E */
  if (dir == N) {

    if (node->y == *ret_y) {
      /* we've hit a node.   */
      newcell = pick_cell_based_on_vec(node, xv, yv);
    }
    else {
      /* find node north of the new y */
      while(node->y < *ret_y)
	node = node->adjacent[N];
      
      if (node->y == *ret_y)
        /* we're at a node */
        newcell = pick_cell_based_on_vec(node, xv, yv);
      else {
        /* return cell in appropriate direction */
        if (xv >= 0)
          newcell = node->cells[SE];
        else
          newcell = node->cells[SW];
      }
    }
  }
  else if (dir == E) {

    if (node->x == *ret_x) {
      newcell = pick_cell_based_on_vec(node, xv, yv);
    }
    else {
      /* find node east of new x */
      while(node->x < *ret_x)
	node = node->adjacent[E];
      
      if (node->x == *ret_x)
        /* we're at a node */
        newcell = pick_cell_based_on_vec(node, xv, yv);
      else {
        /* return cell in appropriate direction */
        if (yv >= 0)
          newcell = node->cells[NW];
        else
          newcell = node->cells[SW];
      }
    }
  }
  else
    GP_PANIC("bad dir!");
  
  return newcell;
}

/* follow vector direction until we hit an edge of cell */
/* and return x,y of that point and a node and direction to search
   to find nearest node on that edge */
void get_new_x_y(double x, double y, double xv, double yv, Gcell *cell,
    double *ret_x, double *ret_y, G_nodes **pnode, char *dir)
{
  double east_y, west_y, north_x, south_x;
  double x0, y0, x1, y1;

  x0 = cell->x0;
  x1 = cell->x1;
  y0 = cell->y0;
  y1 = cell->y1;

  if (xv != 0) {
    east_y = edge_coord(x1, x, y, xv, yv);
    west_y = edge_coord(x0, x, y, xv, yv);
  }
  else
    /* pick bogus val greater than x1 so these edges aren't checked.*/
    west_y = east_y = y1+ fabs(y1);

  if (yv != 0) {
    north_x = edge_coord(y1, y, x, yv, xv);
    south_x = edge_coord(y0, y, x, yv, xv);
  }
  else
    /* pick bogus val greater than x1 so these edges aren't checked.*/
    north_x = south_x = x1 = fabs(x1);

  if (xv >= 0 && east_y <= y1 && east_y >= y0) {
    /* we hit east edge */
    *pnode = cell->bndry.nodes[SE];
    *dir = N;
    *ret_x = x1;
    *ret_y = east_y;
  }
  else if (xv <= 0 && west_y <= y1 && west_y >= y0) {
    /* west edge */
    *pnode = cell->bndry.nodes[SW];
    *dir = N;
    *ret_x = x0;
    *ret_y = west_y;
  }
  else if (yv >= 0 && north_x <= x1 && north_x >= x0) {
    /* north edge */
    *pnode = cell->bndry.nodes[NW];
    *dir = E;
    *ret_x = north_x;
    *ret_y = y1;
  }
  else if (yv <=0 && south_x <= x1 && south_x >= x0) {
    /* south edge */
    *pnode = cell->bndry.nodes[SW];
    *dir = E;
    *ret_x = south_x;
    *ret_y = y0;
  }
  else
    GP_PANIC("get_new_x_y: Vector never hits edge of box?");

}


/* given initial coordinates, and a unit vector direction, find
   the other final coordinate (y_final) */
double edge_coord(double x_final, double x_initial, double y_initial,
    double xv, double yv)
{
  return (x_final - x_initial)/xv*yv + y_initial;
}


/* divide cell until containing cell is less than max_size */
/* (unless it is a hole, then do nothing) */
Gcell *cut_cell(double x, double y, double max_x, double max_y, Gcell *cell,
    Nonuni_gp *gp)
/* double max_x, max_y; max size of cell this point can be in */
{
  double x_width, y_width, x_rat, y_rat;

  if (c_is_hole(cell))
    return cell;

  x_width = cell->x1 - cell->x0;
  y_width = cell->y1 - cell->y0;

  while(max_x < x_width || max_y < y_width) {

    x_rat = x_width/max_x;
    y_rat = y_width/max_y;

    if (x_rat > y_rat)
      break_cell(cell, EW, gp);
    else
      break_cell(cell, NS, gp);

    cell = get_containing_cell(x,y,cell);

    x_width = cell->x1 - cell->x0;
    y_width = cell->y1 - cell->y0;
  }

  return cell;
}

void break_cell(Gcell *cell, char dir, Nonuni_gp *gp)
{
  
  /* allocate and set up two children */
  make_two_kids(cell, dir, gp);

  /* set the children coordinates */
  set_cell_coords(cell, cell->x0, cell->y0, cell->x1, cell->y1);

  /* update node information of cells */
  update_bi_nodes(cell, gp);
}

void update_bi_nodes(Gcell *cell, Nonuni_gp *gp)
{
  Bi *two_kids;
  Gcell *c1, *c2;
  char type;
  G_nodes *node;
  
  two_kids = (Bi *)cell->children;
  c1 = two_kids->child1;
  c2 = two_kids->child2;
  type = two_kids->type;

  /* i'm no longer going to maintain the edge information */
  gp->is_edge_corrupted = TRUE;
  
  if (type == NS) {
    c1->bndry.nodes[NE] = cell->bndry.nodes[NE];
    c1->bndry.nodes[NW] = cell->bndry.nodes[NW];
    
    c2->bndry.nodes[SE] = cell->bndry.nodes[SE];
    c2->bndry.nodes[SW] = cell->bndry.nodes[SW];

    /* for node on west edge, start with SW node and move north looking 
       for a mid node */

    /* find mid node if it exists on west edge (from node[SW] to node[NW]) */
    find_or_make_node(cell, c1, c2, SW, NW, c1->x0, c1->y0, gp);

    /* look for mid node on east edge (from SE to NE) */
    find_or_make_node(cell, c1, c2, SE, NE, c1->x1, c1->y0, gp);

    /* set adjaceny connecting two "new" nodes */
    c1->bndry.nodes[SW]->adjacent[E] = c1->bndry.nodes[SE];
    c1->bndry.nodes[SE]->adjacent[W] = c1->bndry.nodes[SW];
  }
/* SRW - fixed
  else if (type = EW) {
*/
  else if (type == EW) {
    c1->bndry.nodes[NE] = cell->bndry.nodes[NE];
    c1->bndry.nodes[SE] = cell->bndry.nodes[SE];
    
    c2->bndry.nodes[NW] = cell->bndry.nodes[NW];
    c2->bndry.nodes[SW] = cell->bndry.nodes[SW];

    /* find mid node on north boundary (from NW to NE) */
    find_or_make_node(cell, c1, c2, NW, NE, c1->x0, c1->y1, gp);

    /* find mid node on south boundary (from SW to SE) */
    find_or_make_node(cell, c1, c2, SW, SE, c1->x0, c1->y0, gp);

    /* set adjaceny connecting two "new" nodes */
    c1->bndry.nodes[NW]->adjacent[S] = c1->bndry.nodes[SW];
    c1->bndry.nodes[SW]->adjacent[N] = c1->bndry.nodes[NW];
  }

  fix_node_cell_ptrs(c1);
  fix_node_cell_ptrs(c2);
  clear_edge_ptrs(cell);
}

/* set new parent edge pointers to NULL so that we know they are corrupted*/
void clear_edge_ptrs(Gcell *cell)
{
  int i;

  for(i = 0; i < NUMEDGES; i++)
    cell->bndry.edges[i] = NULL;
}
     
/* update the nodes on the boundary of this cell to point in the right place */
void fix_node_cell_ptrs(Gcell *cell)
{
  G_nodes *node, **corners;
  int i;

  corners = cell->bndry.nodes;

  for(i = 0; i < NUM_N_CELLS; i++)
    corners[i]->cells[opposite_dir(i)] = cell;

  /* set nodes on edges to point to this cell */
  set_edge_nodes(corners[SW], corners[NW], N, NE, SE, cell);
  set_edge_nodes(corners[SE], corners[NE], N, NW, SW, cell);
  set_edge_nodes(corners[SW], corners[SE], E, NE, NW, cell);
  set_edge_nodes(corners[NW], corners[NE], E, SE, SW, cell);

}

/* set cell ptrs for nodes between node1 and node2 */
void set_edge_nodes(G_nodes *node1, G_nodes *node2, char adj_dir,
    char dir1, char dir2, Gcell *cell)
{
  G_nodes *node;
  
  for(node = node1->adjacent[adj_dir]; node != node2; 
      node = node->adjacent[adj_dir]) {
    node->cells[dir1] = cell;
    node->cells[dir2] = cell;
  }
}
     
void find_or_make_node(Gcell *cell, Gcell *kid1, Gcell *kid2,
    char node_start, char node_end, double x, double y, Nonuni_gp *gp)
{
  char adj_dir;
  G_nodes *node;

  /* determine adjacency direction to go based on which nodes chosen */
  if (node_start == SW && node_end == NW)
    adj_dir = N;
  else if (node_start == SW && node_end == SE)
    adj_dir = E;
  else if (node_start == NW && node_end == NE)
    adj_dir = E;
  else if (node_start == SE && node_end == NE)
    adj_dir = N;
  else
    GP_PANIC("find_or_make_node: Unknown combo of start and end!");

  /* find mid node if it exists on west edge */
  node = find_mid_node(cell->bndry.nodes[node_start], adj_dir, 
		       cell->bndry.nodes[node_end]);
  if (node != NULL) {
    kid1->bndry.nodes[node_start] = node;
    kid2->bndry.nodes[node_end] = node;
  }
  else {
    /* no mid node.  make a new one.  opposite_dir stuff doesn't matter */
    node = make_new_node(x, y, opposite_dir(node_start), 
			 kid1, ++gp->num_nodes);

    /* coincidentally, adjacent cell dirs are same as node_start, node_end.*/
    /* Other cell dirs updated one level up in fix_node_cell_ptrs() */
    node->cells[node_start] = node->cells[node_end]
      = cell->bndry.nodes[node_start]->cells[node_end];

    kid1->bndry.nodes[node_start] = kid2->bndry.nodes[node_end] = node;
    gp->nodelist = add_to_gnodelist(node, gp->nodelist);
    fix_adjacency(cell->bndry.nodes[node_start], adj_dir, 
		  node, cell->bndry.nodes[node_end]);
  }
}

void fix_adjacency(G_nodes *begin_node, char adj_dir, G_nodes *node,
    G_nodes *end_node)
{
  node->adjacent[adj_dir] = end_node;
  node->adjacent[opposite_dir(adj_dir)] = begin_node;

  begin_node->adjacent[adj_dir] = node;
  end_node->adjacent[opposite_dir(adj_dir)] = node;

}  


G_nodes *find_mid_node(G_nodes *begin_node, char adj_dir, G_nodes *end_node)
{
  G_nodes *node;
  double x,y;
  
  if (begin_node->adjacent[adj_dir] == end_node)
    /* nothing in between */
    return NULL;

  /* mid node coords */
  x = (begin_node->x + end_node->x)/2.0;
  y = (begin_node->y + end_node->y)/2.0;

  /* now search edge for mid node. (it must exist or else the neighbor
     does not have a power-of-two number of children which is an error) */
  
  /* I'll use a brute force search of the edge for the appropriate
     coordinates.  Something more elegant would be nice but
     a pain to implement */

  node = begin_node->adjacent[adj_dir];

  while( !nearzero(node->x - x) || !nearzero(node->y - y) ) 
    node = node->adjacent[adj_dir];

  return node;
}


/* breaks a cell into two.  Turns a leaf into a BI (two kids) */
void make_two_kids(Gcell *cell, char dir, Nonuni_gp *gp)
{
  Gcell *cell1, *cell2;
  Bi *bi_kids;

  cell1 = new_Gcell(NOCLEAR);
  cell2 = new_Gcell(NOCLEAR);

  cell->children_type = BI;
  cell->children = (void *)gp_malloc( sizeof(Bi));
  bi_kids = (Bi *)cell->children;
  bi_kids->type = dir;
  bi_kids->child1 = cell1;
  bi_kids->child2 = cell2;

  cell1->index = ++gp->num_cells;
  cell2->index = ++gp->num_cells;
  cell1->children = cell2->children = NULL;
  cell1->children_type = cell2->children_type = NONE;
  cell1->parent = cell2->parent = cell;

  gp->num_leaves++;
}



/* Insure that all cells withing a rectangle are have max_x and max_y 
   dimensions at most.  Will affect area outside rectangle also */
/* Calls walk_along_line in rectangle regiongs to accomplish it's task */
/* vals[0]-[2] center of rectangel, vals[3] width in x-dir of rect. 
   val[4] width in y-direction.  val[5] = max_x - maximum width of cells
   in x direction.  val[6] - max width in y-dir */

void contact_rect(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double x0, y0, z0;
  double xl, yl, xr, yr;
  Gcell *contain;
  double rect_x_width, rect_y_width, max_cell_x, max_cell_y;

  if (contactp->numvals != 7) 
    contact_error("Exactly 7 values required for a rect contact.","",
		  contactp);

  /* beginning point */
  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &x0, &y0, &z0);

  rect_x_width = vals[3]*units;
  rect_y_width = vals[4]*units;
  max_cell_x = vals[5]*units;
  max_cell_y = vals[6]*units;

  xl = x0 - rect_x_width/2.0;
  xr = x0 + rect_x_width/2.0;

  yl = y0 - rect_y_width/2.0;
  yr = y0 + rect_y_width/2.0;

  if (xl < gp->root_cell->x0) {
    contact_warning("Left edge of rectangle not in plane.  Truncating.",
		    contactp);
    xl = gp->root_cell->x0;
  }

  if (xr > gp->root_cell->x1) {
    contact_warning("Right edge of rectangle not in plane.  Truncating.",
		    contactp);
    xr = gp->root_cell->x1;
  }

  if (yl < gp->root_cell->y0) {
    contact_warning("Bottom edge of rectangle not in plane.  Truncating.",
		    contactp);
    yl = gp->root_cell->y0;
  }

  if (yr > gp->root_cell->y1) {
    contact_warning("Top edge of rectangle not in plane.  Truncating.",
		    contactp);
    yr = gp->root_cell->y1;
  }
  
  x0 = (xl + xr)/2.0;
  y0 = (yl + yr)/2.0;
  rect_x_width = xr - xl;
  rect_y_width = yr - yl;

  if (vals[3] <= 0.0 || vals[4] <= 0.0 || vals[5] <= 0.0 || vals[6] <= 0.0)
    contact_error("Widths in x-dir and y-dir must be greater than 0.","",
		  contactp);

  /* go inside rect and cut up cells */
  cut_inside_rect(x0, y0, z0, rect_x_width, rect_y_width, max_cell_x, 
		   max_cell_y, gp);

}

/* walk around inside of this rectangle at (x,y,z) and with dimensions
   (rect_x_width, rect_y_width), and cut cells to be no bigger than 
   (max_cell_x, max_cell_y) */
/* it calls walk_along_line() for lines spaced appropriately to guarantee
   no big rectangle squeezes through */

void cut_inside_rect(double x, double y, double z, double rect_x_width,
    double rect_y_width, double max_cell_x, double max_cell_y, Nonuni_gp *gp)
{
  double x0, y0, x1, y1; /* coordinates of SW and NE corners of rectangle */
  double xt, yt;

  /*
  fprintf(stderr, "cutting %lg, %lg, %lg, wid=%lg,%lg cellw=%lg,%lg\n",
	  x, y, z, rect_x_width, rect_y_width, max_cell_x, max_cell_y);
          */

  x0 = x - rect_x_width/2.0;
  x1 = x + rect_x_width/2.0;

  y0 = y - rect_y_width/2.0;
  y1 = y + rect_y_width/2.0;

  /* walk in y direction for each value of x */
  for(xt = x0; xt < x1; xt += max_cell_x)
    walk_along_line(xt, y0, z, xt, y1, z, max_cell_x, max_cell_y, gp);

  walk_along_line(x1, y0, z, x1, y1, z, max_cell_x, max_cell_y, gp);

  /* walk in x direction for each value of y */
  for(yt = y0; yt < y1; yt += max_cell_y)
    walk_along_line(x0, yt, z, x1, yt, z, max_cell_x, max_cell_y, gp);

  walk_along_line(x0, y1, z, x1, y1, z, max_cell_x, max_cell_y, gp);
  
}

/* does same as contact_rect, except outside of rectangle,
   segments will gradually get bigger given by some decay function
   until they reach vals[7] = stop_x, vals[8] = stop_y.
   stop_x and stop_y < 0 mean never stop.
*/
void contact_decay_rect(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double x0, y0, z0;
  double xl, yl, xr, yr;
  Gcell *contain;
  double rect_x_width, rect_y_width, max_cell_x, max_cell_y;
  double stop_x, stop_y;
  
  if (contactp->numvals != 9) 
    contact_error("Exactly 9 values required for a rect contact.","",
		  contactp);

  /* beginning point */
  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &x0, &y0, &z0);

  rect_x_width = vals[3]*units;
  rect_y_width = vals[4]*units;
  max_cell_x = vals[5]*units;
  max_cell_y = vals[6]*units;
  stop_x = vals[7]*units;
  stop_y = vals[8]*units;

  xl = x0 - rect_x_width/2.0;
  xr = x0 + rect_x_width/2.0;

  yl = y0 - rect_y_width/2.0;
  yr = y0 + rect_y_width/2.0;

  if (xl < gp->root_cell->x0) {
    contact_warning("Left edge of decay rectangle not in plane.  Truncating.",
		    contactp);
    xl = gp->root_cell->x0;
  }

  if (xr > gp->root_cell->x1) {
    contact_warning("Right edge of decay rectangle not in plane.  Truncating.",
		    contactp);
    xr = gp->root_cell->x1;
  }

  if (yl < gp->root_cell->y0) {
    contact_warning("Bottom edge of decay rectangle not in plane.  Truncating.",
		    contactp);
    yl = gp->root_cell->y0;
  }

  if (yr > gp->root_cell->y1) {
    contact_warning("Top edge of decay rectangle not in plane.  Truncating.",
		    contactp);
    yr = gp->root_cell->y1;
  }
  
  x0 = (xl + xr)/2.0;
  y0 = (yl + yr)/2.0;
  rect_x_width = xr - xl;
  rect_y_width = yr - yl;

  if (vals[3] <= 0.0 || vals[4] <= 0.0 || vals[5] <= 0.0 || vals[6] <= 0.0)
    contact_error("Widths in x-dir and y-dir must be greater than 0.","",
		  contactp);

  if (rect_x_width <= max_cell_x) {
    contact_warning("rectangle x width <= max cell width. Will set max_cell = 0.99 * rect_width\nNote: Almost no gradual decay in x direction will result",
		    contactp);
    max_cell_x = 0.99*rect_x_width;
  }

  if (rect_y_width < max_cell_y) {
    contact_warning("rectangle y width <= max cell width. Will set max_cell = 0.99 * rect_width\nNote: Almost no gradual decay in y direction will result",
		    contactp);
    max_cell_y = 0.99*rect_y_width;
  }

  /* go inside rect and cut up cells */
  do_decay_rect(x0, y0, z0, rect_x_width, rect_y_width, max_cell_x, 
		   max_cell_y, stop_x, stop_y, gp);
     
}

void do_decay_rect(double x, double y, double z, double rect_x_width,
    double rect_y_width, double max_cell_x, double max_cell_y, double stop_x,
    double stop_y, Nonuni_gp *gp)
{

  double x0, y0;
  double temp_width_x, temp_width_y;
  double xl, yl, xr, yr;
  double x_min, x_max, y_min, y_max;

  x_min = gp->root_cell->x0;
  x_max = gp->root_cell->x1;
  y_min = gp->root_cell->y0;
  y_max = gp->root_cell->y1;

  /* change corners of box to fit inside plane */
  limit_box(x,y,rect_x_width, rect_y_width, x_min, x_max, y_min, y_max,
	    &xl, &yl, &xr, &yr);
  
  x0 = (xl + xr)/2.0;
  y0 = (yl + yr)/2.0;
  temp_width_x = xr - xl;
  temp_width_y = yr - yl;
  
  while((max_cell_x < stop_x || max_cell_y < stop_y || stop_x<0 || stop_y < 0)
	&& (xl > x_min || xr < x_max || yl > y_min || yr < y_max)) {

    /* cut up rectangle */
    cut_inside_rect(x0, y0, z, temp_width_x, temp_width_y, max_cell_x, 
		    max_cell_y, gp);

    /* figure out next rectangle and widths */
    compute_new_widths(&rect_x_width, &rect_y_width, &max_cell_x, &max_cell_y);

    /* will infinite loop (with stack overflow if cut_cell not recursive) */
    if (rect_x_width < 0 || rect_y_width < 0 || max_cell_x < 0 
	|| max_cell_y < 0)
      GP_PANIC("do_decay_rect: width less than 0!");

    /* change corners of box to fit inside plane */
    limit_box(x,y,rect_x_width, rect_y_width, x_min, x_max, y_min, y_max,
	      &xl, &yl, &xr, &yr);
    
    x0 = (xl + xr)/2.0;
    y0 = (yl + yr)/2.0;
    temp_width_x = xr - xl;
    temp_width_y = yr - yl;

  }

  /* do one last cut */

  cut_inside_rect(x0, y0, z, temp_width_x, temp_width_y, max_cell_x, 
		  max_cell_y, gp);

}

void limit_box(double x, double y, double x_wid, double y_wid, double xl_min,
    double xr_max, double yl_min, double yr_max, double *xl, double *yl,
    double *xr, double *yr)
{
  /* truncate box if it falls outside of limits */
  *xl = x - x_wid/2.0;
  *xr = x + x_wid/2.0;
  
  *yl = y - y_wid/2.0;
  *yr = y + y_wid/2.0;
  
  if (*xl < xl_min)
    *xl = xl_min;
  if (*xr > xr_max)
    *xr = xr_max;
  if (*yl < yl_min)
    *yl = yl_min;
  if (*yr > yr_max)
    *yr = yr_max;
}


/* computes next bigger rectangle and larger max_widths for cells */
/* It does the following:
   Assume there is a source at the center which
   produces a current density that drops like 1/r from the center.
   Then this chooses new dimensions so that the new segments
   farther away will have the same net current as the one just inside
   the old rectangle.
   i.e. i match the condition  
      int_{rect - cell}^{rect} {1/r dr} = int_{rect}^{new_rect} {1/r dr}

   much easier to do write ratio = cell/rect and then cell = rect/ratio

   then new_rect = rect*(1/(1-ratio)) and new_cell = rect(ratio/(1-ratio))
*/

void compute_new_widths(double *x_rect, double *y_rect,
    double *x_cell, double *y_cell)
/* double *x_rect, *y_rect; dimensions of rectangle */
/* double *x_cell, *y_cell; max width of cells in rectangle */
{
  double x_ratio, y_ratio;

  x_ratio = (*x_cell)/(*x_rect);
  y_ratio = (*y_cell)/(*y_rect);

/*
  printf("before: xr=%lg yr=%lg x_cell=%lg y_cell=%lg x_rect=%lg y_rect=%lg\n",
	 x_ratio, y_ratio, *x_cell, *y_cell, *x_rect, *y_rect);
         */

  *x_cell = *x_rect*(x_ratio/(1 - x_ratio));
  *x_rect = *x_rect*(1.0 / (1.0 - x_ratio));

  *y_cell = *y_rect*(y_ratio/(1 - y_ratio));
  *y_rect = *y_rect*(1.0 / (1.0 - y_ratio));

/*
  printf("after: xr=%lg yr=%lg x_cell=%lg y_cell=%lg x_rect=%lg y_rect=%lg\n",
	 x_ratio, y_ratio, *x_cell, *y_cell, *x_rect, *y_rect);
         */
}

		
/* equivalences (shorts) all the nodes inside rect.     
   center of rectangle: vals[0]-[2], xwidth = [3], ywidth=[4]
*/
void contact_equiv_rect(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double x0, y0, z0;
  double xl, yl, xr, yr;
  Gcell *contain;
  double rect_x_width, rect_y_width;
  
  if (contactp->numvals != 5) 
    contact_error("Exactly 5 values required for a rect contact.","",
		  contactp);

  /* beginning point */
  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &x0, &y0, &z0);

  contain = get_containing_cell(x0, y0, gp->root_cell);
  
  if (!is_in_cell(x0,y0,contain)) 
    contact_error("Contact equiv point not in plane!","",contactp);

  rect_x_width = vals[3]*units;
  rect_y_width = vals[4]*units;

  make_equiv_rect(x0, y0, z0, rect_x_width, rect_y_width, gp, contactp->name);

}

/* equiv all the nodes together on the exterior of the given rectangle */
/* picks the closest node outside the rect along the perimeter. */
/* We will assume the elements on the interior don't matter */
void make_equiv_rect(double x0, double y0, double z0, double x_width,
    double y_width, Nonuni_gp *gp, char *name) 
/* char *name;  node name */
{
  double xl, xr, yb, yt;
  Gcell *cell;
  GROUNDPLANE *grndp = gp->grndp;
  SYS *indsys = gp->grndp->indsys;
  double x,y,z;
  char templine[200];
  G_nodes *center, *gnode;
  double xg,yg,zg;
  NODES *cnode, *node, *realcnode;
  int i, isalready_done;
  Llist *inside_nodes, *endl, *list;

  xl = x0 - x_width/2.0;
  xr = x0 + x_width/2.0;

  yb = y0 - y_width/2.0;
  yt = y0 + y_width/2.0;

  /* put rectangle totally inside plane */
  if (xl < gp->root_cell->x0)
    xl = gp->root_cell->x0;
  if (yb < gp->root_cell->y0)
    yb = gp->root_cell->y0;
  if (xr > gp->root_cell->x1)
    xr = gp->root_cell->x1;
  if (yt > gp->root_cell->y1)
    yt = gp->root_cell->y1;

  inside_nodes = get_nodes_inside_rect(xl, yb, xr, yt, gp->root_cell, &endl);

  if (inside_nodes == NULL) {
    fprintf(stderr, "Error: No nodes inside equiv_rect %s. Rectangle must be made larger,\n     or plane must be refined more at that point\n",name);
    return;
  }

  /* use first node to equiv all the others to (not really the center) */
  center = (G_nodes *)inside_nodes->ptr;
  
  get_global_coords(center->x, center->y, z0, gp, &xg, &yg, &zg);

#if 1==0 
  the old equiv_rect 

  cell = get_containing_cell(x0, y0, gp->root_cell);

  /* the center node */
  center = find_nearest_edge_node(x0,y0,cell);

  get_global_coords(x0, y0, z0, gp, &xg, &yg, &zg);
#endif

  cnode = get_nonuni_node_from_list(center, grndp->usernodes);
  if (cnode != NULL)
    fprintf(stderr, "Warning: equiv_rect %s's center node appears to be part of another contact\n", name);
  else {
    sprintf(templine, "n%s_%s_pseudo_center", grndp->name, name);
    cnode = make_fastH_node(center, grndp, gp, xg, yg, zg, templine);
  }

  /* add user given name to pseudo name list so user can refer to this */
  append_pnlist(create_pn(name, cnode), indsys);

  realcnode = getrealnode(cnode);

  /*short together (equiv) (force to be same potential) all the nodes in list*/
  for(list = inside_nodes->next; list != NULL; list = list->next) {
    gnode = (G_nodes *)list->ptr;
        node = get_nonuni_node_from_list(gnode, grndp->usernodes);
    if (node == NULL) {
      get_global_coords(gnode->x, gnode->y, z0, gp, &x, &y, &z);
      sprintf(templine, "n%s_%s_pseudo_%lg_%lg", grndp->name, name, gnode->x,
              gnode->y);
      node = make_fastH_node(gnode, grndp, gp, x, y, z, templine);
      isalready_done = FALSE;
    }      
    else {
      /* already exists */
      if (getrealnode(node) != realcnode) {
        isalready_done = FALSE;
        fprintf(stderr,"Warning: contact appears to overlap another contact\n");
      }
      else
        isalready_done = TRUE;

    }

    if (isalready_done == FALSE)
      make_equiv(node, realcnode);
  }

  free_Llist(inside_nodes);

#if 1==0
  the old way

  /* now walk along edge, equiv'ing all the nodes to cnode */
  
  /* start in the southwest and go north, walking on edge */
  walk_along_edge(xl, yb, xr, yt, gp, WEST, cnode, name, z0);
  walk_along_edge(xl, yb, xr, yt, gp, NORTH, cnode, name, z0);
  walk_along_edge(xl, yb, xr, yt, gp, EAST, cnode, name, z0);
  walk_along_edge(xl, yb, xr, yt, gp, SOUTH, cnode, name, z0);
#endif

}

/* walk along edge of "equiv_rect", shorting together (equivalencing) 
   all the nodes we hit. */
/* *** This is no longer called *** */

void walk_along_edge(double xl, double yb, double xr, double yt,
    Nonuni_gp *gp, char which_edge, NODES *cnode, char *name, double z0)
/* NODES *cnode;  the node to short to */
/* char *name;  name of equiv_rect */
{
  double x,y;
  Gcell *cell;
  char start_node, end_node, travel_dir; /*travel dir should be N or E always*/
  char scan_dir, equal_dir;
  G_nodes *anode;

  char is_EW; /* are we doing and EAST or WEST edge */

  if (which_edge == WEST) {
    /* for the west edge, start at SW and move N to NW */
    start_node = SW;
    end_node = NW;
    travel_dir = N;
    
    /* at NW node, go E on north edge until we are at node closest to
       the line x=xl */
    scan_dir = E;

    /* if, when scanning along edge, we end up at x=xl, which cell to choose*/
    equal_dir = NE;

    is_EW = TRUE;

    x = xl;
    y = yb;
  }

  if (which_edge == EAST) {
    start_node = SE;
    end_node = NE;
    travel_dir = N;
    
    /* at NE node, go W on north edge until we are at node closest to
       the line x=xr */
    scan_dir = W;

    /* if, when scanning along edge, we end up at x=xr, which cell to choose*/
    equal_dir = NW;

    is_EW = TRUE;

    x = xr;
    y = yb;
  }
  
  if (which_edge == NORTH) {
    start_node = NW;
    end_node = NE;
    travel_dir = E;
    
    scan_dir = S;

    /* if, when scanning along edge, we end up at y=yt, which cell to choose*/
    equal_dir = SE;

    is_EW = FALSE;

    x = xl;
    y = yt;
  }

  if (which_edge == SOUTH) {
    start_node = SW;
    end_node = SE;
    travel_dir = E;
    
    scan_dir = N;

    /* if, when scanning along edge, we end up at y=yb, which cell to choose*/
    equal_dir = NE;

    is_EW = FALSE;

    x = xl;
    y = yb;
  }

  /* find cell for the starting point */
  cell = get_containing_cell(x,y,gp->root_cell);
  
  /* do the first cell */
  equiv_nodes_on_edge(cell->bndry.nodes[start_node], travel_dir, 
                      cell->bndry.nodes[end_node], cnode, gp, name, z0);
  if (is_EW)
    y = cell->bndry.nodes[end_node]->y;
  else 
    x = cell->bndry.nodes[end_node]->x;

  /* while we are still on the edge, keep equiv'ing nodes together */
  while( ((which_edge == WEST || which_edge == EAST) && (y < yt)) 
         || ((which_edge == NORTH || which_edge == SOUTH) && (x < xr)) ) {

    /* go along edge looking for a closer node in case the adjacent
       cell is smaller */
    anode = scan_edge(x, y, cell->bndry.nodes[end_node], scan_dir);

    if (is_EW == TRUE) {
      if (anode->x == x)
        cell = anode->cells[equal_dir];
      else if (anode->x > x)
        /* nearest node is to the E, so adjacent cell is the W */
        /* i choose NW because I assume travel_dir is to the N */
        cell = anode->cells[NW];
      else
        cell = anode->cells[NE];
    }

    if ( !(is_EW == TRUE)) {
      if (anode->y == y)
        cell = anode->cells[equal_dir];
      else if (anode->y > y)
        /* nearest node is to the N, so adjacent cell is the S */
        /* We choose SE because i assume travel_dir is to the E */
        cell = anode->cells[SE];
      else
        cell = anode->cells[NE];
    }

    /* equiv nodes on next cell's boundary also. */
    equiv_nodes_on_edge(cell->bndry.nodes[start_node], travel_dir, 
                        cell->bndry.nodes[end_node], cnode, gp, name, z0);
    if (is_EW)
      y = cell->bndry.nodes[end_node]->y;
    else 
      x = cell->bndry.nodes[end_node]->x;
    
  }

}

/* short cnode to all the nonuni_nodes from node1 to node2 going in
   the "dir" adjacency direction */
/* *** this is no longer called *** */
void equiv_nodes_on_edge(G_nodes *node1, char dir, G_nodes *node2,
    NODES *cnode, Nonuni_gp *gp, char *name, double z0)
{
  NODES *node, *realcnode;
  G_nodes *gnode;
  double x,y,z;
  GROUNDPLANE *grndp = gp->grndp;
  char isalready_done;
  char keepgoing; 
  char templine[200];

  realcnode = getrealnode(cnode);

  keepgoing = TRUE;

  for(gnode = node1; keepgoing == TRUE; gnode = gnode->adjacent[dir]) {
    node = get_nonuni_node_from_list(gnode, grndp->usernodes);
    if (node == NULL) {
      get_global_coords(gnode->x, gnode->y, z0, gp, &x, &y, &z);
      sprintf(templine, "n%s_%s_pseudo_%lg_%lg", grndp->name, name, gnode->x,
              gnode->y);
      node = make_fastH_node(gnode, grndp, gp, x, y, z, templine);
      isalready_done = FALSE;
    }      
    else {
      /* already exists */
      if (getrealnode(node) != realcnode) {
        isalready_done = FALSE;
        fprintf(stderr,"Warning: contact appears to overlap another contact\n");
      }
      else
        isalready_done = TRUE;

    }

    if (isalready_done == FALSE)
      make_equiv(node, realcnode);

    if (gnode == node2)
      keepgoing = FALSE;
  }
        
}

/* makes a new fasthenry NODES structure for this Nonuni_gp node.
   A call to "get_nonuni_node_from_list" must have been shown
   to return NULL before calling this function */
NODES *make_fastH_node(G_nodes *node, GROUNDPLANE *grndp, Nonuni_gp *gp,
    double xg, double yg, double zg, char *pseudoname)
{
  NODES *cnode;

  cnode = make_new_node_with_nonuni(node, pseudoname, 
                                    grndp->indsys->num_nodes,
                                    xg, yg, zg, grndp->indsys, gp);

  /* make a pseudo_seg between other ground plane usernodes */
  /* (to be used by graph searching algs) */
  grndp->fake_seg_list 
    = make_new_fake_segs(cnode, grndp->usernodes, grndp->fake_seg_list);
  grndp->usernodes = add_node_to_list(cnode, grndp->usernodes);
  
  return cnode;
}
  


/* set up an initial 2D grid of cells. */
/*  vals[0] is the number of cells in x, vals[1] is the number in y */
void contact_initial_grid(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double xc, yc, zc;
  Gcell *contain;

  if (contactp->numvals != 2) 
    contact_error("Exactly 2 values required for an initial grid.","",
		  contactp);

  if ( fabs(vals[0] - (int)vals[0]) > EPS 
       || fabs(vals[1] - (int)vals[1]) > EPS )
    contact_error("Both values should be integers: num cells in x, num in y.",
                  "", contactp);

  /* root cell must have no children */
  make_initial_grid(gp, (int)vals[0], (int)vals[1]);
}

void make_initial_grid(Nonuni_gp *gp, int x_cells, int y_cells)
{
  Gcell *root = gp->root_cell;
  Grid_2d *grid;

  if (c_get_children_type(root) != NONE) {
    fprintf(stderr, "Error: attempting to cut an initial grid that has already be divided.\n x_cells=%d, y_cells=%d\n",x_cells,y_cells);
    exit(1);
  }

  make_grid_kids(root, x_cells, y_cells, gp);

  /* set the children coordinates */
  set_cell_coords(root, root->x0, root->y0, root->x1, root->y1);

  /* update node info */
  update_grid_nodes(root, gp);

}

/* make node information for new grid.  Assumes we are gridding the root cell
   and checks to insure it */
void update_grid_nodes(Gcell *cell, Nonuni_gp *gp)
{
  int i, j;
  Grid_2d *grid;
  G_nodes *node, *nodeNW, *nodeNE, *nodeSW, *nodeSE;
  Gcell ***kids, *onekid;
  int x_c, y_c;

  if (c_get_children_type(cell) != GRID_2D) {
    GP_PANIC("update_grid_nodes: not a grid_2d cell!");
  }
  else {
    grid = (Grid_2d *)cell->children;
  }

  /* i'm no longer going to maintain the edge information */
  gp->is_edge_corrupted = TRUE;
  
  x_c = grid->x_cells;
  y_c = grid->y_cells;
  kids = grid->kids;


  /* make nodes for the grid */
  /* assumes [0][0] is in top left (NW corner).  This was a silly
     choice since the SW corner would be better since its the origin
     and set_grid_coords() would be easier to write */
  for(i = 0; i < y_c; i++)
    for(j = 0; j < x_c; j++) {

      onekid = kids[i][j];

      /* do north nodes */
      if (i == 0) {
        /* make nodes for north side */
        /* do NW */
        if (j == 0) {
          nodeNW = cell->bndry.nodes[NW];
          nodeNW->adjacent[S] = nodeNW->adjacent[E] = NULL;
        }
        else 
          nodeNW = kids[i][j-1]->bndry.nodes[NE];

        /* do NE */
        if (j == x_c - 1) {
          nodeNE = cell->bndry.nodes[NE];
          nodeNE->adjacent[S] = nodeNE->adjacent[W] = NULL;
        }
        else
          nodeNE = add_new_node(onekid->x1, onekid->y1, SW, onekid,
                               ++gp->num_nodes, gp);
      }
      else {
        /* nodes already exist from cells to the north */
        nodeNW = kids[i - 1][j]->bndry.nodes[SW];
        nodeNE = kids[i - 1][j]->bndry.nodes[SE];
      }

      /* do south nodes */

      /* do SW */
      if (j == 0) {
        if (i == y_c - 1) {
          nodeSW = cell->bndry.nodes[SW];
          nodeSW->adjacent[N] = nodeSW->adjacent[E] = NULL;
        }
        else
          nodeSW = add_new_node(onekid->x0, onekid->y0, NE, onekid,
                                 ++gp->num_nodes, gp);
      }
      else 
        nodeSW = kids[i][j - 1]->bndry.nodes[SE];
        
      /* do SE */
      if (j == x_c - 1 && i == y_c - 1) {
        nodeSE = cell->bndry.nodes[SE];
        nodeSE->adjacent[N] = nodeSE->adjacent[W] = NULL;
      }
      else
        nodeSE = add_new_node(onekid->x1, onekid->y0, NW, onekid,
                               ++gp->num_nodes, gp);

        
      /* fix node info and cell info */
      set_node_and_cell_info(nodeNW, NW, onekid);
      set_node_and_cell_info(nodeNE, NE, onekid);
      set_node_and_cell_info(nodeSW, SW, onekid);
      set_node_and_cell_info(nodeSE, SE, onekid);

      set_cell_node_adjacency(onekid);
    }

  clear_edge_ptrs(cell);
}

/* have this node and cell point at each other */
void set_node_and_cell_info(G_nodes *node, int dir, Gcell *cell)
{
  cell->bndry.nodes[dir] = node;
  node->cells[opposite_dir(dir)] = cell;
}

/* have nodes on boundary of this cell point at each other */
void set_cell_node_adjacency(Gcell *cell)
{
  G_nodes **nodes = cell->bndry.nodes;

  point_at_each_other(nodes[NE], S, nodes[SE]);
  point_at_each_other(nodes[NW], E, nodes[NE]);
  point_at_each_other(nodes[SE], W, nodes[SW]);
  point_at_each_other(nodes[SW], N, nodes[NW]);

}

void point_at_each_other(G_nodes *node1, int dir, G_nodes *node2)
{

  int opp_dir;

  opp_dir = opposite_dir(dir);

  if (node1->adjacent[dir] == NULL)
    node1->adjacent[dir] = node2;
  else if (node1->adjacent[dir] != node2)
    GP_PANIC("point_at_each_other: node1->adjacent != NULL and != node2");
    
  if (node2->adjacent[opp_dir] == NULL)
    node2->adjacent[opp_dir] = node1;
  else if (node2->adjacent[opp_dir] != node1)
    GP_PANIC("point_at_each_other: node2->adjacent != NULL and != node1");
}


/* calls make_new_node() to make a new node, then adds this one
   to the linked list */
G_nodes *add_new_node(double x, double y, int cell_dir, Gcell *cell,
    int index, Nonuni_gp *gp)
{
  G_nodes *node;
  
  node = make_new_node(x, y, cell_dir, cell, index);
  gp->nodelist = add_to_gnodelist(node, gp->nodelist);

  return node;
}

void make_grid_kids(Gcell *parent, int x_cells, int y_cells, Nonuni_gp *gp)
{
  Grid_2d *grid;
  Gcell *cell;
  int i,j;
  
  grid = (Grid_2d *)gp_malloc(sizeof(Grid_2d));

  grid->x_cells = x_cells;
  grid->y_cells = y_cells;

  grid->kids = (Gcell ***)gp_malloc( y_cells*sizeof(Gcell **));

  /* make 2d array of pointers to cells */
  for(i = 0; i < y_cells; i++) {
    grid->kids[i] = (Gcell **)gp_malloc( x_cells*sizeof(Gcell *));
    for(j = 0; j < x_cells; j++) {
      cell = new_Gcell(NOCLEAR);
      grid->kids[i][j] = cell;
      cell->index = ++gp->num_cells;
      cell->children = NULL;
      cell->children_type = NONE;
      cell->parent = parent;
    }
  }

  parent->children_type = GRID_2D;
  parent->children = (void *)grid;
  
  gp->num_leaves += x_cells*y_cells - 1; /* parent is nolonger a leaf */

}


/* set up an initial 2D grid of cells and then put holes in it. */
/*  vals[0] is the number of cells in x, vals[1] is the number in y */
/* it puts holes every other row and every other column */
void contact_initial_mesh_grid(ContactList *contactp, Nonuni_gp *gp,
    double relx, double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double xc, yc, zc;
  Gcell *contain;

  if (contactp->numvals != 2) 
    contact_error("Exactly 2 values required for an initial grid.","",
		  contactp);

  if ( fabs(vals[0] - (int)vals[0]) > EPS 
       || fabs(vals[1] - (int)vals[1]) > EPS )
    contact_error("Both values should be integers: num cells in x, num in y.",
                  "", contactp);

  /* root cell must have no children */
  make_initial_grid(gp, (int)vals[0], (int)vals[1]);

  poke_holes(gp);
}

/* make a meshed groundplane out of an initially refined plane. */
/*  This pokes holes in every other row, every other column.  You
    could really do any pattern */
void poke_holes(Nonuni_gp *gp)
{
  int i, j;
  Grid_2d *grid;

  if (c_get_children_type(gp->root_cell) != GRID_2D)
    GP_PANIC("poke_holes: ground plane does not have an initial grid!");

  grid = (Grid_2d *)gp->root_cell->children;

  /* put holes in every other row */
  for(i = 1; i < grid->y_cells; i+=2)
    /* put holes in every other column */
    for(j = 1; j < grid->x_cells; j+=2)
      grid->kids[i][j]->ishole = TRUE;
}

  
/* this turns a connection contact into two regular contacts:
 A decay_rect and an equiv_rect.  It takes 6 args. vals[0]-vals[4] are the 
 same as equiv_rect: (x,y,z) are [0]-[2] and the width of rect in x and y are
 [3]and [4].  vals[5] is the ratio of rectangle dimensions to largest allowed
 size of rectangle inside.  ie, arg 5 and 6 to decay_rect are arg 3 / arg 5
 and arg 4 / arg 5 */
ContactList *make_contact_connection(ContactList *head, char *line,
    double units, double relx, double rely, double relz, int *skip)
{
  char name[MAXNAME];
  char nname[80];
  int skip1, skip2, skip3, skip4, skip25;
  ContactList *con_decay, *con_equiv;
  char *linep;
  double val, ratio;
  int i, numvals;

#if 1==0
  name[MAXNAME-1] = '\0';
  sscanf(line, "%19s%n",name,&skip1);
  line+=skip1;

  if (strcmp("contact",name) != 0) {
    fprintf(stderr, "Internal Error:  %s is not the word 'contact'\n",name);
    exit(1);
  }

  name[MAXNAME-1] = '\0';
  sscanf(line, "%19s%n",name,&skip2);
  line+=skip2;

  if (name[MAXNAME-1] != '\0')
    printf("Warning: contact function '%s' truncated to 19 chars\n",name);

  if (strcmp(name,"connection") != 0)
    GP_PANIC("make_contact_connection: How did you get here?");
#endif

  con_decay = (ContactList *)Gmalloc(sizeof(ContactList));
  con_equiv = (ContactList *)Gmalloc(sizeof(ContactList));
  
  con_decay->func = (char *)Gmalloc((strlen("decay_rect") + 1)*sizeof(char));
  con_equiv->func = (char *)Gmalloc((strlen("equiv_rect") + 1)*sizeof(char));
  strcpy(con_decay->func, "decay_rect");
  strcpy(con_equiv->func, "equiv_rect");
  con_decay->units = units;
  con_decay->relx = relx;
  con_decay->rely = rely;
  con_decay->relz = relz;

  con_equiv->units = units;
  con_equiv->relx = relx;
  con_equiv->rely = rely;
  con_equiv->relz = relz;

  /* use name for equiv_rect */
  skip25 = 0;
  sscanf(line,"%s%n",nname,&skip25);
  line += skip25;
  con_equiv->name = (char *)Gmalloc((strlen(nname) + 1)*sizeof(char));
  strcpy(con_equiv->name, nname);
    
  skip3 = 0;
  while(isspace(*line) && *line != '\0') {
    line++;
    skip3++;
  }

  if (*line != '(')
    contact_error2("Values for contact must start with '('",line,"connection"); 
  
  /* let's count how many values we have first */
  linep = line;
  numvals = 0;

  while(*linep != ')' && !eos(*linep)) {
    linep++;
    while(isspace(*linep) || is_one_of(*linep, "1234567890.e+-"))
      linep++;
    numvals++;
  }
  if (*linep != ')')
    contact_error2("Contact values did not end with ')'", line, "connection");

  if (numvals != 6)
    contact_error2("connection type contact must have 6 values",line,"connection");
  
  con_equiv->numvals = 5;
  con_equiv->vals = (double *)Gmalloc(5*sizeof(double));
  con_decay->numvals = 9;
  con_decay->vals = (double *)Gmalloc(9*sizeof(double));

  skip3 += linep - line + 1;

  /* now let's read them in */
  linep = line + 1;  /* skip the ( */
  /*vals = contactp->vals;*/

  /* read the x,y,z and widths of the rectangle */
  i = 0;
  while(*linep != ')') {
    if (sscanf(linep, "%lf%n",&val,&skip4) != 1)
      contact_error2("Couldn't read value starting at:",linep,"connection");
    else {
      linep += skip4;
      if (i != 5) {
        con_equiv->vals[i] = val;
        con_decay->vals[i] = val;
      }
      else
        ratio = val;

      linep += skipspace(linep);
      if (*linep == ',')
	linep++;
      i++;
    }
  }

  *skip += skip25 + skip3;

  /* set the size of the largest cells inside the rect */
  con_decay->vals[5] = con_decay->vals[3]/ratio;
  con_decay->vals[6] = con_decay->vals[4]/ratio;

  /* the size of a rectangle to stop cutting */
  con_decay->vals[7] = -1.0; /* never stop */
  con_decay->vals[8] = -1.0; /* never stop */

  con_decay->next = head;
  con_decay->done = FALSE;

  con_equiv->next = con_decay;
  con_equiv->done = FALSE;

  fprintf(stdout, "Translating contact connection %s into:\n",nname);
  fprintf(stdout, "contact decay_rect (");
  for(i = 0; i < con_decay->numvals; i++)
    fprintf(stdout,"%lg%s",con_decay->vals[i], 
            (i == con_decay->numvals - 1 ? ")\n" : ","));

  fprintf(stdout, "contact equiv_rect %s (",nname);
  for(i = 0; i < con_equiv->numvals; i++)
    fprintf(stdout,"%lg%s",con_equiv->vals[i], 
            (i == con_equiv->numvals - 1 ? ")\n" : ","));
  
  return con_equiv;

}

/* discretize underneath a trace.  Insure that directly underneath the
   region, the cells are at most as wide as half the width of the trace.  And
   adjacent to those right underneath, no bigger than width of trace .
   should we do more?*/
/* vals[0]-[2] first point of line, vals[3]-[5] other end point, 
   vals[6] = width of trace, vals[7] is a computation reducing factor. width
   and height of cells underneath trace at 45 degrees will be vals[7]*vals[6]
   instead of just vals[6] (so set to 1 for no reduction) The reduction
   comes at the price of accuracy of course */
void contact_trace(ContactList *contactp, Nonuni_gp *gp, double relx,
    double rely, double relz, double units)
{
  double *vals = contactp->vals;
  int i1,j1;
  double x0, y0, z0, x1, y1, z1;
  Gcell *contain;

  if (contactp->numvals != 8) 
    contact_error("Exactly 8 values required for a line contact.","",
		  contactp);

  /* beginning point */
  get_nonuni_coords(vals[0]*units + relx, vals[1]*units + rely, 
		    vals[2]*units + relz, gp, &x0, &y0, &z0);

  /* end point */
  get_nonuni_coords(vals[3]*units + relx, vals[4]*units + rely, 
		    vals[5]*units + relz, gp, &x1, &y1, &z1);
  
  contain = get_containing_cell(x0, y0, gp->root_cell);
  
  if (!is_in_cell(x0,y0,contain)) 
    contact_error("First point of contact trace not in plane!","",contactp);

  contain = get_containing_cell(x1, y1, gp->root_cell);
  
  if (!is_in_cell(x1,y1,contain)) 
    contact_error("Second point of contact trace not in plane!","",contactp);

  if (vals[6] <= 0.0 || vals[7] < 1.0)
    contact_error("Width of trace must be greater than zero, and scale factor greater than 1.","",
		  contactp);

  /* walk along line and cut up cells.  don't scale the factor (vals[7]) */
  do_trace(x0, y0, z0, x1, y1, z1, vals[6]*units, vals[7], gp);

  
}

/* discretize underneath a trace.  Insure that directly underneath the
   region, the cells are at most as wide as half the width of the trace.  And
   adjacent to those right underneath, no bigger than width of trace .
   should we do more?*/
void do_trace(double x0, double y0, double z0, double x1, double y1,
    double z1, double width, double factor, Nonuni_gp *gp)
{
  double max_x, max_y;
  double tan_th;
  double len, reallength;
  double half_width = width/2.0;
  double x_orth, y_orth, scale;

  reallength = sqrt( SQUARE(x1-x0) + SQUARE(y1-y0));

  /* the orthogonal direction */
  x_orth = -(y1 - y0)/reallength;
  y_orth = (x1 - x0)/reallength;

  if ( (x1 - x0) == 0.0 ) {
    max_x = half_width;
    max_y = fabs(y1 - y0);
  }
  else if ( (y1 - y0) == 0.0 ) {
    max_y = half_width;
    max_x = fabs(x1 - x0);
  }
  else {
    tan_th = fabs(y1 - y0)/fabs(x1 - x0);
    len = half_width / tan_th;

    /* choose a continuously varying length from infinity 
       at th = 0 or th = 90 to len = width/2.0 at th = 45 */
    /* all scaled by a factor dependent on the angle */
    if (tan_th <= 1) {
      /* theta < 45 */
      scale = pow(factor,tan_th); /* map [0,1]->[1,factor] */
      max_y = half_width*scale;
      max_x = len*scale;
    }
    else {
      scale = pow(1.0/factor,tan_th); /* map [0,1]->[1,factor] */
      max_x = half_width*scale;
      max_y = len*scale;
    }
  }

  max_x = MIN(max_x,reallength);
  max_y = MIN(max_y,reallength);
  
  /* cut up cells directly under trace */
  walk_along_line(x0, y0, z0, x1, y1, z1, max_x, max_y, gp);

  /* this now cuts like decay_rect with ratio = 2 for just a couple steps */

  /* cut up cells on edge of trace */
  walk_along_line(x0 + x_orth*half_width, y0 + y_orth*half_width, z0, 
                  x1 + x_orth*half_width, y1 + y_orth*half_width, z1, 
                  max_x, max_y, gp);
  walk_along_line(x0 - x_orth*half_width, y0 - y_orth*half_width, z0, 
                  x1 - x_orth*half_width, y1 - y_orth*half_width, z1, 
                  max_x, max_y, gp);

  /* cut up cells a bit off edge of trace */
  walk_along_line(x0 + 3*x_orth*half_width, y0 + 3*y_orth*half_width, z0, 
                  x1 + 3*x_orth*half_width, y1 + 3*y_orth*half_width, z1, 
                  max_x*2, max_y*2, gp);
  walk_along_line(x0 - 3*x_orth*half_width, y0 - 3*y_orth*half_width, z0, 
                  x1 - 3*x_orth*half_width, y1 - 3*y_orth*half_width, z1, 
                  3*max_x, 3*max_y, gp);

  
  /* do anymore??? */
}  

/* pick the cell in the direction (xv,yv) from this node */
Gcell *pick_cell_based_on_vec(G_nodes *node, double xv, double yv)
{

  if (xv > 0 && yv > 0)
    return node->cells[NE];
  else if (xv > 0 && yv < 0)
    return node->cells[SE];
  else if (xv < 0 && yv > 0)
    return node->cells[NW];
  else if (xv < 0 && yv < 0)
    return node->cells[SW];
  else if (yv == 0) {
    /*if yv==0 it doesn't matter which we return.  */
    /* but if we are scating the edge of domain, pick one that exists */
    if (xv > 0) {
      if (node->cells[SE] != NULL)
        return node->cells[SE];
      else
        return node->cells[NE];
    }
    else if (xv < 0) {
	if (node->cells[SW] != NULL)
	  return node->cells[SW];
	else
	  return node->cells[NW];
    }
    else
      GP_PANIC("pick_cell_based_on_vec(): xv == yv == 0!")
  }
  else if (xv == 0) {
    /* we've hit a node.  Doesn't matter which we return.  */
    /* but if we are scating the edge of domain, pick one that exists */
    if (yv > 0) {
      if (node->cells[NE] != NULL)
        return node->cells[NE];
      else 
        return node->cells[NW];
    }
    else {
      if (node->cells[SE] != NULL)
        return node->cells[SE];
      else 
        return node->cells[SW];
    }
  }
  else
    GP_PANIC("pick_cell_based_on_vec(): No condition met!");
}


