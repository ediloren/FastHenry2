

#ifndef _GP_H
#define _GP_H

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


#define c_get_x0(cell) (cell->x0)
#define c_get_y0(cell) (cell->y0)
#define c_get_x1(cell) (cell->x1)
#define c_get_y1(cell) (cell->y1)

/* could also check if children_type == NONE */
#define c_is_leaf(cell) (cell->children == NULL)
#define c_get_children_type(cell) (cell->children_type)
#define c_get_bi_type(two) (two->type)
#define c_is_hole(cell) (cell->ishole == TRUE)

#define GP_PANIC(str) { fprintf(stderr,"Internal error in nonuniform plane code: %s\n",str); /* debug_func(); */ exit(1); }

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

/* SRW -- prototypes */

/* contact.c */
ContactList *make_contactlist(ContactList*, char*, double, double, double,
    double, int*);
// void contact_error(char*, char*, ContactList*);
// void contact_error2(char*, char*, char*);
// void contact_warning(char*, ContactList*);
// void regurg_contact(FILE*, ContactList*);
void make_contacts(ContactList*, Nonuni_gp*);
// void contact_point(ContactList*, Nonuni_gp*, double, double, double, double);
// void contact_line(ContactList*, Nonuni_gp*, double, double, double, double);
// void walk_along_line(double, double, double, double, double, double,
//     double, double, Nonuni_gp*);
// Gcell *find_next_cell_along_line(double, double, double, double, Gcell*,
//     double*, double*);
// void get_new_x_y(double, double, double, double, Gcell*, double*, double*,
//     G_nodes**, char*);
// double edge_coord(double, double, double, double, double);
// Gcell *cut_cell(double, double, double, double, Gcell*, Nonuni_gp*);
// void break_cell(Gcell*, char, Nonuni_gp*);
// void update_bi_nodes(Gcell*, Nonuni_gp*);
// void clear_edge_ptrs(Gcell*);
// void fix_node_cell_ptrs(Gcell*);
// void set_edge_nodes(G_nodes*, G_nodes*, char, char, char, Gcell*);
// void find_or_make_node(Gcell*, Gcell*, Gcell*, char, char, double, double,
//     Nonuni_gp*);
// void fix_adjacency(G_nodes*, char, G_nodes*, G_nodes*);
// G_nodes *find_mid_node(G_nodes*, char, G_nodes*);
// void make_two_kids(Gcell*, char, Nonuni_gp*);
// void contact_rect(ContactList*, Nonuni_gp*, double, double, double, double);
// void cut_inside_rect(double, double, double, double, double, double, double,
//     Nonuni_gp*);
// void contact_decay_rect(ContactList*, Nonuni_gp*, double, double, double,
//     double);
// void do_decay_rect(double, double, double, double, double, double, double,
//     double, double, Nonuni_gp*);
// void limit_box(double, double, double, double, double, double, double,
//     double, double*, double*, double*, double*);
// void compute_new_widths(double*, double*, double*, double*);
// void contact_equiv_rect(ContactList*, Nonuni_gp*, double, double, double,
//     double);
// void make_equiv_rect(double, double, double, double, double, Nonuni_gp*,
//     char*); 
// void walk_along_edge(double, double, double, double, Nonuni_gp*, char,
//     NODES*, char*, double);
// void equiv_nodes_on_edge(G_nodes*, char, G_nodes*, NODES*, Nonuni_gp*,
//     char*, double);
// NODES *make_fastH_node(G_nodes*, GROUNDPLANE*, Nonuni_gp*, double, double,
//     double, char*);
// void contact_initial_grid(ContactList*, Nonuni_gp*, double, double, double,
//     double);
// void make_initial_grid(Nonuni_gp*, int, int);
// void update_grid_nodes(Gcell*, Nonuni_gp*);
// void set_node_and_cell_info(G_nodes*, int, Gcell*);
// void set_cell_node_adjacency(Gcell*);
// void point_at_each_other(G_nodes*, int, G_nodes*);
// G_nodes *add_new_node(double, double, int, Gcell*, int, Nonuni_gp*);
// void make_grid_kids(Gcell*, int, int, Nonuni_gp*);
// void contact_initial_mesh_grid(ContactList*, Nonuni_gp*, double, double,
//     double, double);
// void poke_holes(Nonuni_gp*);
// ContactList *make_contact_connection(ContactList*, char*, double, double,
//     double, double, int*);
// void contact_trace(ContactList*, Nonuni_gp*, double, double, double, double);
// void do_trace(double, double, double, double, double, double, double, double,
//     Nonuni_gp*);
// Gcell *pick_cell_based_on_vec(G_nodes*, double, double);

/* find_nonuni_path.c */
SPATH *path_through_nonuni_gp(NODES*, NODES*, GROUNDPLANE*);
// void clear_nonuni_marks(G_nodes*);
// G_nodes *find_nearest_nonuni_node(double, double, double, Nonuni_gp*);
// SPATH *get_a_nonuni_path(G_nodes*, G_nodes*, Nonuni_gp*, Llist*);
// void sort_nonuni_choices(nonuni_choice_list*, int);
// Llist *add_ptr_to_list(void*, Llist*);
// int add_nonuni_choice(nonuni_choice_list*, Llist*, SEGMENT**, G_nodes*,
//     double);
// int is_ptr_in_list(void*, Llist*);
Gcell *get_containing_cell(double, double, Gcell*);
// Gcell *get_containing_grid_cell(double, double, Gcell*);
// Gcell *get_containing_bi_cell(double, double, Gcell*);
int is_in_cell(double, double, Gcell*);
// G_nodes *find_nearest_edge_node(double, double, Gcell*);
G_nodes *scan_edge(double, double, G_nodes*, char);
// double get_node_dist(double, double, G_nodes*);
// double get_dist(double, double);
int make_nonuni_Mlist(GROUNDPLANE*, MELEMENT**);
// void make_children_meshes(Gcell*, MELEMENT**, int*);
// void make_grid_children_meshes(Grid_2d*, MELEMENT**, int*);
// int make_leaf_mesh(Gcell*, MELEMENT**);
// MELEMENT *add_edge_segs_to_list(G_nodes*, G_nodes*, char, int, MELEMENT*);
NODES *get_or_make_nearest_node(char*, int, double, double, double, SYS*,
    Nonuni_gp*, NPATH*);
NODES *make_new_node_with_nonuni(G_nodes*, char*, int, double, double, double,
    SYS*, Nonuni_gp*);
NODES *get_nonuni_node_from_list(G_nodes*, NPATH*);
Llist *get_nodes_inside_rect(double, double, double, double, Gcell*, Llist**);
// Llist *bi_get_nodes_inside_rect(double, double, double, double, Gcell*,
//     Llist**);
// Llist *grid_get_nodes_inside_rect(double, double, double, double, Gcell*,
//     Llist**);
// void get_grid_indices(Gcell*, double, double, int*, int*);
// int intersection(double, double, double, double, double, double, double,
//     double);
// Llist *which_nodes_inside(double, double, double, double, Gcell*, Llist**);
void free_Llist(Llist*);

/* read_tree.c */
int process_plane(GROUNDPLANE*, FILE*, SYS*);
// void set_gp_coord_system(GROUNDPLANE*, Nonuni_gp*);
void get_nonuni_coords(double, double, double, Nonuni_gp*, double*, double*,
    double*);
void get_global_coords(double, double, double, Nonuni_gp*, double*, double*,
    double*);
// void get_global_vec(double, double, double, Nonuni_gp*, double*, double*,
//     double*);
// int readTree(FILE*, Nonuni_gp*);
void set_cell_coords(Gcell*, double, double, double, double);
// void set_bi_coords(Bi*, double, double, double, double);
// void set_grid_coords(Grid_2d*, double, double, double, double);
// void process_tree(Nonuni_gp*);
// void resolve_nodes(Gcell*, Info*);
// void make_nodes(Gcell*, Info*);
G_nodes *add_to_gnodelist(G_nodes*, G_nodes*);
// void resolve_bi_children(Gcell*, Info*);
// G_edges *make_one_edge(Gcell*, Gcell*, char);
// G_edges *make_two_edge(Gcell*, Gcell*, Gcell*, char);
// Gcell *new_Gcells(int, int);
Gcell *new_Gcell(int);
// void init_Gcell(Gcell*);
// G_nodes *new_Gnode(int);
G_nodes *make_new_node(double, double, int, Gcell*, int);
void *gp_malloc(int);
// void Combine_edges(Gcell*, char, Gcell*, char);
// void combine_node_info(Gcell*, char, Gcell*, char);
// void give_cell_adjaceny(Gcell*, char, Gcell*, char, G_nodes**, G_nodes**);
// void combine_nodes(Gcell*, char, G_nodes*, G_nodes*);
// void replace_node(G_nodes*, G_nodes*);
// void kill_node(G_nodes*);
// void delete_first_node(Nonuni_gp*);
// void delete_dead_nodes(Nonuni_gp*);
// void remove_and_free(G_nodes*);
// void free_g_node(G_nodes*);
// void determine_adjaceny(G_nodes*);
// G_nodes *get_adjacent_node(G_nodes*, Gcell*, char, char, Gcell*, char, char);
// G_nodes *get_other_gnode(Gcell*, char, char);
// void compute_z_fils(Nonuni_gp*);
// void generate_segs(Nonuni_gp*, SYS*);
// void get_width_and_shift(char, G_nodes*, Gcell*, Gcell*, double*, double*);
// void get_x_cell_vals(Gcell*, G_nodes*, Gcell*, double*, double*);
// void get_y_cell_vals(Gcell*, G_nodes*, Gcell*, double*, double*);
// void make_segs(char, G_nodes*, G_nodes*, double, double, double, Nonuni_gp*,
//     SYS*);
// SEGMENT *make_one_seg(double, double, double, double, double, double, double,
//     double, double, double, double, Nonuni_gp*, SYS*, NODES**, NODES**);
// void draw_one_seg(char, double, double, double, double, double, double, double,
//     double, double, double, double, double, double, double, int, Nonuni_gp*);
// void print_cell_and_kids(Gcell*);
// void fprint_cell_and_kids(Gcell*, FILE*);
// void fprint_bi_kids(Bi*, FILE*);
// void dump_cell(Gcell*, FILE*);
// void print_bi_addresses(Bi*, FILE*);
// void print_node_list(G_nodes*);
// void fprint_node_list(G_nodes*, FILE*);
// void dump_node(G_nodes*,FILE*);
// void dump_leaf_cells_to_file(Gcell*, char*);
// void dump_leaf_cells(Gcell*, FILE*);
// void dump_grid_leaf_cells(Grid_2d*, FILE*);
// void print_leaf_cell(Gcell*, FILE*);
void dump_nonuni_plane_currents(Nonuni_gp*, CX*, FILE*);
double get_perimeter(G_nodes*);

#endif /* _GP_H */
