/*
Copyright (c) 1994 Massachusetts Institute of Technology, Cambridge, MA.
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

#ifndef _INDUCT_H
#define _INDUCT_H

/* include file for induct.c */

#include "cmplx.h"
#include "mulGlobal.h"

#define FHVERSION "3.0.1"
#define FHDATE "28May12"

#define AVER_MAT_MAX 110        /* These are two constants used to */
#define MAX_DIV_RAT 0.25        /* decide partitioning level      */

#define MAX_PRE_LEVEL 6         /* maximum allowed level if  */
                                /* auto_refine == OFF */

/*#define MAXITERS 200 */           /* Maximum number of GMRES iters */
#define FILS_PER_MESH 2         /* filaments per mesh for breaking big loops.*/
#define PI 3.14159265358979323846
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define K 0.2235
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13
#define QUICK_EVAL 1          /* do one filament approximations if == 1 */
#define MAXsubfils 6          /* maximum points for Gaussian quadrature */
#define NUM_DIMS 10           /*Number of distinguishing parms. table lookup*/

/* for deg_mutual, how small a dimension must be to be degenerate */
#define DEG_TOL 1.1e-4

/* For node and segment types */
#define NORMAL 0
#define GPTYPE 1
#define PSEUDO 2
#define EXTERNTYPE 3

/* for nodes in gp holes */
#define GPHOLE 4

/* m_info.type types */
#define UNCONSTRAINED 0
#define CONSTRAINED 1


/* The following are used for command line options */

/* defined in mulGlobal.h */
/*
  #define TRUE 1
  #define FALSE 0

  #define ON 1
  #define OFF 0
  #define BOTH 3
*/
#define AUTO -1
#define LUDECOMP 0
#define ITERATIVE 1
#define MULTIPOLE 1
#define DIRECT 0

#define SIMPLE 1
#define REFINED (1 << 1)
#define HIERARCHY (1 << 2)
#define BOTH_FCAP (SIMPLE | REFINED)

/* ground plane visualization */
#define THIN 1
#define THICK 2

/* dump options */
#define DUMP_M 1
#define DUMP_RL (1 << 1)
#define MRL (DUMP_M | DUMP_RL)
#define MZMt (1 << 2)
#define GRIDS (1 << 3)
#define PRE (1 << 4)
#define MESHES (1 << 5)
#define DUMP_A (1 << 6)
#define DUMP_Ls (1 << 7)
#define DUMP_ALL (MRL | MZMt | PRE)

#define MATLAB 1
#define TEXT (1 << 1)
#define BOTH_TYPES (MATLAB | TEXT)

/* precond types */
#define LOC 2       /* local inversion preconditioner */
#define SPARSE 3    /* sparsified L preconditioner */

/* precond subtypes */
#define DIAGL 4     /* using only diagonal of partial inductance matrix L */
#define CUBEL 5     /* using blocks of L defined by multipole cubes */
#define SEGML 6     /* using blocks of L defined by fils in a segments */
#define POSDEF_LOC 7  /* local inv only inverting a single cube (no overlap) */
#define OVERLAP 8   /* standard (default) local inversion. Overlapped precond*/
#define SHELLS 9    /* sparse approx to L using B. Krauter's spherical shells*/

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

/* counters for calls made to each filament mutual inductance routine*/
extern int num_exact_mutual;
extern int num_fourfil;
extern int num_mutualfil;
extern int num_found;
extern int num_perp;

/* machine type for savemat_mod, declated in induct.c */
extern int machine;

/* temporary list of gp node references to be resolved.
   (used in addgroundplane() in file readGeom.c */
typedef struct Nodelist{
  char *name;                      /* name of the node */
  double x, y, z;                     /* coordinates of the node */

  struct Nodelist *next;
} NODELIST;

/* temporary list of holes to be made.
   (used in addgroundplane() in file readGeom.c */
typedef struct _holelist {
  char *func;                   /* function to be called (hole type) */
  double *vals;                 /* array of vals (val1,val2,...) */
  int numvals;                  /* length of array */
  double units;                 /* ground plane info that is not saved */
  double relx, rely, relz;      /*   units and relative x,y,z. */
  struct _holelist *next;       /* Next hole in list */
} HoleList;

/* temporary list of holes to be made.
   (used in addgroundplane() in file readGeom.c */
typedef struct _contactlist {
  char *func;                   /* function to be called (hole type) */
  char *name;                   /* a name as the first arg (sometimes) */
  double *vals;                 /* array of vals (val1,val2,...) */
  int numvals;                  /* length of array */
  double units;                 /* ground plane info that is not saved */
  double relx, rely, relz;      /*   units and relative x,y,z. (not used?)*/
  char done;                    /* have you done this already? */
  struct _contactlist *next;       /* Next hole in list */
} ContactList;

typedef struct sseg_ptr {
  char type;   /* NORMAL means seg is a SEGMENT, PSEUDO means seg is PSEUDO */
  void *segp;   /* a pointer to either a SEGMENT or a PSEUDO_SEG */
} seg_ptr;

typedef struct Filament {
  double x[2], y[2], z[2];  /* endpoints */
  double length, area, width, height;
  double lenvect[3];        /* vector along the length of filament */
  int filnumber;
  struct Segment *segm;
  struct charge *pchg;      /* 'charge' to send to multipole routines */
} FILAMENT;

typedef struct Node {
   char *name;
   int number;
   int index;           /* internal number, for placement in A matrix */
   struct Node *equiv;  /* electically equivalence */
   double x, y, z;
   struct spath *to_end;
   int num_to_end;

   struct seglist *connected_segs;

   int type;                   /* NORMAL or (GPTYPE or GPHOLE) */

   struct Groundplane *gp;      /* CMS 6/7/92 ---- pointer to a groundplane */
   int s1, s2;                  /* indices into groundplane node array */
   struct tree *treeptr;
   char examined;      /* 1 = examined or never to be examined */
   int level;          /* number of nodes away from root of tree */
   struct sseg_ptr pred;      /* predecessor, really just the branch of tree */

   struct g_nodes *gp_node;  /* node of nonuni gp that i really correspond to*/

   struct Node *next;
} NODES;

typedef struct Segment {
   char *name;
   double *widthdir;   /*if width is not || to x-y plane and perpendicular to*/
                       /* the length, then this is 3 element vector in       */
                       /* in the direction of width*/
   int number;         /* an arbitrary number for the segment */
   int type;    /* CMS 8/21/92 -- type of structure the segment is in */
   double length;
   double area;        /* area of cross section */
   double width, height;  /*width and height to cross section */
   int hinc, winc;             /* number of filament divisions in each dir */
   NODES *node[2];                /* nodes at the ends */
   double sigma;              /* conductivity */
   double r_width, r_height;  /*ratio of adjacent fil widths(see assignFil())*/
   int num_fils;               /* hinc*winc */
   FILAMENT *filaments;        /* this segment's filaments */
/*   struct npath *conds;  */    /* linked list of conductors which this seg is in */
   struct pathlist *loops;   /* loops in which this segment is a member */
   int is_deleted;           /* has this segment been used already */

   /*struct g_nodes *gp_node[2];*//* nonuni_gp nodes that are really the ends*/

   /*struct _table *table;*/          /* lookup table for mutual terms */
   struct Segment *next;      /* next segment in list */
 } SEGMENT;

  /* a fake segment to represent the voltage source between terminals */
  /*   or the connection between usernodes in a groundplane */
typedef struct pseudo_seg {
  NODES *node[2];
  char type;  /* Voltage source, or ground plane */
  struct pathlist *loops;   /* loops in which this segment is a member */
  int upper_num_segs;  /* an upper bound on the number of real segments for
                          this pseudo_seg.  For estimating num_meshes before
			  big loops are broken into small ones for new
			  preconditioner 			  */
  int is_deleted;
} PSEUDO_SEG;

typedef struct spath {
  seg_ptr seg;
  struct spath *next;
} SPATH;

typedef struct npath {
  NODES *node;
  struct npath *next;
} NPATH;

typedef struct tree {
  struct pathlist *loops;
  int number_of_loops;
  struct tree *next;
} TREE;

typedef struct _int_list {
  int index;
  int sign;
  struct _int_list *next;
} int_list;

typedef struct external {
  PSEUDO_SEG *source;   /* this branch will represent the 1 volt source */
  int_list *indices;    /* indices of the loops into the M matrix */
  int Yindex;           /* index in the final impedance matrix, Y */
  int col_Yindex;       /* column number, in case -x option is used */
  struct pathlist *loops;
  char *name1, *name2;
  char *portname;
  struct external *next;
} EXTERNAL;

typedef struct melement {  /* an element of the M matrix */
  int filindex;     /* filament number (column of M) */
                /* note: filindex is really a mesh index in indsys->Mtrans */
  FILAMENT *fil;   /* pointer to the filament */
  int sign;         /* 1 = assumed direction of filament i same as mesh i */
  struct melement *mnext; /* next filament in mesh */
} MELEMENT;

typedef struct _minfo {  /* info about a mesh */
  int type;               /* UNCONSTRAINED or CONSTRAINED */
  int mesh_num;           /* mesh number for this mesh */
  /* The following apply only if it is CONSTRAINED or is a constraint for
     another mesh */
  int constraining_mesh;  /* if type==CONSTRAINED then this mesh must have
			     the same mesh current as constraining_mesh */
  int other_mesh;      /* we need one more reference mesh that is constained
			   but not the constraining mesh.  For precond */
  int first;            /* Mesh number for first mesh in Mlist for this group*/
  int num_meshes;             /* number of meshes */
} Minfo;

typedef struct precond_element { /* An element in the preconditioner */
                               /* Each row will be saved as a linked list */
                               /* these */
  int meshcol;     /* column to which this element corresponds */
  CX value;        /* value of this element */
  struct precond_element *next;
} PRE_ELEMENT;

typedef struct _ind_opts ind_opts;  /* to be defined later */

/*----------------------------------------------  CMS 7/7/92 -------------------*/
typedef struct Groundplane{
char *name;                      /* name of the plane */
double x[4], y[4], z[4];         /* corner point coordinates of the plane */
int numesh;                      /* number of KVL meshes */
int seg1;                        /* segments from point 0 to point 1 */
int seg2;                        /* segments from point 1 to point 2 */
int num_nodes1, num_nodes2;      /* number of nodes along each direction */
int row[2], col[2];              /* grid dimensions for current distribution*/
int external;                    /* boolean marker if plane is conductor*/
double **grid1;                  /* segment layout of grid (mid to o1) */
double **grid2;                  /* segment layout of grid (mid to o2) */
SEGMENT ***segs1;                /* segments from mid to o1 */
SEGMENT ***segs2;                /* segments from mid to o2 */
NODES ***pnodes;                 /* grid of plane nodes */
NODES *lastnode;                 /* pointer to the last node in the plane */
NODES *innode, *outnode;         /* pointers to the external plane nodes */
NPATH *usernodes;                /* All user defined nodes for this plane */
SPATH *fake_seg_list;            /* list of pseudo_segs for gp */
double d1;                       /* space between nodes along side 1 */
double d2;                       /* space between nodes along side 2 */
double length1, length2;         /* length of each side */
double ux1, uy1, uz1;            /* unit vector along side 1 */
double ux2, uy2, uz2;            /* unit vector along side 2 */
double unitdiag;                 /* magnitude of (dx1,dy1,dz1)+(dx2,dy2,dz2) */
double sigma;                    /* conductivity of all the segments */
int hinc;                        /* number of filaments to stack in a seg */
double rh;                       /* ratio of adjacent fil heights in a seg */
double thick;                    /* thickness of plane */
double segwid1, segwid2;         /* width of segs in each direction */
HoleList *list_of_holes;         /* list of holes */
ContactList *list_of_contacts;      /* list of contacts (not holes) */
NODELIST *usernode_coords;       /* user referenced coordinates */

/* Nonuniform ground planes */
char *filename;
struct nonuni_gp *nonuni;

struct indsystem *indsys;

struct Groundplane *next;
} GROUNDPLANE;
/*----------------------------------------------------------------------------*/

typedef struct indsystem {
  GROUNDPLANE *planes;         /* CMS 7/1/92  head of plane list */
  GROUNDPLANE *endplane;       /* CMS 7/1/92 end of plane list */

  SEGMENT *segment;            /* head of node list */
  SEGMENT *endseg;             /*end of node list */
  NODES *nodes;                /* head of node list */
  NODES *endnode;              /* end of node list */
  struct pseudo_node *pseudo_nodes;/*lists of names that aren't defined nodes*/
  int num_fils;                /* number of filaments */
  int num_nodes;               /* num of nodes */
  int num_real_nodes;          /* number of real nodes (without equivalences)*/
  int num_segs;                /* number of conductor segments */
  int num_planes;      /* CMS 7/1/92 number of ground planes */
  int num_mesh;                /* number of meshes */
  int num_trees;               /* # of trees (physically separate conductors)*/
  int tree_meshes;             /* big loops from graph */
  int extra_meshes;            /* upper bound on number of meshes from
				  breaking big loops for new precond */
  int num_extern;             /* number of external nodes (also # conductors)*/
  int num_sub_extern;         /* no. of ports (conductors) requested with -x*/
  int *meshsect;               /* array of row indices for beginning of */
                              /* section of mesh matrix for a given conductor*/

  double fmin, fmax;           /* start and end of frequency anaysis */
  double logofstep;           /* log of the frequency step size (sort of)*/
                             /* It is 1/ndec */
  /*double r_height, r_width;*/
  int precond_type;           /* type of precond: Local inv. or sparse */
  int precond_subtype;        /* which type of local inversion or sparse */
  double *diagL;              /* diagonal of the L matrix for SPARSE precond*/
  char *sparMatrix;           /* sparse matrix */

  /* double **A; */    /* incidence matrix  (KCL) */
  MELEMENT **Alist; /* incidence matrix (KCL), only computed if requested */
  double **M;     /* KVL mesh matrix */
  MELEMENT **Mlist; /* KVL list of meshes (to replace M). Each linked list is a row*/
  Minfo *m_info;    /* Information about each mesh */
  MELEMENT **Mtrans;  /* transpose of Mlist.  Each linked list is a row */
  PRE_ELEMENT **Precond;
  struct tree *trees;
  EXTERNAL *externals;    /* for each .external statement. Each element */
          /* in this list will correspond to a column in the final Y matrix */

  double **Z;     /* impedance matrix */
  double *R;      /* Resistance array */
  CX **MtZM;      /* transpose(M)*Z*M */
  CX **FinalY;    /* The conductor admittance matrix */
  double units;   /* factor to convert units to meters */

  double **resids; /* place to store the residuals for each conductor */
  double **resid_real;
  double **resid_imag;
  double *niters;     /* number of iterations for each conductor */

  ind_opts *opts; /* command line options */
  int dont_form_Z;   /* if fmin=0, don't from the L matrix */
  double **Ar, **Br, **Cr, **Dr;  /* reduced order model */

  char *title;  /* title for this input file */
} SYS;

typedef struct dups { /* duplicate mesh information for indPrecond() */
  int sign;          /* zero if not a duplicate, 1 or -1 if a duplicate */
  int dup;           /* mesh number which is the duplicate. */
                     /* sign is 1 if dup is in the same direction,-1 othrwse */
} DUPS;

/***** for graph searching routines ******/
typedef struct pathlist {
  SPATH *path;
  struct pathlist *next;
} PATHLIST;

/* list of segments connect to a node */
typedef struct seglist {
  seg_ptr seg;
  NODES *original;
  struct seglist *next;
} SEGLIST;

typedef struct gplist {
  GROUNDPLANE *gp;
  struct gplist *next;
} GPLIST;

/* for get_a_path() */
typedef struct _choice_list {
  NODES *node;
  SEGMENT *seg;
  int rank;                    /* weight for this choice */
  int s1_mom, s2_mom;  /* 'momentum' for this direction */
} choice_list;

/****** end graph search stuff */

typedef struct pseudo_node {
  NODES *node;
  char *name;
  struct pseudo_node *next;
} PSEUDO_NODE;

/* a list of strings */
typedef struct _strlist {
  char *str;
  struct _strlist *next;
} strlist;

struct _ind_opts {
  int soln_technique;    /* LUDECOMP or ITERATIVE */
  int mat_vect_prod;     /* DIRECT or MULTIPOLE */
  int precond;           /* ON or OFF */
  int order;             /* multipole expansion order */
  int level;             /* multipole partition level. AUTO,0,1...*/
  int makeFastCapFile;   /* Make a fastcap file of the structure.
			    OFF, SIMPLE, REFINED, BOTH */
  int gp_draw;           /* Draw ground planes in fastcap file.
			    OFF - just make outline. ON - Draw all segs */
  int auto_refine;       /* Refine discretization. ON or OFF */
  int init_refine;       /* Initial refinement */
  int dumpMats;          /* Dump matrices. OFF,ON,MRL,MZMt,GRIDS,PRE,MESHES */
  int orderROM;          /* order for reduced order model via block Arnoldi */
  int onlyROM;           /* only do the ROM creation */
  int kind;              /* Kind of dump. MATLAB, TEXT, BOTH */
  double tol;            /* gmres convergence tolerance */
  double abs_tol;        /* gmres absolute tolerance */
  int maxiters;          /* gmres maximum number of iterations */
  int limit;             /* number of filaments for a cube to be exact */
                         /* (used for auto_refine) AUTO,1,2,3... */
  int debug;             /* debug information. ON or OFF */
  strlist *portlist;     /* list of portnames to compute admittance for */
  char *suffix;          /* suffix to append to all output file names*/
  double shell_r0;       /* radius for shell preconditioner */
  int regurgitate;        /* whether or not to spit input file back out */
  char *fname;           /* input filename */
} /* ind_opts */;

/* stuff for mutual terms lookup table */
typedef struct _table {
  double val;
  struct _table *next_val;
  union {
    struct _table *next_dim;
    double *mut_term;
  } u;
} Table;

typedef struct _alloc_list {
  char *ptr;
  struct _alloc_list *next;
} AllocList;

typedef struct _alloc_info {
  int size;
  int blocksize;      /* how many elements to allocate at once */
  int elems_left;
  char *next_elem;
  AllocList *head;
} AllocInfo;

enum degen_type {brick = 0, flat = 1, skinny = 2, too_long = 3, too_short = 4,
		   short_flat = 5, short_skinny = 6, impossible = 7};

PATHLIST *add_to_front();
char *Gmalloc();
PATHLIST *make_new_path();
NODES *getrealnode();
NODES *getothernode();
NPATH *add_node_to_list();
GPLIST *add_to_gplist();
SPATH *copypath();
SPATH *lastelem();
NODES *get_node_from_name();
PSEUDO_NODE *create_pn();
NODES *find_next_external();
double mag();
double magsq();
double dotp();
double vdotp();
double dotprod();
NODES *find_nearest_gpnode();

PATHLIST *add_path_to_list();
SEGLIST *get_next_branch();
int count_tree_meshes();
int count_externals();
SPATH *add_seg_to_list();
PSEUDO_SEG *make_pseudo_seg();
SPATH *make_new_fake_segs();
EXTERNAL *add_to_external_list();
EXTERNAL *make_external();
NODES *get_next_treeless_node();
TREE *make_new_tree();
TREE *add_tree_to_list();
EXTERNAL *get_external_from_portname();
EXTERNAL *get_next_ext();
NODES *find_nearest_node();
double dist_between_nodes();

NPATH *find_surrounding();
NPATH *get_next_unexamined_node();
SEGMENT *get_next_gp_seg();
SPATH *make_gp_loop();
NODES *get_next_gphole_node();
SEGMENT *get_next_around_hole();

/* from fillM.c */
MELEMENT *make_mesh_from_path();
MELEMENT *insert_in_list();
NODES *getnode();
MELEMENT *make_melement();
int_list *make_int_list();
int_list *add_to_int_list();
char *get_a_name();
NODES *find_first_node();

SPATH *path_through_gp();

/* from readGeom */
SEGMENT *makeseg();

/* from parse_command_line.c */
char *Pmalloc();
ind_opts *Parse_Command_Line();

/* from Prec_cost.c */
double OneCubeCost();
double ratio_of_divided_segs();

/* from hole.c */
HoleList *make_holelist();

/* from contact.c */
ContactList *make_contactlist();
NODES *get_or_make_nearest_node();

/* from dist_betw_fils.c */
double dist_betw_fils();
double aspectratio();

/* addgroundplane.c with path search alg */
SPATH *old_path_through_gp();
SPATH *get_a_path();
void **Make_C_array();

/* table lookup stuff */
char *AllocAnEntry();

/* mutual.c and induct.c */
void savemat_mod();
double **MatrixAlloc();
char *MattAlloc();
void fillA();
void fillM();
void fillZ();
double resistance();
double selfterm();
double mutual();
double mutualfil();
double exact_mutual();
/*double magdiff();*/
double magdiff2();
double mut_rect();
extern int gmres();
extern int matvec();
extern int directmatvec();
extern CX inner();
extern int SetupComputePsi();
charge *assignFil();
double fourfil();

double brick_to_brick();
double parallel_fils();
double flat_to_flat_tape();
double flat_to_skinny_tape();
double do_tape_to_brick();
double tape_to_fil();  /* not implemented */
double brick_to_fil(); /* not implemented */
double compute_for_degenerate();

extern int realMatVect();
extern int realComputePsi();

/* SHELLS preconditioner shift */
double shift_mutual();

#endif
