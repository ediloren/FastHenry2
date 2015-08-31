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

#define FHVERSION "3.0wr"
#define FHDATE "29Sep96, mod 082514"

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
#if SUPERCON == ON
   double lambda;             /* London penetration depth for superconductor */
   double r1;                 /* resistive part of 1/sigma */
   double r2;                 /* reactive part of 1/sigma  */
#endif
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
#if SUPERCON == ON
double lambda;                /* London penetration depth for superconductor */
#endif
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

/* SRW -- prototypes */

/* BreakupSeg.c */
void DivideSegs(double, charge*, SYS*, int);
// void BreakupSeg(SEGMENT*, double, charge*, SYS*);
// void alterFils(SEGMENT*, NODES*, double, double, double);

/* Prec_cost.c */
double OneCubeCost(cube*****, int, int, int, int, int, double*);
double ratio_of_divided_segs(double, charge*, SYS*);
int is_gp_charge(charge*);
// void add_to_counts(cube*, int, int*****, int*****);

/* Precond.c */
// void indPrecond(ssystem*, SYS*, double);
// void multPrecond(PRE_ELEMENT**, CX*, CX*, int);
// MELEMENT *getnext(MELEMENT*, int*);
// void cx_invert_dup(CX**, int, DUPS*);
// void mark_dup_mesh(MELEMENT**, int*, int, DUPS*, int*);
// void dumpPrecond(PRE_ELEMENT**, int, char*);
// void indPrecond_direct(ssystem*, SYS*, double);

/* SetupComputePsi.c */
void SetupComputePsi(CX*, ssystem*, CX*, int, charge*, double, double*, SYS*);
// void realmatCXvec(CX*, double**, CX*, int);
// void fixEvalDirect(charge**, int, int*, charge**, int, double**);

/* SetupMulti.c */
ssystem *SetupMulti(charge*, SYS*);

/* addgroundplane.c */
// double square(double);
// double lengthof(int, int, double*, double*, double*);
// double lengthof2(double*, double, double, double);
// double innerproduct(int, int, int, int, double*, double*, double*);
// double findangle(int, int, int, int, double*, double*, double*);
// void order(int*, int*);
// void bubblesort(int*, int);
int checkmiddlepoint(double*, double*, double*, int, int, int);
int checkplaneformula(double*, double*, double*, double, double, double,
    int, int, int);
double findsegmentwidth(double*, double*, double*, int, int, int, int);
void doincrement(double, double, double, double, double, double, int,
    double*, double*, double*);
void dounitvector(double, double, double, double, double, double,
    double*, double*, double*);
// void fillgrids(GROUNDPLANE*);
void make_nodelist(NODELIST*, char*, double, double, double);
// double find_coordinate(GROUNDPLANE*, double, double, double, int);
// void findrefnodes(GROUNDPLANE*, NODES*, NODES*, int*, int*, int*, int*);
// SPATH *old_path_through_gp(NODES*, NODES*, GROUNDPLANE*);
// MELEMENT **old_makeMlist(GROUNDPLANE*, int*, double**, int);
SPATH *path_through_gp(NODES*, NODES*, GROUNDPLANE*);
// SPATH *get_a_path(NODES*, NODES*, GROUNDPLANE*, NPATH*, int, int);
// void sort_choices(choice_list*, int);
void clear_marks(SYS*);
// void increment_usage(SEGMENT*);
void dump_mesh_coords(SYS*);
void dump_ascii_mesh_coords(SYS*);
// void **Make_C_array(void*, int, int, int);
// int add_choice(choice_list*, NPATH*, SEGMENT*, NODES*, int, int, int, int, int);
// void clear_plane_marks(GROUNDPLANE*);

/* barnoldi.c */
typedef int (*barnoldi_cb)(double**, ssystem*, double**, charge*, SYS*, int,
    int, int);
int ArnoldiROM(double**, double**, double**, char**, int, int, int, int,
    barnoldi_cb, SYS*, ssystem*, charge*);
// int qr(double**, double**, double**, int, int, int);
// int qr_P(double**, double**, double**, double**, int, int, int, char*);
int dumpROM(FILE*, double**, double**, double**, double**, int, int, int);
void dumpROMequiv_circuit(FILE*, double**, double**, double**, double**,
    int, int, int, char*, char*, SYS*);
// int dumpROMbin(FILE*, double**, double**, double**, double**, int, int, int);
int createMRMt(char**, SYS*);
// int createMRMtinvMLMt(double***, SYS*, char*);
int realComputePsi(double**, ssystem*, double**, charge*, SYS*, int, int, int);
int realMatVect(double**, ssystem*, double**, charge*, SYS*, int, int, int);
// int printRowCol(double**, int, int, int, int);
void formMLMt(SYS*);
// void ZeroMatrix(double**, int, int);

/* bigmeshPre.c */
// void bigmeshPre(ssystem*, SYS*, double);
// cube *getcube(charge*, ssystem*);
// int is_in_list(cube*, cube**, int*);
// void addtolist(cube*, cube**, int*, int*);

/* bigmeshPre_direct.c */
// void bigmesh_direct(ssystem *sys, SYS *indsys, double w)
// int is_in_Precond(PRE_ELEMENT *prelist, int col, PRE_ELEMENT **last)

/* calcp.c */
double calcp(charge*, charge*, double*);
void dumpnums(int, int);
double tilelength(charge*);

/* cx_ludecomp.c */
CX **cx_ludecomp(CX**, int, int);
void cx_lu_solve(CX**, CX*, CX*, int);

/* default_opts.c */
void default_opts(ind_opts*);

/* deg_mutual.c */
enum degen_type find_deg_dims(FILAMENT*);
double compute_for_degenerate(FILAMENT*, FILAMENT*, int, double*, double*,
    enum degen_type, enum degen_type, double);
// void setup_tape_to_tape(FILAMENT*, FILAMENT*, int, double*, double*,
//     enum degen_type, enum degen_type, FILAMENT*, FILAMENT*, double**, double**);
// double do_tape_to_brick(FILAMENT*, FILAMENT*, int, double*, double*,
//     enum degen_type, enum degen_type);

/* dist_betw_fils.c */
double dist_betw_fils(FILAMENT*, FILAMENT*, int*);
// void getD(FILAMENT*, double*);
// void getr(double*, double*, double*, double*, double, double*);
double vdotp(double*, double*);
// double dist_between(double, double, double, double, double, double);
// double min_endpt_sep(FILAMENT*, FILAMENT*);
// double dist_betw_pt_and_fil(FILAMENT*, double*, double*, double, FILAMENT*,
//     double);
// double aspectratio(FILAMENT*);
void fill_Gquad(void);
// void findnfils(FILAMENT*, FILAMENT*, int);
// void gquad_weights(int, double*, double*);

/* fillM.c */
void fillM(SYS*);
// MELEMENT *make_mesh_from_path(SPATH*, int, SYS*);
// int is_next_seg_in_gp(SPATH*, NODES*);
MELEMENT *insert_in_list(MELEMENT*, MELEMENT*);
// NODES *getnode(int, seg_ptr);
// void bad_seg_type(char*, seg_ptr);
MELEMENT *make_melement(int, FILAMENT*, int);
// void add_to_external(PSEUDO_SEG*, int, int, SYS*);
// int_list *make_int_list(int, int);
// int_list *add_to_int_list(int_list*, int_list*);
// int makeMlist(GROUNDPLANE*, MELEMENT**, Minfo*, int);
void fill_b(EXTERNAL*, CX*);
void extractYcol(CX**, CX*, EXTERNAL*, EXTERNAL*);
char *get_a_name(PSEUDO_SEG*);
// NODES *find_first_node(SPATH*);
void makegrids(SYS*, CX*, int, int);

/* findpaths.c */
// PATHLIST *add_to_front(PATHLIST*, PATHLIST*);
char *Gmalloc(int);
// PATHLIST *make_new_path(void);
NODES *getrealnode(NODES*);
NODES *getothernode(NODES*, seg_ptr);
int is_normal(NODES*);
// int is_gp(seg_ptr);
int is_gp_node(NODES*);
int is_node_in_list(NODES*, NPATH*);
int is_orignode_in_list(NODES*, NPATH*);
NPATH *add_node_to_list(NODES*, NPATH*);
// int is_gp_in_list(GROUNDPLANE*, GPLIST*);
// GPLIST *add_to_gplist(GROUNDPLANE*, GPLIST*);
// void free_nodelist(NPATH*);
// void insert_path(SPATH*, PATHLIST*);
// SPATH *copypath(SPATH*);
// SPATH *lastelem(SPATH*);
NODES *get_node_from_name(char*, SYS*);
int equivnodes(char*, SYS*);
PSEUDO_NODE *create_pn(char*, NODES*);
void make_equiv(NODES*, NODES*);
void append_pnlist(PSEUDO_NODE*, SYS*);
// NODES *find_next_external(NODES*);
void add_to_connected_segs(NODES*, SEGMENT*, PSEUDO_SEG*);
void remove_from_connected_segs(NODES*, SEGMENT*, PSEUDO_SEG*);
double mag(double, double, double);
// double magsq(double, double, double);
double dotp(double, double, double, double, double, double);
NODES *find_nearest_gpnode(double, double, double, GROUNDPLANE*, int*, int*);
// int is_real_node(NODES*);
SPATH *add_seg_to_list(seg_ptr, SPATH*);
PSEUDO_SEG *make_pseudo_seg(NODES*, NODES*, char);
SPATH *make_new_fake_segs(NODES*, NPATH*, SPATH*);
EXTERNAL *add_to_external_list(EXTERNAL*, EXTERNAL*);
EXTERNAL *make_external(PSEUDO_SEG*, int, char*, char*, char*);
EXTERNAL *get_external_from_portname(char*, SYS*);
EXTERNAL *get_next_ext(EXTERNAL*);
// NODES *get_next_treeless_node(NODES*);
// TREE *make_new_tree(void);
// TREE *add_tree_to_list(TREE*, TREE*);
// NODES *pop_node(NPATH**);
// void push_node(NODES*, NPATH**);
void make_trees(SYS*);
// SEGLIST *get_next_branch(SEGLIST*);
// int is_marked(seg_ptr);
// void mark_seg(seg_ptr);
// void unmark_seg(seg_ptr);
// void mark_node(NODES*);
// void unmark_node(NODES*);
// int is_node_marked(NODES*);
// PATHLIST *add_path_to_list(SPATH*, PATHLIST*);
// void make_loop(NODES*, NODES*, seg_ptr, TREE*);
int count_tree_meshes(TREE*);
int count_externals(EXTERNAL*);
void find_hole_meshes(SYS*);
// NODES *get_next_gphole_node(NODES*);
// NPATH *find_surrounding(NPATH*);
// void clear_marks_and_level(NPATH*);
// void make_gp_trees(NPATH*, TREE*);
// NPATH *get_next_unexamined_node(NPATH*);
// SEGMENT *get_next_around_hole(NODES*, int*, NPATH*);
// SEGMENT *get_next_gp_seg(NODES*, int*);
// SPATH *make_gp_loop(NODES*, NODES*, seg_ptr);
// void clear_used_segs(PATHLIST*);
// void mark_used_segs(PATHLIST*);
// NODES *find_nearest_node(NODES*, NPATH*);
// double dist_between_nodes(NODES*, NODES*);

/* gmres.c */
typedef CX (*gmres_cb1)(CX*, CX*, int);
typedef void (*gmres_cb2)(CX*, ssystem*, CX*, int, charge*, double, double*,
    SYS*);
// void blockPrecond(CX**, CX*, int, CX*, CX**, int, int*);
void cx_invert(CX**, int);
// void matMult(CX**, CX**, int, int, int, CX*);
// void matMultVect(CX**, CX*, CX*, int);
int gmres(CX**, CX*, CX*, gmres_cb1, gmres_cb2, int, int, double, ssystem*,
     charge*, double, double*, SYS*, int);
// void sub_3(CX*, CX*, CX*, int);
// void set_scalar_mult_1(CX*, double, int);
// void set_scalar_mult_2(CX*, CX*, double, int);
// void add_scalar_mult_2(CX*, CX*, double, int);
// void add_cx_mult_2(CX*, CX*, CX, int);
// void sub_cx_mult_2(CX*, CX*, CX, int);
CX inner(CX*, CX*, int);
// void veryoldmatvec(CX*, CX**, CX*, int);
// void disMat(CX**, int);
// void multvec2(CX**, CX*, CX*, int*, int);
void directmatvec(CX*, ssystem*, CX*, int, charge*, double, double*, SYS*);
// void calc_other_norm(CX*, CX*, int, double, double*, double*, double*, double*);

/* hole.c */
int is_next_word(char*, char*);
int is_hole(NODES*);
HoleList *make_holelist(HoleList*, char*, double, double, double, double, int*);
int skipspace(char*);
int eos(char);
// void hole_error(char*, char*, HoleList*);
int is_one_of(char, char*);
// void delete_node(NODES*);
void make_holes(HoleList*, GROUNDPLANE*);
// void hole_point(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_rect(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_circle(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user1(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user2(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user3(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user4(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user5(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user6(HoleList*, GROUNDPLANE*, double, double, double, double);
// void hole_user7(HoleList*, GROUNDPLANE*, double, double, double, double);

/* induct.c */
charge *assignFil(SEGMENT*, int*, charge*);
double **MatrixAlloc(int, int, int);
// void fillA(SYS*);
// void old_fillM(SYS*);
// void fillZ(SYS*);
#if SUPERCON == ON
// void fillZ_diag(SYS*, double);
// void set_rvals(SYS*, double);
#endif
// double resistance(FILAMENT*, double);
// int countlines(FILE*);
// void savemats(SYS*);
void savecmplx(FILE*, char*, CX**, int, int);
// void savecmplx2(FILE*, char*, CX**, int, int);
// void formMZMt(SYS*);
// void oldformMZMt(SYS*);
char* MattAlloc(int, int);
// void formMtrans(SYS*);
// void compare_meshes(MELEMENT*, MELEMENT*);
// void cx_dumpMat_totextfile(FILE*, CX**, int, int);
// void dumpMat_totextfile(FILE*, double**, int, int);
// void dumpVec_totextfile(FILE*, double*, int);
// void fillMrow(MELEMENT**, int, double*);
// void dump_to_Ycond(FILE*, int, SYS*);
// void saveCarray(FILE*, char*, double**, int, int);
// int nnz_inM(MELEMENT**, int);
// void dump_M_to_text(FILE*, MELEMENT**, int, int);
// void dump_M_to_matlab(FILE*, MELEMENT**, int, int, char*);
// void pick_ground_nodes(SYS*);
// int pick_subset(strlist*, SYS*);
void concat4(char*, char*, char*, char*);

/* joelself.c */
double self(double, double, double);
int edges_parallel(FILAMENT*, FILAMENT*, double*, int*);
void get_wid(FILAMENT*, double*);
void get_height(FILAMENT*, double*, double*);
double exact_mutual(FILAMENT*, FILAMENT*, int, double*, double*,
    enum degen_type, enum degen_type);
// void fill_4(double*, double, double, double);
// double eval_eq(double, double, double, double);
// double log_term(double, double, double, double, double);
// double tan_term(double, double, double, double, double);
int lookup(FILAMENT*, FILAMENT*, int, double*, double*, double*, double*,
    int*, Table***, int*);
// void fill_dims(FILAMENT*, FILAMENT*, double*, double*, double*, int);
// void fill_dims_seg(FILAMENT*, FILAMENT*, double*, double*, double*, int);
// int find_dims(double*, int, Table**, double*, int*, Table***);
void put_in_table(FILAMENT*, FILAMENT*, int, double, double*, int, Table**,
    int);
void init_table(void);
int get_table_mem(void);
void destroy_table(void);
// char *AllocAnEntry(AllocInfo*);
// void DestroyEntries(AllocInfo*);
// int MemoryForEntries(AllocInfo*);
// double brick_to_brick(double, double, double, double, double, double, double,
//     double, double);
// double flat_to_flat_tape(double, double, double, double, double, double,
//     double);
// double eval_eq_tape(double, double, double, double);
// double flat_to_skinny_tape(double, double, double, double, double, double,
//     double);
// double eval_eq_tape2(double, double, double, double);
// double tan_tape(double, double, double, double);
// double tape_to_fil(double, double, double, double, double, double);
// double brick_to_fil(double, double, double, double, double, double, double);

/* mulSetup.c */
ssystem *mulInit(int, int, int, charge*, SYS*, int);
// int placeq(int, ssystem*, charge*, SYS*, int);
// other protos in mulGlobal.h

/* mutual.c */
double mutual(FILAMENT*, FILAMENT*);
// void print_infinity_warning(FILAMENT*, FILAMENT*);
// void findfourfils(FILAMENT*, FILAMENT*);
double selfterm(FILAMENT*);
double mutualfil(FILAMENT*, FILAMENT*);
// double magdiff2(FILAMENT*, int, FILAMENT*, int);
// double mut_rect(double, double);
double dotprod(FILAMENT*, FILAMENT*);
double fourfil(FILAMENT*, FILAMENT*);
// double parallel_fils(FILAMENT*, FILAMENT*, int, double*, double*, double);

/* newPrecond.c */
void choose_and_setup_precond(SYS*);
// void get_selfs(SYS*);
void fill_spPre(ssystem*, SYS*, double);
void create_sparMatrix(SYS*);
// void fill_bySegment(ssystem*, SYS*, double);
// void fill_diagL(ssystem*, SYS*, double);
void fill_diagR(SYS*);
double shift_mutual(FILAMENT*, FILAMENT*, ssystem*);

/* newPrecond_cube.c */
void indPrecond(ssystem*, SYS*, double);
void multPrecond(PRE_ELEMENT**, CX*, CX*, int);
// MELEMENT *getnext(MELEMENT*, int*);
// void cx_invert_dup(CX**, int, DUPS*);
// void mark_dup_mesh(MELEMENT**, int*, int, DUPS*, int*);
// void dumpPrecond(PRE_ELEMENT**, int);
void indPrecond_direct(ssystem*, SYS*, double);

/* parse_command_line.c */
ind_opts *Parse_Command_Line(int, char**);
// Option *gather_opts(int, char**, char*);
// char *Pmalloc(int);
// int is_in_optstring(char, char*, int*);
// int checkarg(int, int, char**);
// void Describe_Usage(char*)
// int read_on_off(char*, int*);
// void add_to_subset_of_columns(char*, ind_opts*);
// void fix_and_print_opts(ind_opts*);

/* readGeom.c */
int readGeom(FILE*, SYS*);
// int dodot(char*, SYS*);
// int changeunits(char*, SYS*);
// int addexternal(char*, SYS*);
// int choosefreqs(char*, SYS*);
// int old_equivnodes(char*, SYS*);
// int dodefault(char*);
// int addnode(char*, SYS*, NODES**, int);
NODES *makenode(char*, int, double, double, double, int, SYS*);
// int addseg(char*, SYS*, int, SEGMENT**);
SEGMENT *makeseg(char*, NODES*, NODES*, double, double, double,
#if SUPERCON == ON
    double,
#endif  
    int, int, double, double, double*, int, int, SYS*);
// int addgroundplane(char*, SYS*, GROUNDPLANE**);
// int nothing(char*);
// char *getaline(FILE*);
// char *plusline(FILE*);
// char *getoneline(FILE*);
// void savealine(char*);
int notblankline(char*);
void tolowercase(char*);
int is_nonuni_gp(GROUNDPLANE*);

/* regurgitate.c */
// typedef void (*regurg_cb)(double, double, double, double*, double*, double*);
void regurgitate(SYS*);
// void do_end_stuff(SYS*);
// void set_translate(double, double, double);
// void translate(double, double, double, double*, double*, double*);
// void reflect_x(double, double, double, double*, double*, double*);
// void reflect_y(double, double, double, double*, double*, double*);
// void reflect_y(double, double, double, double*, double*, double*);
// void do_nothing(double, double, double, double*, double*, double*);
// void spit(SYS*, regurg_cb, char*);

/* savemat_mod.c */
extern int machine_type(void);
// void savemat(FILE*, int, char*, int, int, int, double*, double*);
void savemat_mod(FILE*, int, char*, int, int, int, double*, double*, int, int);

/* uglieralloc.c */
// void fill_it(void*, char, int);
// void ufree(void);
// char *ualloc(unsigned int);
// void ualloc_verify(void);
void uallocEfcy(long);

/* writefastcap.c */
void writefastcap(char*, char*, SYS*);
// void make_between(FILE*, GROUNDPLANE*, double, double, double, double, int,
//     int);
// void unit_cross_product(double, double, double, double, double, double,
//     double*, double*, double*);
// void assign_shades(SYS*);

#endif
