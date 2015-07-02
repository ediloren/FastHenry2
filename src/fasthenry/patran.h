/* # ***** sort to /src/header
   # ***** */
/***************************************************************************

  This is the header file for patran.c.  
  
  Written by Songmin Kim, July 24, 1990.

***************************************************************************/
/* Refer to Chapter 29, Neutral System, of PATRAN manual for explanation.
   The variable names are identical to those that appear in the manual. */

typedef struct node { 
  int ID;
  double coord[3];
} NODE;

typedef struct element { 
  int ID, shape, num_nodes;
  int corner[4];
} ELEMENT;

typedef struct grid {
  int ID, *equiv_ID, number_equiv_grids;
  double coord[3];
  struct grid *next, *prev;
} GRID;

typedef struct cfeg {
  int ID, NELS, LPH, LPH_ID, LSHAPE, NODES, ICONF, NDIM;
  struct cfeg *next, *prev;
  int *element_list;
} CFEG;

typedef struct patch {
  int ID, corner[4], conductor_ID;
  struct patch *next, *prev;
} PATCH;

typedef struct sm_patch {
  int ID, conductor_ID;
  struct sm_patch *next;
} SM_PATCH;

/* intermediate name struct; used for compatability with patran i/f */
typedef struct name {
  char *name;
  SM_PATCH *patch_list;
  struct name *next;
} NAME;

/* used to build linked list of conductor names */
struct Name {
  char *name;
  struct Name *next;
  struct Name *alias_list;
};
typedef struct Name Name;

/* used to make linked lists of iteration or conductor #s */
struct ITER {
  int iter;
  struct ITER *next;
};
typedef struct ITER ITER;
