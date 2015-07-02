/*
Copyright (c) 1990 Massachusetts Institute of Technology, Cambridge, MA.
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
