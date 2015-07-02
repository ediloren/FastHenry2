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
