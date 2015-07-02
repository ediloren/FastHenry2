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

/**************************************************************************

  This is the subroutine that handles PATRAN to FastCap interface.
  Patfront returns a pointer to a linked list of charge structures.

  Written by Songmin Kim, July 24, 1990.

**************************************************************************/
#include "mulGlobal.h"
// Enrico
#include <string.h>

#define BIG 35000              /* Size of element and node serach table. */
#define SMALL_NUMBER 0.005     /* See functions if_same_coord() and
				 grid_equiv_check(). */

int type_number, ID, IV, KC, N1, N2, N3, N4, N5;
int number_nodes, number_elements, number_grids=0, number_patches=0;
NODE *list_nodes, *current_node, *node_search_table[BIG];
ELEMENT *list_elements, *current_element, *element_search_table[BIG];
GRID *start_grid;
PATCH *start_patch;
CFEG *start_cfeg;

/* these are now only used for temporary name storage for patran surfaces */
extern NAME *start_name;       	/* name data linked list (1 entry/component) */
extern NAME *current_name;	/* tail pointer to above list */
extern NAME *start_name_this_time;     /* name structs made on this call */
char *title;       	/* overall name given in the title card */
int conductor_count;

/* these flags added to allow multiple calls; used to reset static variables */
int first_grid;			/* note that current_name static is not */
int first_patch;		/*   reset since the name list must */
int first_cfeg;			/*   be preserved as new files are read */

charge *patfront(stream, file_is_patran_type, surf_type, trans_vector,
		 name_list, num_cond, name_suffix)
  FILE *stream;
int *file_is_patran_type, surf_type, *num_cond;
double *trans_vector;
Name **name_list;
char *name_suffix;
{
  int *patch_patch_table, numq=0;
  static char *line = NULL;
  charge *make_charges_all_patches(), *firstq, *quickif();
  double *corner0, *corner1, *corner2, *corner3;

  if(line == NULL) CALLOC(line, BUFSIZ, char, ON, AMSC);

  start_name_this_time = NULL;
  first_grid = first_patch = first_cfeg = TRUE;
  number_grids = number_patches = 0;

  /* figure out input file type and read it in */
  fgets(line, BUFSIZ, stream);
  if(line[0] == '0') {
    *file_is_patran_type = FALSE;
    firstq = quickif(stream, line, title, surf_type, trans_vector,
		     num_cond, name_list, name_suffix);
  }
  else {
    *file_is_patran_type = TRUE;

    input(stream, line, surf_type, trans_vector);

    grid_equiv_check();

    /*********************************************************************
      This section of patfront is for assigning conductor numbers to patches
      depending on their connectivity.                                    */

    if(surf_type == CONDTR || surf_type == BOTH) {

      CALLOC(patch_patch_table, number_patches*number_patches, int, ON, AMSC);

      fill_patch_patch_table(patch_patch_table);

      assign_conductor(patch_patch_table);

    /*********************************************************************/

      assign_names();
    }

    firstq = make_charges_all_patches(name_list, num_cond, surf_type,
				      name_suffix);
  }

  return (firstq);
}


/****************************************************************************

  This part of code is for reading in Patran Neutral files.

****************************************************************************/

input(stream, line, surf_type, trans_vector)
  char *line;
  FILE *stream;
int surf_type;
double *trans_vector;
{
  int END=0;

  /* Reads in the first line from each card, and branch off.  ID, IV, KC,
     N1, N2, N3, N4 and N5 are global variables accessible by subroutines. */

  while (!END) {
    if(line[0] == '2') {	/* if first line */
      sscanf(line,"%d %d %d %d %d %d %d %d %d",
	   &type_number, &ID, &IV, &KC, &N1, &N2, &N3, &N4, &N5);
      line[0] = '0';
    }
    else fscanf(stream,"%d %d %d %d %d %d %d %d %d",
		&type_number, &ID, &IV, &KC, &N1, &N2, &N3, &N4, &N5);


    switch (type_number) {
    case 25:
      file_title(stream);
      break;
    case 26:
      summary_data(stream);
      break;
    case 1:
      node_data(stream, trans_vector);
      break;
    case 2:
      element_data(stream);
      break;
    case 31:
      grid_data(stream, trans_vector);
      break;
    case 33:
      patch_data(stream);
      break;
    case 45:
      CFEG_table(stream);
      break;
    case 21:
      if(surf_type == CONDTR || surf_type == BOTH) name_data(stream);
      else waste_line(KC, stream);
      break;
    case 99:
      END = 1;
      break;
    default:
      waste_line(KC,stream);
      break;
    }
  }
}

/* Simply read in 'num_line' lines from stream and dump. */

waste_line(num_line,stream)
  int num_line;
  FILE *stream;
{
  int c, tmp;
  tmp=num_line+1;
  while (tmp) {
    c=getc(stream);
    if (c=='\n') tmp--;}
}


/* Save the title of the Neutral file. */

file_title(stream)
  FILE *stream;
{
  char line[BUFSIZ], *delcr();

  fgets(line, sizeof(line), stream);
  if(title[0] == '\0') strcpy(title, delcr(line));
}


/* Since the summary card has informations on the number of elements and
   nodes, this function allocates spaces for nodes and elements, and sets up
   the global pointers to these arrays. */

summary_data(stream)
  FILE *stream;
{
  number_nodes = N1; number_elements = N2;

  CALLOC(list_nodes, number_nodes, NODE, ON, AMSC);
  CALLOC(list_elements, number_elements, ELEMENT, ON, AMSC);

  current_node = list_nodes;
  current_element = list_elements;

  waste_line(1,stream);
}


/* Current_node is the global variable that points to the next entry in
   node array, list_nodes, which is preallocated by summary_data function.
   Node_search_table is sorted by node ID to make indexing of a node easier. */

node_data(stream, trans_vector)
  FILE *stream;
double *trans_vector;
{
  double tmp_coord[3];
  int i;

  fscanf(stream,"%lf %lf %lf",tmp_coord,tmp_coord+1,tmp_coord+2);
  waste_line(1,stream);

  for (i=0; i<3; i++)
      current_node->coord[i] = tmp_coord[i] + trans_vector[i];
  node_search_table[ID] = current_node;
  current_node++;
}


/* Current_element is the global variable that points to the next entry in
   element array, list_elements, which is preallocated by summary_data
   function.  Element_search_table is sorted by element ID to make indexing
   of an element easier.  */

element_data(stream)
  FILE *stream;
{
  int num_nodes, corner[4], i, tmp;
  float tmp1;

  current_element->shape = IV;

  if ((IV != 3) && (IV != 4)) waste_line(KC,stream);
  else {
    fscanf(stream,"%d %d %d %d %f %f %f",
           &num_nodes,&tmp,&tmp,&tmp,&tmp1,&tmp1,&tmp1);
    current_element->num_nodes = num_nodes;

    /* IV==3 and 4 imply triangular and quad elements, respectively. */
    if (IV==3) fscanf(stream,"%d %d %d",corner,corner+1,corner+2);
    else fscanf(stream,"%d %d %d %d",corner,corner+1,corner+2,corner+3);

    for (i=0; i<num_nodes; i++) current_element->corner[i] = corner[i];
    element_search_table[ID] = current_element;
    current_element++;

    if (N1) waste_line(1,stream);
  }
}

/* Grid data are linked together by next and prev pointers within GRID
   structure.  Start_grid is the global variable that points to the very
   first GRID structure created.  */

grid_data(stream, trans_vector)
  FILE *stream;
double *trans_vector;
{
  static GRID *prev_grid=0;
  GRID *current_grid;
  double coord[3];
  int i;

  if(first_grid) {
    prev_grid = NULL;
    first_grid = FALSE;
  }

  CALLOC(current_grid, 1, GRID, ON, AMSC);
  if (number_grids==0) start_grid=current_grid;
  current_grid->ID = ID;
  current_grid->prev = prev_grid;
  if (prev_grid) prev_grid->next = current_grid;

  fscanf(stream, "%lf %lf %lf", coord, coord+1, coord+2);
  for (i=0; i<3; i++) current_grid->coord[i] = coord[i] + trans_vector[i];
  prev_grid = current_grid;
  current_grid->next=0;
  number_grids++;
}


/* Patch data are linked together by next and prev pointers within PATCH
   structure.  Start_patch is the global variable that points to the very
   first PATCH structure created.  */

patch_data(stream)
  FILE *stream;
{
  static PATCH *prev_patch=0;
  PATCH *current_patch;
  double tmp;
  int i, corner[4];

  if(first_patch) {
    prev_patch = NULL;
    first_patch = FALSE;
  }

  CALLOC(current_patch, 1, PATCH, ON, AMSC);
  if (number_patches==0) start_patch=current_patch;
  current_patch->ID = ID;
  current_patch->prev = prev_patch;
  if (prev_patch) prev_patch->next = current_patch;

  waste_line(9,stream);
  fscanf(stream, "%f %f %f %d %d %d %d",
	 &tmp, &tmp, &tmp, corner, corner+1, corner+2, corner+3);
  for (i=0; i<4; i++) current_patch->corner[i] = corner[i];
  prev_patch = current_patch;
  current_patch->next=0;
  number_patches++;
}


/* CFEG data are linked together with next and prev pointers within CFEG
   structure.  Start_cfeg is the global variable that points to the very
   first CFEG structure created.  CFEG table has the result from meshing
   a patch. */

CFEG_table(stream)
  FILE *stream;
{
  static CFEG *prev_cfeg=0;
  CFEG *current_cfeg;
  int tmp, NELS, LPH, LPH_ID, LSHAPE, NDIM, NODES, ICONF;
  int i, *element_list, element_num1, element_num2;

  if(first_cfeg) {
    prev_cfeg = NULL;
    first_cfeg = FALSE;
  }

  waste_line(1,stream);
  fscanf(stream,"%d %d %d %d %d %d %d %d",
         &NDIM, &LSHAPE, &NODES, &ICONF, &LPH, &LPH_ID, &tmp, &tmp);

  if (LPH != 3) waste_line(KC-2,stream);
  else {
    CALLOC(current_cfeg, 1, CFEG, ON, AMSC);
    if (!prev_cfeg) start_cfeg=current_cfeg;
    current_cfeg->ID = ID;
    current_cfeg->NELS = IV; NELS = IV;
    current_cfeg->prev = prev_cfeg;
    if (prev_cfeg) prev_cfeg->next = current_cfeg;

    /* This is the list of elements associated with this particular patch. */
    CALLOC(element_list, NELS, int, ON, AMSC);
    current_cfeg->element_list = element_list;

    current_cfeg->LPH = LPH;
    current_cfeg->LPH_ID = LPH_ID;
    current_cfeg->LSHAPE = LSHAPE;
    current_cfeg->NDIM = NDIM;
    current_cfeg->NODES = NODES;
    current_cfeg->ICONF = ICONF;

    if (LSHAPE==3) {                  /* Triangular elements. */
      for (i=1; i<=NELS/2; i++) {
        fscanf(stream, "%d %d %d %d %d %d %d %d %d %d",
          &tmp,&tmp,&tmp,&tmp,&element_num1,&tmp,&tmp,&tmp,&tmp,&element_num2);
        *element_list++ = element_num1;
        *element_list++ = element_num2;
      }
      if (NELS%2) {
        fscanf(stream, "%d %d %d %d %d",&tmp,&tmp,&tmp,&tmp,&element_num1);
        *element_list++ = element_num1;
      }
    }
    else if (LSHAPE==4) {             /* Quad elements. */
      for (i=1; i<=NELS/2; i++) {
        fscanf(stream, "%d %d %d %d %d %d %d %d %d %d",
          &tmp,&tmp,&tmp,&tmp,&element_num1,&tmp,&tmp,&tmp,&tmp,&element_num2);
        *element_list++ = element_num1;
        *element_list++ = element_num2;
      }
      if (NELS%2) {
        fscanf(stream, "%d %d %d %d %d",&tmp,&tmp,&tmp,&tmp,&element_num1);
        *element_list++ = element_num1;
      }
    }
    prev_cfeg = current_cfeg;
  }
}

/*
  reads in name data cards and puts information into NAME struct linked list
  - for every component (each must contain at least 1 conductor surface patch)
    stores the component name and patch ID numbers in that component
  - later Song's patch list is used in assign_names() to set the patch
    conductor ID numbers
  - the output routine looks at the first sm_patch struct associated with
    each NAME struct to determine the number of the corresponding cond name
*/
name_data(stream)
FILE *stream;
{
  int len, iv, i, j, ntype, id, patch_cnt = 0;
  char line[BUFSIZ], *delcr();
  SM_PATCH *current_patch = NULL;

  if(start_name == NULL) {	/* if first time on first patfront() call */
    CALLOC(start_name, 1, NAME, ON, AMSC);
    current_name = start_name_this_time = start_name;
  }
  else{
    CALLOC(current_name->next, 1, NAME, ON, AMSC);
    current_name = current_name->next;
    if(start_name_this_time == NULL) {	/* if 1st time on this patfront call */
      start_name_this_time = current_name;
    }
  }

  /* get conductor name and store */
  fgets(line, sizeof(line), stream); /* eat CR */
  fgets(line, sizeof(line), stream);
  len = strlen(line);
  CALLOC(current_name->name, len+1, char, ON, AMSC);
  delcr(line);
  strcpy(current_name->name, line);

  /* input NTYPE ID pair lines until no more, save patch id's that come in */
  for(i = iv = 0; i < KC-1; i++) {	/* loop on lines */
    for(j = 0; j < 5 && iv < IV/2; j++, iv++) { /* loop on items */
      fscanf(stream, "%d %d", &ntype, &id);
      if(ntype == 3) {		/* if its a patch, save ID */
	if(current_patch == NULL) { /* if 1st patch */
	  CALLOC(current_name->patch_list, 1, SM_PATCH, ON, AMSC);
	  current_patch = current_name->patch_list;
	}
	else {
	  CALLOC(current_patch->next, 1, SM_PATCH, ON, AMSC);
	  current_patch = current_patch->next;
	}
	current_patch->ID = id;
	patch_cnt++;
      }
    }
  }
  if(patch_cnt == 0) {
    fprintf(stderr,
	    "\nname_data: conductor '%s'\n  has no patch - redo naming so that one is included.\n",
	    current_name->name);
    exit(0);
  }
}

/* This function checks for coordinate-wise equivalent grid points.  Each
   grid structure has a list of equivalent grids.  If all three coordinates
   from two grid points are within SMALL_NUMBER, defined in patran.h, then
   they are equivalent.  */

grid_equiv_check()
{
  GRID *grid_ptr_1, *grid_ptr_2;
  int i;

  /* First, allocate spaces for equivalent grid arrays. */
  grid_ptr_1 = start_grid;
  while (grid_ptr_1) {
    CALLOC(grid_ptr_1->equiv_ID, number_grids, int, ON, AMSC);
    grid_ptr_1->number_equiv_grids = 0;
    grid_ptr_1 = grid_ptr_1->next;
  }

  /* Begin search.  Grid N is compared with grids from N+1 through the end
     of the list.  */
  grid_ptr_1 = start_grid;
  while (grid_ptr_1) {
    grid_ptr_2 = grid_ptr_1->next;
    while (grid_ptr_2) {
      if (if_same_coord(grid_ptr_1->coord,grid_ptr_2->coord)) {
	*(grid_ptr_1->equiv_ID + grid_ptr_1->number_equiv_grids)
	  = grid_ptr_2->ID;
	*(grid_ptr_2->equiv_ID + grid_ptr_2->number_equiv_grids)
	  = grid_ptr_1->ID;
	(grid_ptr_1->number_equiv_grids)++;
	(grid_ptr_2->number_equiv_grids)++;
      }
      grid_ptr_2 = grid_ptr_2->next;
    }
    grid_ptr_1 = grid_ptr_1->next;
  }

  /* Print the equivalent grid information.
  grid_ptr_1 = start_grid;
  while (grid_ptr_1) {
    printf("\nGrid %d : (%d)", grid_ptr_1->ID, grid_ptr_1->number_equiv_grids);
    for (i=0; i<grid_ptr_1->number_equiv_grids; i++)
      printf (" %d ", *(grid_ptr_1->equiv_ID+i));
    grid_ptr_1 = grid_ptr_1->next;
  } */
}


int if_same_coord(coord_1, coord_2)
  double coord_1[3], coord_2[3];
{
  int i;

  for (i=0; i<3; i++)
    if (fabs(coord_1[i] - coord_2[i]) > SMALL_NUMBER) return 0;
  return 1;
}

/*
  makes 1st \n in a string = \0 and then deletes all trail/leading wh space
*/
char *delcr(str)
char *str;
{
  int i, j, k;
  for(k = 0; str[k] != '\0'; k++) if(str[k] == '\n') { str[k] = '\0'; break; }
  for(i = 0; str[i] == ' ' || str[i] == '\t'; i++); /* count leading spaces */
  if(i > 0) {
    for(j = 0; str[j+i] != '\0'; j++) str[j] = str[j+i];
    str[j] = '\0';
  }
  for(k--; str[k] == ' ' || str[k] == '\t'; k--) str[k] = '\0';
  return(str);
}

/****************************************************************************

  This section of code is responsible for assigning conductor numbers to
  all the patches.

****************************************************************************/

/* This function fills the table that shows the connectivity of patches.
   If two patches share at least one common corner point, then they are
   connected.  It is done by going through all the grid points and finding
   patches that are connected by the grid point.  The end result table is
   symmetric.  */

fill_patch_patch_table(patch_patch_table)
  int *patch_patch_table;
{
  int patch_count, patch_count_save, *current_table_ptr, *corner, i;
  GRID *grid_ptr;
  PATCH *patch_ptr;

  grid_ptr = start_grid;
  while (grid_ptr) {

    /* Patch_count is generic counter of current position in the patch array,
       start_patch.  Patch_count_save is index of the last patch that had
       the current grid as its corner.  */
    patch_count = 0;
    patch_count_save = 0;
    current_table_ptr = 0;
    patch_ptr = start_patch;

    while (patch_ptr) {
      corner = patch_ptr->corner;
      for (i=0; i<4; i++)
	if (if_same_grid(*corner++,grid_ptr)) {
	  if (current_table_ptr) {  /* Have we already found another patch
				       with the same grid as its corner?  */
	    *(current_table_ptr+patch_count)=1;
	    *(patch_patch_table + (patch_count * number_patches)
	      + patch_count_save)=1;
	  }
	  current_table_ptr = patch_patch_table + patch_count*number_patches;
	  patch_count_save = patch_count;
	}
      patch_ptr = patch_ptr->next;
      patch_count++;
    }
    grid_ptr = grid_ptr->next;
  }
}


/* Return 1 if ID matches grid_ptr's ID or IDs of its equivalent grids,
   and 0 otherwise. */

int if_same_grid(ID,grid_ptr)
  int ID;
  GRID *grid_ptr;
{
  int *equiv_ID, i;

  if ((grid_ptr->ID)==ID) return 1;
  else {
    equiv_ID = grid_ptr->equiv_ID;
    for (i=0; i<grid_ptr->number_equiv_grids; i++)
      if (ID == equiv_ID[i]) return 1;
    return 0;
  }
}


/* This function searches through the patch_patch_table and finds groups of
   patches that are connected only among themselves. */

assign_conductor(patch_patch_table)
  int *patch_patch_table;
{
  PATCH *patch_ptr;
  int patch_count=0, *current_table_ptr;

  conductor_count=1;

  /* Sets all the patches to conductor 0, meaning that it is yet to be
     assigned a conductor_ID.  */
  patch_ptr = start_patch;
  while (patch_ptr) {
    patch_ptr->conductor_ID = 0;
    patch_ptr = patch_ptr->next;
  }

  /* Current_table_ptr points the row that needs to be searched through.
     That row is associated with the current patch in need of a conductor
     number.  */
  current_table_ptr = patch_patch_table;
  patch_ptr = start_patch;
  while (patch_ptr) {
    if ((patch_ptr->conductor_ID) == 0) {  /* If the patch is not assigned
					      a conductor number. */
      patch_ptr->conductor_ID = conductor_count;
      depth_search(patch_patch_table,current_table_ptr,conductor_count);
      conductor_count++;
    }
    patch_count++;
    current_table_ptr = patch_patch_table + patch_count*number_patches;
    patch_ptr = patch_ptr->next;
  }

  /* Prints the conductor information.
  patch_ptr = start_patch;
  while (patch_ptr) {
    printf("\nPatch %d   Conductor %d",
	   patch_ptr->ID, patch_ptr->conductor_ID);
    patch_ptr = patch_ptr->next;
  } */
}


/* This function searches through patch_patch_table recursively to
   find all patches that are somehow connected the current patch. */

depth_search(patch_patch_table,current_table_ptr,conductor_count)
  int *patch_patch_table, *current_table_ptr,conductor_count;
{
  PATCH *patch_ptr;
  int i, *new_table_ptr;

  patch_ptr=start_patch;
  new_table_ptr=patch_patch_table;
  for (i=0; i<number_patches; i++) {
    if ((*(current_table_ptr+i)) != 0) {  /* If the current patch is connected
					     to i'th patch. */
      if (patch_ptr->conductor_ID == 0) {  /* If the patch is yet to be
					      assigned a conductor number. */
	patch_ptr -> conductor_ID = conductor_count;
	new_table_ptr = patch_patch_table+i*number_patches;

	/* Call depth_search recursively to continue searching for
	   connected patches. */
	depth_search(patch_patch_table,new_table_ptr,conductor_count);
      }
    }
    patch_ptr=patch_ptr->next;
  }
}

/*
  used with new naming functions---finds the patran name in the patran list
  - this code used to be in mksCapDump()
*/
char *getPatranName(cond_num)
int cond_num;
{
  NAME *cname;

  cname = start_name_this_time;
  while(cname != NULL) {
    if((cname->patch_list)->conductor_ID == cond_num) return(cname->name);
    else cname = cname->next;
  }

  /*strcpy(cond_name, "??UNKNOWN??");*/
  fprintf(stdout, "getPatranName: conductor %d has no name\n", cond_num);
  return(NULL);

}

/****************************************************************************

  The following functions create the linked list of charges that can be
  used in Keith's FastCap program.

****************************************************************************/

charge *make_charges_all_patches(name_list, num_cond, surf_type, name_suffix)
Name **name_list;		/* master list of conductor names */
int *num_cond;			/* master conductor counter */
int surf_type;
char *name_suffix;
{
  CFEG *cfeg_ptr;
  int NELS, LPH_ID,conductor_ID,*element_list;
  char cond_name[BUFSIZ];
  PATCH *patch_ptr;
  charge *first_pq=0,*current_pq,*make_charges_patch();

  cfeg_ptr = start_cfeg;
  while (cfeg_ptr) {
    if (cfeg_ptr->LPH == 3) {
      NELS = cfeg_ptr->NELS;
      LPH_ID = cfeg_ptr->LPH_ID;

      /* Find the patch structure that is associated with the current cfeg
	 pointer in order to find the conductor number. */
      patch_ptr = start_patch;
      while (patch_ptr) {
	if (patch_ptr->ID == LPH_ID) {
	  if(surf_type == CONDTR || surf_type == BOTH) {
	    strcpy(cond_name, getPatranName(patch_ptr->conductor_ID));
	    strcat(cond_name, name_suffix);
	    conductor_ID = getConductorNum(cond_name, name_list, num_cond);
	  }
	  else conductor_ID = 0;
	  break;
	}
	patch_ptr = patch_ptr->next;
      }
/*      printf("\nCEFG %d  LPH %d LPH_ID %d Conductor_ID %d",
	     cfeg_ptr->ID,cfeg_ptr->LPH,cfeg_ptr->LPH_ID,conductor_ID); */

      /* For each patch, call the subroutine to handle the detail.
         Make sure all the lists of charges are linked. */
      element_list = cfeg_ptr->element_list;
      if (!first_pq) {
	first_pq = make_charges_patch(NELS,element_list,conductor_ID);
	current_pq = first_pq + NELS - 1;
      }
      else {
	current_pq->next =
	  make_charges_patch(NELS,element_list,conductor_ID);
	current_pq = (current_pq->next) + NELS - 1;
      }
    }
    cfeg_ptr=cfeg_ptr->next;
  }
  /* Put a nil pointer at the end.  */
  current_pq->next = 0;
  return first_pq;
}


/* This function creates the linked list of charges for a single patch. */

charge *make_charges_patch(NELS,element_list,conductor_ID)
  int NELS, *element_list, conductor_ID;
{
  charge *pq, *current_pq;
  int i,element_number,*element_corner_ptr;
  ELEMENT *element_ptr;
  NODE *node_ptr;

  CALLOC(pq,NELS,charge, ON, AMSC);

  /* Make sure that they are linked. */
  current_pq = pq;
  for (i=0; i<NELS-1; i++) {
    current_pq = pq + i;
    current_pq->next = current_pq+1;
  }

  /* NELS stands for number of elements. */
  for (i=0; i<NELS; i++) {
    (pq+i)->cond = conductor_ID;

    /* Original element number in Neutral file can be a negative number.  */
    if ((element_number= *(element_list+i))<0) element_number= -element_number;
    element_ptr = element_search_table[element_number];
    element_corner_ptr = element_ptr->corner;

    /* Pointers to the corner points' coordinates are set.  */

    if ((element_ptr->shape) == 4) {  /* Quadrilateral panels. */
      (pq+i)->shape = 4;
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[0], node_ptr->coord);
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[1], node_ptr->coord);
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[2], node_ptr->coord);
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[3], node_ptr->coord);
    }
    else {  /* Triangular panels. */
/*printf("\nTTT\n");*/
      (pq+i)->shape = 3;
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[0], node_ptr->coord);
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[1], node_ptr->coord);
      node_ptr = node_search_table[*(element_corner_ptr++)];
      VCOPY((pq+i)->corner[2], node_ptr->coord);
    }
  }
  return pq;
}


/*
  assigns correct conductor number to all patches in name structs
  - really already done implicitly in name_data by setting up sm_patch lists
  - conductor_ID of first patch is used as the number associated w/the name
  - checks for no names and names including several conductor's panels
  - checks one linked list against another, potentially n^2 => named
    regions should be kept small (as few patches as possible)
*/
assign_names()
{
  int quit, current_conductor, cnt = 0;
  PATCH *current_patch;
  SM_PATCH *current_name_patch;
  NAME *cur_name = start_name_this_time;

  if(start_name_this_time == NULL) {
    fprintf(stderr, "\nassign_names: no conductor names specified\n");
    exit(0);
  }

  /* for each name struct, find cond no of each patch (can be n^2 loop) */
  while(cur_name != NULL) {
    current_name_patch = cur_name->patch_list;
    current_conductor = 0;
    while(current_name_patch != NULL) {
      current_patch = start_patch;
      quit = 0;
      while(current_patch != NULL && quit == 0) {
	if(current_patch->ID == current_name_patch->ID) {
	  current_name_patch->conductor_ID = current_patch->conductor_ID;
	  if(current_conductor == 0) { /* if this is 1st name struct patch */
	    current_conductor = current_patch->conductor_ID;
	  }
	  else if(current_conductor != current_patch->conductor_ID) {
	    fprintf(stderr,
		    "\nassign_names: alleged conductor '%s'\n  has patches from more than one conductor - rename more carefully\n", cur_name->name);
	    exit(0);
	  }
	  quit = 1;
	}
	current_patch = current_patch->next;
      }
      if(quit == 0) {
	fprintf(stderr, "\nassign_names: in conductor '%s'\n  can't find named patch in master list\n", cur_name->name);
	exit(0);
      }
      current_name_patch = current_name_patch->next;
    }
    cur_name = cur_name->next;
    cnt++;
  }

  /* check to see if all conductors have a name and if too many names */
  if(cnt < conductor_count - 1) {
    fprintf(stderr, "\nassign_names: %d conductors have no names\n",
	    conductor_count - 1 - cnt);
    exit(0);
  }
  if(cnt > conductor_count - 1) {
    fprintf(stderr, "\nassign_names: %d names given for %d conductors\n",
	    cnt, conductor_count - 1);
    exit(0);
  }


}
