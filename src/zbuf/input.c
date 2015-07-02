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

#include "mulGlobal.h"
#include "zbufGlobal.h"
// Enrico
#include <string.h>

#ifdef SOLARIS
#include <sys/systeminfo.h>
#endif

/*
  reads an input file list file (a list of dielectric i/f and conductor
    surface files with permittivities)
  returns linked list of file pointers and permittivities in surface structs
  surface list file is specified on the command line with `-l<filename>'
  each line in the list file specifies a surface filename and its permittivites
  if a list file line has the filename string `stdin' then stdin will be
    read (note that more than one `stdin' is not allowed)

  list file line formats
  conductor surface:
  C <filename> <outer rel permittivity> <tx> <ty> <tz> [+]

  dielectric surface:
  D <file> <outer rel perm> <inner rel perm> <tx> <ty> <tz> <rx> <ry> <rz> [-]

  thin conductor on dielectric boundary (B => "both"):
  B <file> <outer rel perm> <inner rel perm> <tx> <ty> <tz> <rx> <ry> <rz> [+-]

  group name specification line:
  G <group name>

  comment line:
  * <comment line>

  the <tx> <ty> <tz> are 3 components of a translation vector that is applied
    to all the panels in the file
  the <rx> <ry> <rz> specify a reference point on the outside of the
    interface surface (all surface normals should point towards the point)
    the optional `-' indicates that the reference point is inside
    the surface (all surface normals should point away from the point)
    the reference point is used to figure which permittivity is on which side
    for type B surfaces, there must be no space between the `+' and `-' if
      both are used - note that D surfaces must never have a `+'
    since the reference point must be on one side, each file must contain
      a convex surface (ie any surface must be broken into such subsurfaces)
  the optional `+' indicates that the panels in the next conductor line file
    should be grouped together for the purposes of renumbering
    - if two files have two distinct conductors but both sets of panels
      are numbered `1' inside the two files then use something like
      C first_file 1.5 <tx> <ty> <tz>
      C second_file 1.5 <tx> <ty> <tz>
    - on the other hand, if parts of the same conductor are split into
      two files (say because the second part borders a different dielectric)
      then use something like
      C first_file 1.5 <tx> <ty> <tz> +
      C second_file 2.0 <tx> <ty> <tz>
      in this case it is up to the user to make sure first_file's panels
      and second_file's panels all have the same conductor number
    - to disable the renumbering entirely, use the `+' on all the
      conductor lines:
      C first_file 3.0 <tx> <ty> <tz> +
      C second_file 4.0 <tx> <ty> <tz> +
      C last_file 3.0 <tx> <ty> <tz> +
    - files grouped together with the + option have their conductor names
      appended with the string ` (GROUP<number>)'
      - for example, the conductor name `BIT_LINE' shows up as
        `BIT_LINE (GROUP3)' if it is in the third group
      - a string other than `GROUP<number>' may be specified for the
        group name using G line `G <group name>' just before the group to
	be renamed; this is helpful when idenifying conductors to omit
	from capacitance calculations using the -k option
*/
void read_list_file(surf_list, num_surf, list_file, read_from_stdin)
int *num_surf, read_from_stdin;
char *list_file;
surface **surf_list;
{
  int linecnt, end_of_chain, ref_pnt_is_inside, group_cnt;
  FILE *fp, *fopen();
  char tline[BUFSIZ], file_name[BUFSIZ], plus[BUFSIZ], group_name[BUFSIZ];
  double outer_perm, inner_perm, tx, ty, tz, rx, ry, rz;
  surface *cur_surf;

  /* find the end of the current surface list */
  if(*surf_list != NULL) {
    for(cur_surf = *surf_list; cur_surf->next != NULL;
	cur_surf = cur_surf->next);
  }

  /* attempt to open file list file */
  if((fp = fopen(list_file, "r")) == NULL) {
    fprintf(stderr, "read_list_file: can't open list file\n  `%s'\nto read\n",
	    list_file);
    exit(0);
  }

  /* read file names and permittivities, build linked list */
  linecnt = 0;
  group_cnt = read_from_stdin + 1;
  sprintf(group_name, "GROUP%d", group_cnt);
  while(fgets(tline, sizeof(tline), fp) != NULL) {
    linecnt++;
    if(tline[0] == 'C' || tline[0] == 'c') {
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf",
		file_name, &outer_perm, &tx, &ty, &tz) != 5) {
	fprintf(stderr,
	       "read_list_file: bad conductor surface format, tline %d:\n%s\n",
		linecnt, tline);
	exit(0);
      }

      /* check if end of chain of surfaces with same conductor numbers */
      end_of_chain = TRUE;
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf %s",
		file_name, &outer_perm, &tx, &ty, &tz, plus) == 6) {
	if(!strcmp(plus, "+")) end_of_chain = FALSE;
      }

      /* allocate and load surface struct */
      if(*surf_list == NULL) {
	CALLOC(*surf_list, 1, surface, ON, AMSC);
	cur_surf = *surf_list;
      }
      else {
	CALLOC(cur_surf->next, 1, surface, ON, AMSC);
	cur_surf->next->prev = cur_surf;
	cur_surf = cur_surf->next;
      }

      cur_surf->type = CONDTR;
      cur_surf->trans[0] = tx;
      cur_surf->trans[1] = ty;
      cur_surf->trans[2] = tz;
      cur_surf->end_of_chain = end_of_chain;
      CALLOC(cur_surf->name, strlen(file_name)+1, char, ON, AMSC);
      strcpy(cur_surf->name, file_name);
      cur_surf->outer_perm = outer_perm;

      /* set up group name */
      CALLOC(cur_surf->group_name, strlen(group_name)+1, char, ON, AMSC);
      strcpy(cur_surf->group_name, group_name);

      /* update group name if end of chain */
      if(end_of_chain) {
	sprintf(group_name, "GROUP%d", ++group_cnt);
      }

      (*num_surf)++;
    }
    else if(tline[0] == 'B' || tline[0] == 'b') {
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf %lf %lf %lf %lf",
		file_name, &outer_perm, &inner_perm, &tx, &ty, &tz,
		&rx, &ry, &rz) != 9) {
	fprintf(stderr,
		"read_list_file: bad thin conductor on dielectric interface surface format, line %d:\n%s\n",
		linecnt, tline);
	exit(0);
      }

      /* check if end of chain of surfaces with same conductor numbers */
      end_of_chain = TRUE;
      ref_pnt_is_inside = FALSE;
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf %lf %lf %lf %lf %s",
		file_name, &outer_perm, &inner_perm, &tx, &ty, &tz,
		&rx, &ry, &rz, plus)
	 == 10) {
	if(!strcmp(plus, "+")) end_of_chain = FALSE;
	if(!strcmp(plus, "+-") || !strcmp(plus, "-+")) {
	  end_of_chain = FALSE;
	  ref_pnt_is_inside = TRUE;
	}
	if(!strcmp(plus, "-")) ref_pnt_is_inside = TRUE;
      }

      /* allocate and load surface struct */
      if(*surf_list == NULL) {
	CALLOC(*surf_list, 1, surface, ON, AMSC);
	cur_surf = *surf_list;
      }
      else {
	CALLOC(cur_surf->next, 1, surface, ON, AMSC);
	cur_surf->next->prev = cur_surf;
	cur_surf = cur_surf->next;
      }

      cur_surf->type = BOTH;
      cur_surf->trans[0] = tx;
      cur_surf->trans[1] = ty;
      cur_surf->trans[2] = tz;
      cur_surf->ref[0] = rx;
      cur_surf->ref[1] = ry;
      cur_surf->ref[2] = rz;
      cur_surf->ref_inside = ref_pnt_is_inside;
      cur_surf->end_of_chain = end_of_chain;
      CALLOC(cur_surf->name, strlen(file_name)+1, char, ON, AMSC);
      strcpy(cur_surf->name, file_name);
      cur_surf->outer_perm = outer_perm;
      cur_surf->inner_perm = inner_perm;

      /* set up group name */
      CALLOC(cur_surf->group_name, strlen(group_name)+1, char, ON, AMSC);
      strcpy(cur_surf->group_name, group_name);

      /* update group name if end of chain */
      if(end_of_chain) {
	sprintf(group_name, "GROUP%d", ++group_cnt);
      }

      (*num_surf)++;
    }
    else if(tline[0] == 'D' || tline[0] == 'd') {
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf %lf %lf %lf %lf",
		file_name, &outer_perm, &inner_perm, &tx, &ty, &tz,
		&rx, &ry, &rz) != 9) {
	fprintf(stderr,
		"read_list_file: bad dielectric interface surface format, line %d:\n%s\n",
		linecnt, tline);
	exit(0);
      }

      /* check to see if reference point is negative side of surface */
      ref_pnt_is_inside = FALSE;
      if(sscanf(&(tline[1]), "%s %lf %lf %lf %lf %lf %lf %lf %lf %s",
		file_name, &outer_perm, &inner_perm, &tx, &ty, &tz,
		&rx, &ry, &rz, plus)
	 == 10) {
	if(!strcmp(plus, "-")) ref_pnt_is_inside = TRUE;
      }

      /* allocate and load surface struct */
      if(*surf_list == NULL) {
	CALLOC(*surf_list, 1, surface, ON, AMSC);
	cur_surf = *surf_list;
      }
      else {
	CALLOC(cur_surf->next, 1, surface, ON, AMSC);
	cur_surf->next->prev = cur_surf;
	cur_surf = cur_surf->next;
      }

      cur_surf->type = DIELEC;
      cur_surf->trans[0] = tx;
      cur_surf->trans[1] = ty;
      cur_surf->trans[2] = tz;
      cur_surf->ref[0] = rx;
      cur_surf->ref[1] = ry;
      cur_surf->ref[2] = rz;
      cur_surf->ref_inside = ref_pnt_is_inside;
      cur_surf->end_of_chain = TRUE;
      CALLOC(cur_surf->name, strlen(file_name)+1, char, ON, AMSC);
      strcpy(cur_surf->name, file_name);
      cur_surf->outer_perm = outer_perm;
      cur_surf->inner_perm = inner_perm;

      /* set up group name */
      CALLOC(cur_surf->group_name, strlen(group_name)+1, char, ON, AMSC);
      strcpy(cur_surf->group_name, group_name);

      /* update group name (DIELEC surface is always end of chain) */
      sprintf(group_name, "GROUP%d", ++group_cnt);

      (*num_surf)++;
    }
    else if(tline[0] == 'G' || tline[0] == 'g') {
      if(sscanf(&(tline[1]), "%s", group_name) != 1) {
	fprintf(stderr,"read_list_file: bad group name format, line %d:\n%s\n",
		linecnt, tline);
	exit(0);
      }
    }
    else if(tline[0] == '%' || tline[0] == '*' ||
	    tline[0] == '#'); /* ignore comments */
    else {
      fprintf(stderr, "read_list_file: bad line format, line %d:\n%s\n",
		linecnt, tline);
      exit(0);
    }
  }
  fclose(fp);

}

#if 1 == 0			/* now done one panel at a time in initcalcp */

/*
  TEMPORARY: works best for closed surfaces in which case ref should
    be a point inside with ref_inside = TRUE
  can also be used with open surfaces: ref = point on - side, ref_inside = TRUE
  will not work if the surface folds back on itself (and other ways too)
*/
void align_normals(panel_list, surf)
surface *surf;
charge *panel_list;
{
  int i, flip_normal;
  char *hack_path();
  charge *nc;
  double ctr_minus_n[3], ctr_plus_n[3], norm_minus, norm_plus, norm, norm_sq;
  double x, y, z, *normal, *direction, *tempd;
  int ref_inside = surf->ref_inside;
  double *ref = surf->ref, temp;
  char *surf_name = surf->name;

  for(nc = panel_list; nc != NULL; nc = nc->next) {

    /* get panel position (relative to reference point) and normal */
    x = nc->x - ref[0]; y = nc->y - ref[1]; z = nc->z - ref[2];
    norm_sq = x*x + y*y + z*z;
    norm = sqrt(norm_sq);
    normal = nc->Z;

    /* add the (scaled) normal and negative normal to the panel center */
    /* negative normal result should be closer to ref point if(ref_inside) */
    ctr_minus_n[0] = x - 0.1*norm*normal[0];
    ctr_minus_n[1] = y - 0.1*norm*normal[1];
    ctr_minus_n[2] = z - 0.1*norm*normal[2];
    ctr_plus_n[0] = x + 0.1*norm*normal[0];
    ctr_plus_n[1] = y + 0.1*norm*normal[1];
    ctr_plus_n[2] = z + 0.1*norm*normal[2];

    /* get norms of test points, one inside (minus) other out (plus) */
    norm_minus = ctr_minus_n[0]*ctr_minus_n[0];
    norm_plus = ctr_plus_n[0]*ctr_plus_n[0];
    for(i = 1; i < 3; i++) {
      norm_minus += ctr_minus_n[i]*ctr_minus_n[i];
      norm_plus += ctr_plus_n[i]*ctr_plus_n[i];
    }

    flip_normal = FALSE;
    if(norm_minus > norm_sq) {
      if(norm_plus > norm_sq) {
	fprintf(stderr,
		"align_normals: both test points on non-reference side\n");
	fprintf(stderr, "  Surface: %s\n", hack_path(surf_name));
	fprintf(stderr, "  Translation: (%g %g %g)\n", surf->trans[0],
		surf->trans[1], surf->trans[2]);
	fprintf(stderr, "  Reference point: (%g %g %g)\n",
		ref[0], ref[1], ref[2]);
	fprintf(stderr, "  Panel cntr: (%g %g %g)\n",
		nc->x, nc->y, nc->z);
	fprintf(stderr, "  Normal: (%g %g %g)\n",
		normal[0], normal[1], normal[2]);
	exit(0);
      }
      if(ref_inside) flip_normal = TRUE;
    }
    else if(norm_plus < norm_sq) {
      if(norm_minus < norm_sq) {
	fprintf(stderr,
		"align_normals: both test points on reference point side\n");
	fprintf(stderr, "  Surface: %s\n", hack_path(surf_name));
	fprintf(stderr, "  Translation: (%g %g %g)\n", surf->trans[0],
		surf->trans[1], surf->trans[2]);
	fprintf(stderr, "  Reference point: (%g %g %g)\n",
		ref[0], ref[1], ref[2]);
	fprintf(stderr, "  Panel cntr: (%g %g %g)\n",
		nc->x, nc->y, nc->z);
	fprintf(stderr, "  Normal: (%g %g %g)\n",
		normal[0], normal[1], normal[2]);
	exit(0);
      }
      if(!ref_inside) flip_normal = TRUE;
    }

    if(flip_normal) {
      for(i = 0; i < 3; i++) {
	normal[i] = -normal[i];	/* flip the normal */
	nc->X[i] = -(nc->X[i]);	/* flip the x direction */
	/* interchange points 0 and 2 so that corner order will be
	   consistent with X flip (note that this is OK for quads and tris) */
	temp = nc->corner[0][i];
	nc->corner[0][i] = nc->corner[2][i];
	nc->corner[2][i] = temp;
      }
    }
  }
}

#endif

/*
  add dummy panel structs to the panel list for electric field evaluation
  - assumes its handed a list of DIELEC or BOTH type panels
*/
void add_dummy_panels(panel_list)
charge *panel_list;
{
  double h;
  charge *dummy_list = NULL;
  charge *cur_panel, *cur_dummy;

  for(cur_panel = panel_list; cur_panel != NULL; cur_panel = cur_panel->next) {
    cur_panel->dummy = FALSE;

    /* make 2 dummy panels for evaluation points needed to do div difference */
    /* make the first */
    if(dummy_list == NULL) {
      CALLOC(dummy_list, 1, charge, ON, AMSC);
      cur_dummy = dummy_list;
    }
    else {
      CALLOC(cur_dummy->next, 1, charge, ON, AMSC);
      cur_dummy = cur_dummy->next;
    }

    cur_dummy->dummy = TRUE;
    h = HPOS;
    cur_dummy->x = cur_panel->x + cur_panel->Z[0]*h;
    cur_dummy->y = cur_panel->y + cur_panel->Z[1]*h;
    cur_dummy->z = cur_panel->z + cur_panel->Z[2]*h;
    /* note ABUSE OF area field - used to store div dif distance */
    cur_dummy->area = h;

    cur_panel->pos_dummy = cur_dummy; /* link dummy to its real panel */

    /* make the second dummy struct */
    CALLOC(cur_dummy->next, 1, charge, ON, AMSC);
    cur_dummy = cur_dummy->next;

    cur_dummy->dummy = TRUE;
    h = HNEG;
    cur_dummy->x = cur_panel->x - cur_panel->Z[0]*h;
    cur_dummy->y = cur_panel->y - cur_panel->Z[1]*h;
    cur_dummy->z = cur_panel->z - cur_panel->Z[2]*h;
    /* note ABUSE OF area field - used to store div dif distance */
    cur_dummy->area = h;

    cur_panel->neg_dummy = cur_dummy; /* link dummy to its real panel */
  }

  /* put the dummies in the list */
  for(cur_panel = panel_list; cur_panel->next != NULL;
      cur_panel = cur_panel->next);
  cur_panel->next = dummy_list;

}

/* returns a pointer to a file name w/o the path (if present) */
char *hack_path(str)
char *str;
{
  int i;
  int last_slash;

  for(i = last_slash = 0; str[i] != '\0'; i++) {
    if(str[i] == '/') last_slash = i;
  }

  if(str[last_slash] == '/') return(&(str[last_slash+1]));
  else return(str);
}

/*
  reassigns conductor numbers to a list of panels so that they'll
    be numbered contiguously from 1
  - also changes conductor numbers associated with conductor name structs
  - dummy panels are skipped
  - dielectric panels, with conductor number 0, are also skipped
*/
void reassign_cond_numbers(panel_list, name_list, surf_name)
char *surf_name;
NAME *name_list;
charge *panel_list;
{
  int i, j, cond_nums[MAXCON], num_cond, cond_num_found, temp;
  char str[BUFSIZ], *hack_path();
  charge *cur_panel;
  NAME *cur_name;

  /* get the conductor numbers currently being used */
  num_cond = 0;
  for(cur_panel = panel_list; cur_panel != NULL; cur_panel = cur_panel->next) {
    if(cur_panel->dummy || cur_panel->cond == 0) continue;

    cond_num_found = FALSE;
    for(i = 0; i < num_cond; i++) {
      if(cur_panel->cond == cond_nums[i]) {
	cond_num_found = TRUE;
	break;
      }
    }
    if(!cond_num_found) cond_nums[num_cond++] = cur_panel->cond;
  }

  /* rewrite all the conductor numbers to be their position in the array */
  for(cur_panel = panel_list; cur_panel != NULL; cur_panel = cur_panel->next) {
    if(cur_panel->dummy || cur_panel->cond == 0) continue;

    for(i = 0; i < num_cond && cur_panel->cond != cond_nums[i]; i++);
    if(i == num_cond) {
      fprintf(stderr,
       "reassign_cond_numbers: cant find conductor number that must exist\n");
      exit(0);
    }
    cur_panel->cond = i+1;
  }

  /* do the same for the name structs */
  for(cur_name = name_list; cur_name != NULL; cur_name = cur_name->next) {
    for(i = 0;
	i < num_cond && cur_name->patch_list->conductor_ID != cond_nums[i];
	i++);
    if(i == num_cond) {
      fprintf(stderr,
	   "reassign_cond_numbers: cant find conductor number in name list\n");
      exit(0);
    }
#if 1 == 0
    /* change the name given to this conductor if it was derived from the ID */
    /* check the number */
    if(sscanf(&(cur_name->name[9]), "%d", &temp) == 1) {
      if(temp == cur_name->patch_list->conductor_ID) {
	strcpy(str, cur_name->name);
	str[9] = '\0';
	/* check the rest of the string, replace if necessary */
	if(!strcmp(str, "CONDUCTOR")) {
	  sprintf(str, "COND%d (%s)", i+1, hack_path(surf_name));
	  if(strlen(str) > strlen(cur_name->name))
	      CALLOC(cur_name->name, strlen(str)+1, char, ON, AMSC);
	  strcpy(cur_name->name, str);
	}
      }
    }
    cur_name->patch_list->conductor_ID = i+1;
#endif
  }


}

/*
  negates all the conductor numbers - used to make a panel list's conds unique
    just before renumbering
*/
void negate_cond_numbers(panel_list, name_list)
NAME *name_list;
charge *panel_list;
{
  charge *cur_panel;
  NAME *cur_name;

  for(cur_panel = panel_list; cur_panel != NULL; cur_panel = cur_panel->next) {
    if(cur_panel->dummy) continue;

    cur_panel->cond = -cur_panel->cond;
  }

  for(cur_name = name_list; cur_name != NULL; cur_name = cur_name->next) {
    cur_name->patch_list->conductor_ID = -cur_name->patch_list->conductor_ID;
  }
}

/*
  for debug - dumps the iter list
*/
int dump_ilist()
{
  ITER *cur_iter;
  extern ITER *qpic_num_list;

  /* check the list for the iter number passed in */
  fprintf(stdout, "Iter list:");
  for(cur_iter = qpic_num_list; cur_iter != NULL; cur_iter = cur_iter->next) {
    fprintf(stdout, "%d ", cur_iter->iter);
  }
  fprintf(stdout, "\n");
  return(TRUE);
}

#if 1 == 0
/*
  adds an iteration number to the list that get shaded .ps file dumps
  - list is built on global variable q_iter
*/
void add_iter(iter_num)
int iter_num;
{
  ITER *cur_iter, *tail_iter;
  extern ITER *q_iter;

  /* check the list for the iter number passed in */
  for(cur_iter = q_iter; cur_iter != NULL;
      tail_iter = cur_iter, cur_iter = cur_iter->next) {
    if(cur_iter->iter == iter_num) {
      return;
    }
  }

  /* not in list; create a new iter struct to store the new iter number */
  if(q_iter == NULL) {
    CALLOC(q_iter, 1, ITER, ON, AMSC);
    tail_iter = q_iter;
  }
  else {
    CALLOC(tail_iter->next, 1, ITER, ON, AMSC);
    tail_iter = tail_iter->next;
  }
  tail_iter->iter = iter_num;
  tail_iter->next = NULL;
}
#endif

/*
  checks if a particular iter is in the list; returns TRUE if it is
*/
int want_this_iter(iter_list, iter_num)
ITER *iter_list;
int iter_num;
{
  ITER *cur_iter;

  for(cur_iter = iter_list; cur_iter != NULL; cur_iter = cur_iter->next) {
    if(cur_iter->iter == iter_num) {
      return(TRUE);
    }
  }

  return(FALSE);
}
/*
  sets up the ps file base string
*/
void get_ps_file_base(argv, argc)
char *argv[];
int argc;
{
  int i, j;
  char temp[BUFSIZ], *hack_path();
  extern char *ps_file_base, *in_file_name;

  /* - if no list file, use input file; otherwise use list file */
  /* - if neither present, use "stdin" */
  /*   check for list file */
  for(i = 1; i < argc; i++) {
    if(argv[i][0] == '-' && argv[i][1] == 'l') {
      strcpy(temp, &(argv[i][2]));
      /* go to end of string, walk back to first period */
      for(j = 0; temp[j] != '\0'; j++);
      for(; temp[j] != '.' && j >= 0; j--);
      if(temp[j] == '.') temp[j] = '\0';
      /* save list file base */
      CALLOC(ps_file_base, strlen(temp)+1, char, ON, AMSC);
      strcpy(ps_file_base, hack_path(temp));
      break;
    }
    else if(argv[i][0] != '-') { /* not an option, must be input file */
      /* modified for fasthenry to allow ../../fname.inp */
      strcpy(temp, argv[i]);
      for(j = strlen(temp); j != -1 && temp[j] != '.' && temp[j] != '/'; j--);
      if (j == -1 || temp[j] == '/')
        j = strlen(temp) + 1;
      temp[j] = '\0';
      /* save list file base */
      CALLOC(in_file_name, strlen(temp)+1, char, ON, AMSC);
      strcpy(in_file_name, temp);
      CALLOC(ps_file_base, strlen(temp)+1, char, ON, AMSC);
      strcpy(ps_file_base, hack_path(temp));
      break;
    }
  }

  if(ps_file_base == NULL) {	/* input must be stdin */
    CALLOC(ps_file_base, strlen("stdin")+1, char, ON, AMSC);
    strcpy(ps_file_base, "stdin");
  }
}

/*
  open all the surface files and return a charge (panel) struct list
  set up pointers from each panel to its corresponding surface struct
  align the normals of all the panels in each surface so they point
    towards the same side as where the ref point is (dielectric files only)
*/
charge *read_panels(surf_list, name_list, num_cond)
Name **name_list;
int *num_cond;
surface *surf_list;
{
  int patran_file, num_panels, stdin_read, num_dummies, num_quads, num_tris;
  charge *panel_list = NULL, *cur_panel, *patfront(), *panel_group, *c_panel;
  surface *cur_surf;
  extern NAME *start_name, *start_name_this_time;
  extern char *title;
  NAME *name_group;
  FILE *fp, *fopen();
  char surf_name[BUFSIZ], *hack_path();
  int patran_file_read;

  /*title[0] = '\0';*/
  stdin_read = FALSE;
  for(cur_surf = surf_list; cur_surf != NULL; cur_surf = cur_surf->next) {
    if(!strcmp(cur_surf->name, "stdin")) {
      if(stdin_read) {
	fprintf(stderr, "read_panels: attempt to read stdin twice\n");
	exit(0);
      }
      else {
	stdin_read = TRUE;
	fp = stdin;
      }
    }
    else if((fp = fopen(cur_surf->name, "r")) == NULL) {
      fprintf(stderr, "read_panels: can't open\n  `%s'\nto read\n",
	      cur_surf->name);
      exit(0);
    }

    /* input the panel list */
    if(panel_list == NULL) {
      /* group names are set up in read_list_file() */
      sprintf(surf_name, "%%%s", cur_surf->group_name);
      panel_list = cur_panel
	  = patfront(fp, &patran_file, cur_surf->type, cur_surf->trans,
		     name_list, num_cond, surf_name);
      patran_file_read = patran_file;
      panel_group = cur_panel;
      name_group = start_name;
    }
    else {
      if(cur_surf->prev->end_of_chain) {
	sprintf(surf_name, "%%%s", cur_surf->group_name);
	patran_file_read = FALSE;
      }
      cur_panel->next
	  = patfront(fp, &patran_file, cur_surf->type, cur_surf->trans,
		     name_list, num_cond, surf_name);
      if(!patran_file && patran_file_read) {
	fprintf(stderr, "read_panels: generic format file\n  `%s'\nread after neutral file(s) in same group---reorder list file entries\n", cur_surf->name);
	exit(0);
      }
      patran_file_read = patran_file;
      cur_panel = cur_panel->next;
      if(cur_surf->prev->end_of_chain) {
	/* if previous surface was the end of a chain, set up new group */
	panel_group = cur_panel;
	name_group = start_name_this_time; /* not really used anymore */
      }
    }
    if(strcmp(cur_surf->name, "stdin") != 0) fclose(fp);

    cur_surf->panels = cur_panel;

    /* save the surface file's title */
    CALLOC(cur_surf->title, strlen(title)+1, char, ON, AMSC);
    strcpy(cur_surf->title, title);
    title[0] = '\0';		/* not sure if needed */

    /* if the surface is a DIELEC, make sure all conductor numbers are zero */
    /* - also link each panel to its surface */
    for(c_panel = cur_panel; c_panel != NULL; c_panel = c_panel->next) {
      if(cur_surf->type == DIELEC) c_panel->cond = 0;
      c_panel->surf = cur_surf;
    }

    /* align the normals and add dummy structs if dielec i/f */
    initcalcp(cur_surf->panels);/* get normals, edges, perpendiculars */
    if(cur_surf->type == DIELEC || cur_surf->type == BOTH) {
      /* if(patran_file) align_normals(cur_surf->panels);
      align_normals(cur_surf->panels, cur_surf); /* now done in calcp */
      add_dummy_panels(cur_surf->panels); /* add dummy panels for field calc */
    }

    /* make cur_panel = last panel in list, count panels */
    num_panels = num_dummies = num_tris = num_quads = 0;
    for(cur_panel = cur_surf->panels ; ; cur_panel = cur_panel->next) {
      num_panels++;
      if(cur_panel->dummy) num_dummies++;
      else if(cur_panel->shape == 3) num_tris++;
      else if(cur_panel->shape == 4) num_quads++;
      else {
	fprintf(stderr, "read_panels: bad panel shape, %d\n",
		cur_panel->shape);
	exit(0);
      }
      if(cur_panel->next == NULL) break;
    }

    /*fprintf(stdout, "Surface %s has %d quads and %d tris\n",
	    cur_surf->name, num_quads, num_tris);*/

    cur_surf->num_panels = num_panels;
    cur_surf->num_dummies = num_dummies;

#if 1 == 0			/* now done implicitly with suffix names */
    if(cur_surf->type == CONDTR || cur_surf->type == BOTH) {
      if(cur_surf->end_of_chain) {
	/* if the current surface is the end of a group to be #ed together, */
	/*   renumber the conductors in the newly input surface so numbering */
	/*   is unique and contiguous over the whole list */
	if(cur_surf != surf_list) { /* if this is not the first surface */
	  negate_cond_numbers(panel_group, name_group);
	  reassign_cond_numbers(panel_list, start_name, cur_surf->name);
	}
	else reassign_cond_numbers(panel_list, start_name, cur_surf->name);
      }
    }
#endif

  }
  return(panel_list);
}

/*
  returns either conductor number or one of two error codes
  NOTUNI => no group name given and name by itself is not unique
  NOTFND => neither name by itself nor with group name is not in list
  - any unique leading part of the name%group_name string may be specified
*/
int getUniqueCondNum(name, name_list)
char *name;
Name *name_list;
{
  int nlen, cond;
  char name_frag[BUFSIZ], *last_alias(), *cur_alias;
  Name *cur_name, *prev_name;
  int i, j, alias_match_name(), times_in_list;

  nlen = strlen(name);
  times_in_list = 0;

  /* fish through name list for name---check first nlen chars for match */
  for(cur_name = name_list, i = 1; cur_name != NULL && times_in_list < 2;
      cur_name = cur_name->next, i++) {
    cur_alias = last_alias(cur_name);
    for(j = 0; j < nlen; j++) name_frag[j] = cur_alias[j];
    name_frag[j] = '\0';
    if(!strcmp(name_frag, name)) {
      times_in_list++;	/* increment times name in list count */
      cond = i;
    }
    prev_name = cur_name;
  }

  /* name can't be dealt with; return appropriate error code */
  if(times_in_list > 2) return(NOTUNI);
  else if(times_in_list == 1) return(cond);
  else return(NOTFND);
}


/*
  called after all conductor names have been resolved to get list
    of conductor numbers that whose columns will not be calculated
  parses the conductor kill list spec from command line arg saved before
  (conds that get no solve); puts result in kill_num_list
  list string format:
  [<cond name>[%<group name>]],[<cond name>[%<group name>]]...
  - no spaces; group name may be omitted if conductor name is unique
    across all groups
  - conductor names can't have any %'s
  - redundant names are detected as errors
*/
ITER *get_kill_num_list(name_list, kill_name_list)
char *kill_name_list;
Name *name_list;
{
  int i, j, start_token, end_token, end_name, cond;
  char name[BUFSIZ], group_name[BUFSIZ];
  ITER *kill_num_list;
  ITER *cur_cond;

  /* check for no name list given */
  if(kill_name_list == NULL) return(NULL);

  start_token = 0;
  kill_num_list = NULL;
  while(kill_name_list[start_token] != '\0') {

    /* loop until next comma or end of list */
    for(i = start_token; kill_name_list[i] != '\0' && kill_name_list[i] != ',';
	i++);
    end_token = i;

    /* extract the name%group_name string */
    /*   copy the name */
    for(i = start_token, j = 0; i < end_token; i++, j++)
	name[j] = kill_name_list[i];
    name[j] = '\0';

    /* attempt to get conductor number from name and group_name */
    cond = getUniqueCondNum(name, name_list);
    if(cond == NOTUNI) {
      fprintf(stderr,
   "get_kill_num_list: cannot find unique conductor name starting `%s'\n",
	      name);
      exit(0);
    }
    else if(cond == NOTFND) {
      fprintf(stderr,
	      "get_kill_num_list: cannot find conductor name starting `%s'\n",
	      name);
      exit(0);
    }

    /* add conductor name to list of conductors to omit */
    if(kill_num_list == NULL) {
      CALLOC(kill_num_list, 1, ITER, ON, AMSC);
      cur_cond = kill_num_list;
    }
    else {
      CALLOC(cur_cond->next, 1, ITER, ON, AMSC);
      cur_cond = cur_cond->next;
    }
    cur_cond->iter = cond;

    if(kill_name_list[end_token] == ',') start_token = end_token+1;
    else start_token = end_token;
  }
  return(kill_num_list);
}

/*
  command line parsing routine
*/
void parse_command_line(argv, argc, autmom, autlev, relperm, numMom, numLev,
			input_file, surf_list_file, read_from_stdin)
int argc, *autmom, *autlev, *numMom, *numLev, *read_from_stdin;
double *relperm;
char *argv[], **input_file, **surf_list_file;
{
  int cmderr, i;
  char **chkp, *chk;
  long strtol();
  extern char *kill_name_list, *kinp_name_list;
  extern ITER *kill_num_list, *kinp_num_list;
  extern double iter_tol;

  extern int s_, n_, g_, c_, x_, k_, rc_, rd_, rb_, q_, rk_, m_, f_, dd_;
  extern double view[], moffset[], rotation, distance, linewd, scale, axeslen;
  extern double elevation, azimuth;
  extern int up_axis;
  extern char *line_file;
  extern char *ps_file_base;
  extern ITER *qpic_num_list;
  extern char *qpic_name_list;
  extern ITER *kq_num_list;
  extern char *kq_name_list;
  /* load default parameters */
  azimuth = DEFAZM;             /* azimuth */
  elevation = DEFELE;           /* elevation */
  rotation = DEFROT;            /* rotation relative to image of z axis */
  distance = DEFDST;            /* distance to view pnt = (1+distance)radius */
  moffset[0] = OFFSETX;         /* puts the origin this dist from lower left */
  moffset[1] = OFFSETY;
  scale = DEFSCL;               /* master scaling - applied to 2d image */
  linewd = DEFWID;              /* line width used in ps file */
  axeslen = DEFAXE;             /* length of axes lines in 3d */
  up_axis = DEFUAX;             /* upward-pointing axis in 2d image */
  line_file = NULL;             /* file of lines/arrows in .fig format */
  ps_file_base = NULL;		/* base used to form .ps file names */
  qpic_num_list = NULL;		/* list of cond nums to get shaded plots for */
  qpic_name_list = NULL;	/* list of cond names to get shaded plots */
  kq_num_list = NULL;		/* list of cond nums in shaded plots */
  kq_name_list = NULL;		/* list of cond names in shaded plots */
  s_ = n_ = g_ = c_ = x_ = k_ = rc_ = rd_ = rb_ = q_ = rk_ = m_ = f_ = FALSE;
  dd_ = FALSE;

  /* for fasthenry, all we ever do is visualization */
  /* m_ = TRUE; */

  iter_tol = ABSTOL;
  kill_num_list = kinp_num_list = NULL;
  kill_name_list = kinp_name_list = NULL;
  cmderr = FALSE;
  chkp = &chk;			/* pointers for error checking */

  for(i = 1; i < argc && cmderr == FALSE; i++) {
    if(argv[i][0] == '-') {
      if(argv[i][1] == 'd' && argv[i][2] == 'c') {
	dd_ = TRUE;
      }
      else if(argv[i][1] == 'l') {
	*surf_list_file = &(argv[i][2]);
      }
      else if(argv[i][1] == 'r' && argv[i][2] == 's') {
	kill_name_list = &(argv[i][3]);
      }
      else if(argv[i][1] == 'r' && argv[i][2] == 'i') {
	kinp_name_list = &(argv[i][3]);
      }
      else if(argv[i][1] == '\0') {
	*read_from_stdin = TRUE;
      }
/*#if CAPVEW == ON*/
      else if(argv[i][1] == 'f') {
	f_ = TRUE;
      }
      else if(argv[i][1] == 'b') {
	line_file = &(argv[i][2]);
      }
      else if(argv[i][1] == 'a') {
        if(sscanf(&(argv[i][2]), "%lf", &azimuth) != 1) {
	  fprintf(stderr, "%s: bad view point azimuth angle '%s'\n",
		  argv[0], &argv[i][2]);
	  cmderr = TRUE;
	  break;
	}
      }
      else if(argv[i][1] == 'e') {
        if(sscanf(&(argv[i][2]), "%lf", &elevation) != 1) {
	  fprintf(stderr, "%s: bad view point elevation angle '%s'\n",
		  argv[0], &argv[i][2]);
	  cmderr = TRUE;
	  break;
	}
      }
      else if(argv[i][1] == 'r' && argv[i][2] == 'c') {
	kq_name_list = &(argv[i][3]);
	rc_ = TRUE;
      }
      else if(!strcmp(&(argv[i][1]), "rd")) rd_ = TRUE;
      /*else if(!strcmp(&(argv[i][1]), "rb")) rb_ = TRUE;*/
      else if(!strcmp(&(argv[i][1]), "rk")) rk_ = TRUE;
      else if(argv[i][1] == 'r') {
        if(sscanf(&(argv[i][2]), "%lf", &rotation) != 1) {
	  fprintf(stderr, "%s: bad image rotation angle '%s'\n",
		  argv[0], &argv[i][2]);
	  cmderr = TRUE;
	  break;
	}
      }
      else if(argv[i][1] == 'h') {
        if(sscanf(&(argv[i][2]), "%lf", &distance) != 1) cmderr = TRUE;
	else if(distance <= 0.0) cmderr = TRUE;
	if(cmderr) {
	  fprintf(stderr, "%s: bad view point distance '%s'\n",
		  argv[0], &argv[i][2]);
	  break;
	}
      }
      else if(argv[i][1] == 's') {
        if(sscanf(&(argv[i][2]), "%lf", &scale) != 1) cmderr = TRUE;
	else if(scale <= 0.0) cmderr = TRUE;
	if(cmderr) {
	  fprintf(stderr, "%s: bad image scale factor '%s'\n",
		  argv[0], &argv[i][2]);
	  break;
	}
      }
      else if(argv[i][1] == 'w') {
        if(sscanf(&(argv[i][2]), "%lf", &linewd) != 1) {
				/* no check for < 0 so dash (-1) is pos. */
	  fprintf(stderr, "%s: bad line width '%s'\n",
		  argv[0], &argv[i][2]);
	  cmderr = TRUE;
	  break;
	}
      }
      /* -x sets up axes of default length, -x<len> uses len as length */
      else if(argv[i][1] == 'x') {
	if(argv[i][2] == '\0') x_ = TRUE;
	else {
	  if(sscanf(&(argv[i][2]), "%lf", &axeslen) != 1) {
				/* no check for < 0 so axes can flip */
	    fprintf(stderr, "%s: bad axes length '%s'\n",
		    argv[0], &argv[i][2]);
	    cmderr = TRUE;
	    break;
	  }
	  else x_ = TRUE;
	}
      }
      else if(argv[i][1] == 'v') s_ = TRUE;
      else if(argv[i][1] == 'n') n_ = TRUE;
      else if(argv[i][1] == 'g') g_ = TRUE;
      else if(argv[i][1] == 'c') c_ = TRUE;
      else if(argv[i][1] == 'm') m_ = TRUE;
      else if(argv[i][1] == 'q') {
	get_ps_file_base(argv, argc); /* set up the output file base */
	qpic_name_list = &(argv[i][2]);
	q_ = TRUE;
      }
      else if(argv[i][1] == 'u') {
	if(!strcmp(&(argv[i][2]), "x") || !strcmp(&(argv[i][2]), "X"))
	    up_axis = XI;
	else if(!strcmp(&(argv[i][2]), "y") || !strcmp(&(argv[i][2]), "Y"))
	    up_axis = YI;
	else if(!strcmp(&(argv[i][2]), "z") || !strcmp(&(argv[i][2]), "Z"))
	    up_axis = ZI;
	else {
	  fprintf(stderr, "%s: bad up axis type `%s' -- use x, y or z\n");
          cmderr = TRUE;
          break;
        }
      }
/*#endif*/
      else {
	fprintf(stderr, "%s: illegal option -- %s\n", argv[0], &(argv[i][1]));
	cmderr = TRUE;
	break;
      }
    }
    else {			/* isn't an option, must be the input file */
      *input_file = argv[i];
    }
  }

  if(cmderr == TRUE) {
/*#if CAPVEW == ON*/
    fprintf(stderr,
	    "Usage: '%s [<input file>]\n                [-rs<cond list>] [-ri<cond list>]\n                [-] [-l<list file>] [-a<azimuth>] [-e<elevation>]\n                [-r<rotation>] [-h<distance>] [-s<scale>] [-w<linewidth>]\n                [-u<upaxis>] [-q<cond list>] [-rc<cond list>] [-x<axeslength>]\n                [-b<.figfile>] [-m] [-rk] [-rd] [-dc] [-c] [-v] [-n] [-f] [-g]\n", argv[0]);
    fprintf(stderr, "DEFAULT VALUES:\n");
    fprintf(stderr, "  azimuth = %g\n  elevation = %g\n  rotation = %g\n",
	    DEFAZM, DEFELE, DEFROT);
    fprintf(stderr,
	    "  distance = %g (0 => 1 object radius away from center)\n",
	    DEFDST);
    fprintf(stderr, "  scale = %g\n  linewidth = %g\n",
	    DEFDST, DEFSCL, DEFWID);
    if(DEFUAX == XI) fprintf(stderr, "  upaxis = x\n");
    else if(DEFUAX == YI) fprintf(stderr, "  upaxis = y\n");
    else if(DEFUAX == ZI) fprintf(stderr, "  upaxis = z\n");
    fprintf(stderr, "  axeslength = %g\n", DEFAXE);
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "  -   = force conductor surface file read from stdin\n");
    fprintf(stderr, "  -m  = DUMP in MATLAB FORMAT instead of postscript\n");
    fprintf(stderr, "  -rs = remove conductors from solve list\n");
    fprintf(stderr, "  -ri = remove conductors from input\n");
    fprintf(stderr,
     "  -q  = select conductors for at-1V charge distribution .ps pictures\n");
    fprintf(stderr,
     "  -rc = remove conductors from all charge distribution .ps pictures\n");
    fprintf(stderr,
"  -b  = superimpose lines, arrows and dots in .figfile on all .ps pictures\n");
    fprintf(stderr,
      "  -rk = remove key in shaded .ps picture file (use with -q option)\n");
    fprintf(stderr,
      "  -rd = remove DIELEC type surfaces from all .ps picture files\n");
    fprintf(stderr,
"  -dc = display total charges in shaded .ps picture file (use with -q option)\n");
    fprintf(stderr, "  -c  = print command line in .ps picture file\n");
    fprintf(stderr, "  -v  = suppress showpage in all .ps picture files\n");
    fprintf(stderr, "  -n  = number faces with input order numbers\n");
    fprintf(stderr, "  -f  = do not fill in faces (don't rmv hidden lines)\n");
    fprintf(stderr, "  -g  = dump depth graph and quit\n");
    fprintf(stderr, "  <cond list> = [<name>],[<name>],...,[<name>]\n");
    /*dumpConfig(stderr, argv[0]);*/
    exit(0);
  }
}

/*
  surface information input routine - panels are read by read_panels()
*/
surface *read_all_surfaces(input_file, surf_list_file, read_from_stdin, infile,
			   relperm)
int read_from_stdin;
char *input_file, *surf_list_file, *infile;
double relperm;
{
  int num_surf, i;
  char group_name[BUFSIZ];
  surface *surf_list, *cur_surf;

  /* get the surfaces from stdin, the list file or the file on cmd line */
  /* the `- ' option always forces the first cond surf read from stdin */
  /* can also read from stdin if there's no list file and no cmd line file */
  infile[0] = '\0';
  num_surf = 0;
  surf_list = NULL;
  strcpy(group_name, "GROUP1");
  if(read_from_stdin || (input_file == NULL && surf_list_file == NULL)) {
    CALLOC(surf_list, 1, surface, ON, AMSC);
    surf_list->type = CONDTR;	/* only conductors can come in stdin */
    CALLOC(surf_list->name, strlen("stdin")+1, char, ON, AMSC);
    strcpy(surf_list->name, "stdin");
    surf_list->outer_perm = relperm;
    surf_list->end_of_chain = TRUE;

    /* set up group name */
    CALLOC(surf_list->group_name, strlen(group_name)+1, char, ON, AMSC);
    strcpy(surf_list->group_name, group_name);
    strcpy(group_name, "GROUP2");

    cur_surf = surf_list;

    strcpy(infile, "stdin");
    num_surf++;
  }

  /* set up to read from command line file, if necessary */
  if(input_file != NULL) {
    if(surf_list == NULL) {
      CALLOC(surf_list, 1, surface, ON, AMSC);
      cur_surf = surf_list;
    }
    else {
      CALLOC(cur_surf->next, 1, surface, ON, AMSC);
      cur_surf = cur_surf->next;
    }
    cur_surf->type = CONDTR;
    CALLOC(cur_surf->name, strlen(input_file)+1, char, ON, AMSC);
    strcpy(cur_surf->name, input_file);
    cur_surf->outer_perm = relperm;
    cur_surf->end_of_chain = TRUE;

    /* set up group name */
    CALLOC(cur_surf->group_name, strlen(group_name)+1, char, ON, AMSC);
    strcpy(cur_surf->group_name, group_name);

    for(i = 0; infile[i] != '\0'; i++);
    if(infile[0] != '\0') sprintf(&(infile[i]), ", %s", input_file);
    else sprintf(&(infile[i]), "%s", input_file);
    num_surf++;
    read_from_stdin++;
  }

  /* read list file if present */
  if(surf_list_file != NULL) {
    read_list_file(&surf_list, &num_surf, surf_list_file, read_from_stdin);
    for(i = 0; infile[i] != '\0'; i++);
    if(infile[0] != '\0') sprintf(&(infile[i]), ", %s", surf_list_file);
    else sprintf(&(infile[i]), "%s", surf_list_file);
  }

  return(surf_list);
}

/*
  surface input routine and command line parser
  - inputs surfaces (ie file names whose panels are read in read_panels)
  - sets parameters accordingly
*/
surface *input_surfaces(argv, argc, autmom, autlev, relperm,
			numMom, numLev, infile)
int argc, *autmom, *autlev, *numMom, *numLev;
double *relperm;
char *argv[], *infile;
{
  int read_from_stdin, num_surf;
  surface *read_all_surfaces();
  char *surf_list_file, *input_file;

  /* initialize defaults */
  surf_list_file = input_file = NULL;
  read_from_stdin = FALSE;

  parse_command_line(argv, argc, autmom, autlev, relperm, numMom, numLev,
		     &input_file, &surf_list_file, &read_from_stdin);

  return(read_all_surfaces(input_file, surf_list_file,
			   read_from_stdin, infile, *relperm));
}

/*
  dump the data associated with the input surfaces
*/
void dumpSurfDat(surf_list)
surface *surf_list;
{
  surface *cur_surf;
  char *hack_path();

  fprintf(stdout, "  Input surfaces:\n");
  for(cur_surf = surf_list; cur_surf != NULL; cur_surf = cur_surf->next) {

    /* possibly write group name */
    if(cur_surf == surf_list) fprintf(stdout, "   %s\n", cur_surf->group_name);
    else if(cur_surf->prev->end_of_chain)
	fprintf(stdout, "   %s\n", cur_surf->group_name);

    /* write file name */
    fprintf(stdout, "    %s", hack_path(cur_surf->name));
    if(cur_surf->type == CONDTR) {
      fprintf(stdout, ", conductor\n");
      fprintf(stdout, "      title: `%s'\n", cur_surf->title);
      fprintf(stdout, "      outer permittivity: %g\n",
	      cur_surf->outer_perm);
    }
    else if(cur_surf->type == DIELEC) {
      fprintf(stdout, ", dielectric interface\n");
      fprintf(stdout, "      title: `%s'\n", cur_surf->title);
      fprintf(stdout, "      permittivities: %g (inner) %g (outer)\n",
	      cur_surf->inner_perm, cur_surf->outer_perm);
    }
    else if(cur_surf->type == BOTH) {
      fprintf(stdout, ", thin conductor on dielectric interface\n");
      fprintf(stdout, "      title: `%s'\n", cur_surf->title);
      fprintf(stdout, "      permittivities: %g (inner) %g (outer)\n",
	      cur_surf->inner_perm, cur_surf->outer_perm);
    }
    else {
      fprintf(stderr, "dumpSurfDat: bad surface type\n");
      exit(0);
    }
    fprintf(stdout,"      number of panels: %d\n",
	    cur_surf->num_panels - cur_surf->num_dummies);
    fprintf(stdout,"      number of extra evaluation points: %d\n",
	    cur_surf->num_dummies);
    fprintf(stdout,"      translation: (%g %g %g)\n",
	    cur_surf->trans[0], cur_surf->trans[1], cur_surf->trans[2]);

  }
}

/*
  replaces name (and all aliases) corresponding to "num" with unique string
*/
void remove_name(name_list, num)
int num;
Name **name_list;
{
  static char str[] = "%`_^#$REMOVED";
  Name *cur_name, *cur_alias;
  int i, slen;

  slen = strlen(str);

  for(i = 1, cur_name = *name_list; cur_name != NULL;
      cur_name = cur_name->next, i++) {
    if(i == num) {

      /* overwrite name */
      if(strlen(cur_name->name) < slen) {
	CALLOC(cur_name->name, slen+1, char, ON, AMSC);
      }
      strcpy(cur_name->name, str);

      /* overwrite aliases */
      for(cur_alias = cur_name->alias_list; cur_alias != NULL;
	  cur_alias = cur_alias->next) {
	if(strlen(cur_alias->name) < slen) {
	  CALLOC(cur_alias->name, slen+1, char, ON, AMSC);
	}
	strcpy(cur_alias->name, str);
      }
    }
  }

}

/*
  removes (unlinks from linked list) panels that are on conductors to delete
*/
void remove_conds(panels, num_list, name_list)
charge **panels;
ITER *num_list;
Name **name_list;
{
  ITER *cur_num;
  charge *cur_panel, *prev_panel;

  for(cur_panel = prev_panel = *panels; cur_panel != NULL;
      cur_panel = cur_panel->next) {
    if(cur_panel->dummy) continue;
    if(cur_panel->surf->type == CONDTR || cur_panel->surf->type == BOTH) {
      if(want_this_iter(num_list, cur_panel->cond)) {
	/* panel's conductor is to be removed, so unlink the panel */
	/* - if panel to be removed is first panel, rewrite head pointer */
	if(cur_panel == *panels) *panels = cur_panel->next;
	/* - otherwise bypass cur_panel with next pointers */
	else prev_panel->next = cur_panel->next;
      }
      else prev_panel = cur_panel;
    }
  }

  /* remove all -ri'd conductor names from master name list
     - required to get rid of references in capsolve()
     - actually, name and all its aliases are replaced by ugly string
       (not the cleanest thing) */
  for(cur_num = num_list; cur_num != NULL; cur_num = cur_num->next) {
    remove_name(name_list, cur_num->iter);
  }
}

/*
  checks for kill lists with inconsistent demands
  -rs list: can't remove a conductor physically removed from computation w/-ri
  -q list: can't dump q plot for cond physically rmed or rmed from comp
  -rc list: no restrictions
  -ri/-rs: can't exhaust all conductors with combination of these lists
*/
void resolve_kill_lists(rs_num_list, q_num_list, ri_num_list, num_cond)
ITER *rs_num_list, *q_num_list, *ri_num_list;
int num_cond;
{
  int i, lists_exhaustive;
  ITER *cur_num;
  extern int m_;

  /* check for anything in -rs list in -ri list */
  for(cur_num = ri_num_list; cur_num != NULL; cur_num = cur_num->next) {
    if(want_this_iter(rs_num_list, cur_num->iter)) {
      fprintf(stderr,
 "resolve_kill_lists: a conductor removed with -ri is in the -rs list\n");
      exit(0);
    }
  }

  /* check for anything in -q list in -ri or -rs list
     - recall that -q by itself means plot for all active,
       so null q_num_list always ok */
  for(cur_num = q_num_list; cur_num != NULL; cur_num = cur_num->next) {
    if(want_this_iter(rs_num_list, cur_num->iter)
       || want_this_iter(ri_num_list, cur_num->iter)) {
      fprintf(stderr,
"resolve_kill_lists: a conductor removed with -ri or -rs is in the -q list\n");
      exit(0);
    }
  }

  /* check that -rs and -ri lists don't exhaust all conductors */
  lists_exhaustive = TRUE;
  for(i = 1; i <= num_cond; i++) {
    if(!want_this_iter(rs_num_list, i) && !want_this_iter(ri_num_list, i)) {
      lists_exhaustive = FALSE;
      break;
    }
  }
  if(lists_exhaustive && !m_) {
    fprintf(stderr,
"resolve_kill_lists: all conductors either in -ri or -rs list\n");
    exit(0);
  }
}

/*
  main input routine, returns a list of panels in the problem
*/
charge *input_problem(argv, argc, autmom, autlev, relperm,
		      numMom, numLev, name_list, num_cond)
int argc, *autmom, *autlev, *numMom, *numLev, *num_cond;
double *relperm;
char *argv[];
Name **name_list;
{
  surface *surf_list, *input_surfaces();
  char infile[BUFSIZ], *ctime(), hostname[BUFSIZ];
  charge *read_panels(), *chglist;
  long clock;
  extern ITER *kill_num_list, *qpic_num_list, *kinp_num_list, *kq_num_list;
  extern char *kill_name_list, *qpic_name_list, *kinp_name_list;
  extern char *kq_name_list;

  /* read the conductor and dielectric interface surface files, parse cmds */
  surf_list = input_surfaces(argv, argc, autmom, autlev, relperm,
			     numMom, numLev, infile);

  if(*autmom == ON) *numMom = DEFORD;

#if DIRSOL == ON || EXPGCR == ON
  /*fprintf(stderr, "DIRSOL and EXPGCR compile options not implemented\n");
  exit(0);*/
  *numLev = 0;	       	/* put all the charges in first cube */
  *autlev = OFF;
#endif

  strcpy(hostname, "19Sep96");
  fprintf(stdout, "Running %s %.1f (%s)\n  Input: %s\n",
	  argv[0], VERSION, hostname, infile);

  /* input the panels from the surface files */
  *num_cond = 0;		/* initialize conductor count */
  chglist = read_panels(surf_list, name_list, num_cond);

#if 1 == 0
 removed for fasthenry: Not used for anything with visualization, hopefully
  /* set up the lists of conductors to remove from solve list */
  kill_num_list = get_kill_num_list(*name_list, kill_name_list);

  /* remove the panels on specified conductors from input list */
  kinp_num_list = get_kill_num_list(*name_list, kinp_name_list);
  remove_conds(&chglist, kinp_num_list, name_list);

  /* set up the lists of conductors to dump shaded plots for */
  qpic_num_list = get_kill_num_list(*name_list, qpic_name_list);

  /* set up the lists of conductors to eliminate from shaded plots */
  kq_num_list = get_kill_num_list(*name_list, kq_name_list);

  /* check for inconsistencies in kill lists */
  resolve_kill_lists(kill_num_list, qpic_num_list, kinp_num_list, *num_cond);
#endif

#if DISSRF == ON
  dumpSurfDat(surf_list);
#endif

  time(&clock);
  fprintf(stdout, "  Date: %s", ctime(&clock));
#ifndef NO_GETHOSTNAME
#ifndef SOLARIS
  if(gethostname(hostname, BUFSIZ) != -1)
      fprintf(stdout, "  Host: %s\n", hostname);
  else fprintf(stdout, "  Host: ? (gethostname() failure)\n");
#else
  if (sysinfo(SI_HOSTNAME,hostname,BUFSIZ) != -1)
      fprintf(stdout, "  Host: %s\n", hostname);
  else fprintf(stdout, "  Host: ? (sysinfo() failure)\n");
#endif
#endif

#if CFGDAT == ON
  dumpConfig(stdout, argv[0]);
#endif

  /* return the panels from the surface files */
  return(chglist);
}
