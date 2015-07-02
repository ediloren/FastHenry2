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

// Enrico
#include <string.h>

main(argc, argv)
int argc;
char *argv[];
{
  int ttliter, i, j, num_cond;
  charge *chglist, *nq, *input_problem();
  double *get_q();
  ssystem *sys, *mulInit();
  double **capmat, dirtimesav, mulsetup, initalltime, ttlsetup, ttlsolve;
  double relperm;
  int autmom, autlev, numMom, numLev;
  char *fname, *concat3();

  extern int fulldirops, fullPqops;
  extern int num_dummy_panels, num_dielec_panels;
  extern int num_both_panels, num_cond_panels, up_size, eval_size;
  extern char *title, *ps_file_base, *in_file_name;
  extern long memcount;
  extern double prectime, conjtime, dirtime, multime, uptime, downtime;
  extern double evaltime, lutime, fullsoltime, prsetime;
  extern char *kill_name_list;

  Name *name_list;

#if DUMPPS == ON || DUMPPS == ALL
  char filename[BUFSIZ];
#endif

/*#if CAPVEW == ON*/
  extern char **argvals;
  extern int argcnt, m_, q_, dd_;
  extern double ***axes;
/*#endif*/

  /* initialize memory and time counters, etc. */
  fulldirops = fullPqops = 0;
  prectime = conjtime = dirtime = multime = uptime = downtime = 0.0;
  evaltime = lutime = fullsoltime = mulsetup = 0.0;
  memcount = 0;
  CALLOC(title, BUFSIZ, char, ON, AMSC);

  /* initialize defaults, etc */
  autmom = autlev = ON;
  relperm = 1.0;
  argvals = argv;
  argcnt = argc;
  CALLOC(axes, 10, double **, ON, AMSC);
  for(i = 0; i < 10; i++) {
    CALLOC(axes[i], 2, double *, ON, AMSC);
    for(j = 0; j < 2; j++) {
      CALLOC(axes[i][j], 3, double, ON, AMSC);
    }
  }

  /* get the list of all panels in the problem */
  /* - many command line parameters having to do with the postscript
       file dumping interface are passed back via globals (see mulGlobal.c) */
  chglist = input_problem(argv, argc, &autmom, &autlev, &relperm,
			  &numMom, &numLev, &name_list, &num_cond);

  /* just dump the psfile */
    get_ps_file_base(argv, argc);
    if (!q_) {
      if (m_)
        /* dump in matlab format */
        dump_struct(chglist, NULL);
      else {
        /* dump a postscript file */
        fprintf(stdout,"Creating postscript file\n");
        dump_ps_geometry(chglist, NULL, 0, dd_);
      }
    }
    else {
      fname = concat3(in_file_name,"_","shadings");
      if (m_)
        /* dump in matlab format */
        dump_struct(chglist, get_q(fname,chglist));
      else {
        /* dump in postscript format */
        fprintf(stdout,"Creating postscript file with shading\n");
        dump_ps_geometry(chglist, get_q(fname,chglist), 0, dd_);
      }
    }
    exit(0);

}

char *concat3(st1, st2, st3)
char *st1, *st2, *st3;
{
  int length = 0;
  char *allthree;

  length = strlen(st1);
  length += strlen(st2);
  length += strlen(st3);

  CALLOC(allthree, length+1, char, ON, AMSC);

  allthree[0] = '\0';
  strcat(allthree,st1);
  strcat(allthree,st2);
  strcat(allthree,st3);

  return allthree;
}

/*  For FastHenry, this reads the file of shadings for each of the faces */

double *get_q(fname, chglist)
char *fname;
charge *chglist;
{
  FILE *fp;
  int numchgs = 0;
  charge *chg;
  int error, i;
  double *q;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("Couldn't open %s for -q option\n", fname);
    exit(1);
  }

  for(chg = chglist; chg != NULL; chg = chg->next)
    numchgs++;

  CALLOC(q, numchgs, double, ON, AMSC);
  error = 0;

  for(i = 0; i < numchgs && error == 0; i++) {
    if (error == 0) {
      if (fscanf(fp, "%lf", &q[i]) != 1) {
        error = 1;
        fprintf(stderr, "Error reading shading file. Rest of panels: q = 0\n");
      }
    }
    else
      q[i] = 0.0;
  }

  return q;

}
