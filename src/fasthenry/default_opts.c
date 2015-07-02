/* This fills an ind_opts with the default options */
/*  There is little error checking that these are valid */
#include "induct.h"

default_opts(opts)
ind_opts *opts;
{
  opts->soln_technique = ITERATIVE;      /* -s */
  opts->mat_vect_prod  = MULTIPOLE;      /* -m */
  opts->precond = ON;                    /* -p */
  opts->order = 2;                       /* -o */
  opts->level = AUTO;                    /* -l */
  opts->makeFastCapFile = OFF;           /* -f */
  opts->gp_draw = OFF;                   /* -g */
  opts->auto_refine = ON;                /* -a */
  opts->init_refine = 0;                 /* -i */ 
  opts->dumpMats = OFF;                  /* -d */
  opts->orderROM = -1;                   /* -r */
  opts->onlyROM = 0;                     /* -M */
  opts->kind = MATLAB;                   /* -k */
  opts->tol = 1e-3;                      /* -t */
  opts->abs_tol = 1e-2;                  /* -b */
  opts->maxiters = 200;                  /* -c */ 
  opts->limit = AUTO;                    /* -e */
  opts->debug = OFF;                     /* -D */
  opts->portlist = NULL;                 /* -x */
  opts->suffix = "";                     /* -S */
  opts->shell_r0 = 0.87;                 /* -R */
  opts->regurgitate = FALSE;             /* -v */
  opts->fname = NULL;             
}

