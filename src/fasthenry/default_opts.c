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
 
*//* This fills an ind_opts with the default options */
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

