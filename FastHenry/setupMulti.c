/* this is 'fastcap.c' */

#include "FHWindow.h"

#include "induct.h"

ssystem *SetupMulti(chglist, indsys)
charge *chglist;
SYS *indsys;
{
  int i, j;
  charge *nq, *input_problem();
  ssystem *sys, *mulInit();
  double dirtimesav, mulsetup, initalltime, ttlsetup, ttlsolve;
  double relperm;
  int autmom, autlev, numMom, numLev;

  int pseudo_lev;   /* level for breaking up only */

  extern int fulldirops, fullPqops;
  extern int num_dummy_panels, num_dielec_panels; 
  extern int num_both_panels, num_cond_panels, up_size, eval_size;
  extern char *title, *ps_file_base;
  extern long long memcount;
  extern double prectime, conjtime, dirtime, multime, uptime, downtime;
  extern double evaltime, lutime, fullsoltime, prsetime;
  extern char *kill_name_list;

#if DUMPPS == ON || DUMPPS == ALL
  char filename[BUFSIZ];
#endif

#if CAPVEW == ON
  extern char **argvals;
  extern int argcnt, m_, q_, dd_;
  extern double ***axes;
#endif

  /* initialize memory and time counters, etc. */
  fulldirops = fullPqops = 0;
  prectime = conjtime = dirtime = multime = uptime = downtime = 0.0;
  evaltime = lutime = fullsoltime = mulsetup = 0.0;
/*  memcount = 0;   Do this in induct.c */
/*  CALLOC(title, BUFSIZ, char, ON, AMSC); */

  /* initialize defaults, etc */
  autmom = autlev = ON;
  relperm = 1.0;

#if CAPVEW == ON
/*  argvals = argv;
  argcnt = argc;
*/
  CALLOC(axes, 10, double **, ON, AMSC);
  for(i = 0; i < 10; i++) {
    CALLOC(axes[i], 2, double *, ON, AMSC);
    for(j = 0; j < 2; j++) {
      CALLOC(axes[i][j], 3, double, ON, AMSC);
    }
  }
#endif

  /* get the list of all panels in the problem */
  /* - many command line parameters having to do with the postscript
     file dumping interface are passed back via globals (see mulGlobal.c) */
  /* Name *name_list; 
     int num_cond;
    chglist = input_problem(argv, argc, &autmom, &autlev, &relperm, 
                              &numMom, &numLev, &name_list, &num_cond);
*/
  if(autmom == ON) numMom = DEFORD;
  else numMom = DEFORD;
  
#if DIRSOL == ON || EXPGCR == ON
  /*viewprintf(stderr, "DIRSOL and EXPGCR compile options not implemented\n");
    FHExit(FH_NORMAL_END);*/
  *numLev = 0;	       	/* put all the charges in first cube */
  *autlev = OFF;
#endif

#if CAPVEW == ON
  /* if no fastcap run is to be done, just dump the psfile */
/*
  if(m_) {
    if(!q_) get_ps_file_base(argv, argc);
    dump_ps_geometry(chglist, NULL, 0, dd_);
    FHExit(FH_NORMAL_END);
  }
*/
#endif

  /* the following is added for FastHenry to alter settings for mulInit()*/

  autlev = OFF;
  numMom = indsys->opts->order;
  if (indsys->opts->level == AUTO) 
    autlev = ON;
  else
    numLev = indsys->opts->level;

/*
  if (numLev == 0) {
    printf("Level to break up segments:");
    scanf("%d",&pseudo_lev);
  }
  else 
*/
     pseudo_lev = numLev;

  /*  starttimer; */
  sys = mulInit(autlev, numLev, numMom, chglist, indsys, pseudo_lev); 
                         /*Set up cubes, charges.*/
  /* stoptimer */
  initalltime = dtime;
 
  numLev = sys->depth;

/*
  sys->num_cond = num_cond;
  sys->cond_names = name_list;
  */

  viewprintf(stdout, "\nMultipole Summary\n");

#if CMDDAT == ON
  viewprintf(stdout, "  Expansion order: %d\n", numMom);
  if (sys->depth == sys->real_depth)
    viewprintf(stdout, "  Number of partitioning levels: %d\n", numLev);
  else {
    viewprintf(stdout, "  Number of partitioning levels: %d\n", sys->real_depth);
    viewprintf(stdout, "  Number of levels for Preconditioner: %d\n", numLev);
  }
/*  viewprintf(stdout, "  Overall permittivity factor: %.3g\n", relperm); */
#endif

  /* Figure out number of panels and conductors. */
  eval_size = up_size = num_dummy_panels = num_dielec_panels = 0;
  num_both_panels = num_cond_panels = 0;
  for(nq = chglist; nq != NULL; nq = nq->next) {
    if(nq->dummy) num_dummy_panels++;
    else if(nq->surf->type == CONDTR) {
      num_cond_panels++;
    }
    else if(nq->surf->type == DIELEC) num_dielec_panels++;
    else if(nq->surf->type == BOTH) num_both_panels++;
  }
  up_size = num_cond_panels + num_both_panels + num_dielec_panels;
  eval_size = up_size + num_dummy_panels;

#if DISSRF == OFF
/*  viewprintf(stdout, "Title: `%s'\n", title); */
#endif
  viewprintf(stdout, "  Total number of filaments: %d\n", up_size);
#if 1==0
  viewprintf(stdout, "    Number of conductor panels: %d\n", num_cond_panels);
  viewprintf(stdout, "    Number of dielectric interface panels: %d\n", 
	  num_dielec_panels);
  viewprintf(stdout, 
	  "    Number of thin conductor on dielectric interface panels: %d\n", 
	  num_both_panels);
  /*viewprintf(stdout, "  Number of extra evaluation points: %d\n", 
	  num_dummy_panels);*/
  viewprintf(stdout, "  Number of conductors: %d\n", num_cond);
#endif

  /*********************************************************************/
  /* exit now if we only want to do multipole setup and never evaluate */
  if (indsys->opts->mat_vect_prod == DIRECT) {
    /* viewprintf(stdout, "  Inductance mem allocated: %7.d kilobytes\n",
	    memIND/1024); */
    return sys;
  }
  /*********************************************************************/

#if NAMDAT == ON
  dumpCondNames(stdout, name_list);
#endif

  if(num_both_panels > 0) {
    viewprintf(stderr, 
	    "Thin cond panels on dielectric interfaces not supported\n");
    FHExit(FH_NORMAL_END);
  }

#if CKCLST == ON
  viewprintf(stdout, "Checking panels...");
  if(has_duplicate_panels(stdout, chglist)) {
    viewprintf(stdout, "charge list has duplicates\n");
    FHExit(FH_GENERIC_ERROR);
  }
  viewprintf(stdout, "no duplicates\n");
#endif

#if MULDAT == ON
  dumpMulSet(sys, numLev, numMom);
#endif
  fflush(stdout);

  /* starttimer; */
  mulMultiAlloc(MAX(sys->max_eval_pnt, sys->max_panel), numMom, sys->depth);
  /* stoptimer */
  initalltime += dtime;		/* save initial allocation time */

#if DUMPPS == ON || DUMPPS == ALL
  strcpy(filename, "psmat.ps");
  dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, OPEN);
#endif

  mulMatDirect(sys);		/* Compute the direct part matrices. */

  if (bFHContinue == FALSE) {		// Enrico 
	return sys;
  }

#if DIRSOL == OFF		/* with DIRSOL just want to skip to solve */

#if PRECOND == BD
  /* starttimer; */
  bdmulMatPrecond(sys);
  /* stoptimer */
  prsetime = dtime;		/* preconditioner set up time */
#endif

#if PRECOND == OL
  /* starttimer; */
  olmulMatPrecond(sys);
  /* stoptimer */
  prsetime = dtime;		/* preconditioner set up time */
#endif

#if DMPREC == ON
  dump_preconditioner(sys, chglist, 1);	/* dump prec. and P to matlab file */
#endif

#if DPSYSD == ON
  dissys(sys);
#endif

#if CKDLST == ON
  chkList(sys, DIRECT);
#endif

#endif				/* DIRSOL == OFF */
  //dumpnums(ON, eval_size, up_size); /* save num/type of pot. coeff calcs */
  dumpnums(ON, eval_size); /* save num/type of pot. coeff calcs */

  dirtimesav = dirtime;		/* save direct matrix setup time */
  dirtime = 0.0;		/* make way for direct solve time */

#if DIRSOL == OFF

#if DUMPPS == ON
  dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, CLOSE);
#endif

  /* starttimer; */
  mulMatUp(sys);		/* Compute the upward pass matrices. */

#if DNTYPE == NOSHFT
  mulMatDown(sys);		/* find matrices for no L2L shift dwnwd pass */
#endif

#if DNTYPE == GRENGD
  mulMatDown(sys);		/* find matrices for full Greengard dnwd pass*/
#endif

#if CKDLST == ON
  chkList(sys, DIRECT);
  chkLowLev(sys, DIRECT);
  chkEvalLstD(sys, DIRECT);
#endif

  mulMatEval(sys);		/* set up matrices for evaluation pass */

  if (bFHContinue == FALSE) {		// Enrico 
	return sys;
  }
  /* stoptimer */
  mulsetup = dtime;		/* save multipole matrix setup time */

  //dumpnums(OFF, eval_size, up_size);	/* dump num/type of pot. coeff calcs */
  dumpnums(OFF, eval_size);	/* dump num/type of pot. coeff calcs */

#if DUMPPS == ALL
  dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, CLOSE);
#endif

#if DISSYN == ON
  dumpSynop(sys);
#endif

#if DMTCNT == ON
  dumpMatBldCnts(sys);
#endif

#endif				/* DIRSOL == ON */
/*
double **capmat;
int ttliter;
  viewprintf(stdout, "\nITERATION DATA");
  ttliter = capsolve(&capmat, sys, chglist, eval_size, up_size, num_cond,
		     name_list);
*/
#if MKSDAT == ON		/* dump symmetrized, 4 pi eps scaled matrix */
/*  mksCapDump(capmat, num_cond, relperm, &name_list); */
#endif

#if TIMDAT == ON 
  ttlsetup = initalltime + dirtimesav + mulsetup;
  multime = uptime + downtime + evaltime;
  ttlsolve = dirtime + multime + prectime + conjtime;

  if (indsys->opts->debug == ON) 
    viewprintf(stdout, "\nTIME AND MEMORY USAGE SYNOPSIS\n");
#endif

#ifdef OTHER
  if(TIMDAT == ON) {
    viewprintf(stdout, 
	    "Warning: compilation with OTHER flag gives incorrect times\n");
  }
#endif

#if TIMDAT == ON
#if 1==0
  viewprintf(stdout, "Total time: %g\n", ttlsetup + ttlsolve);
  viewprintf(stdout, "  Total setup time: %g\n", ttlsetup);
  viewprintf(stdout, "    Direct matrix setup time: %g\n", dirtimesav);
  viewprintf(stdout, "    Multipole matrix setup time: %g\n", mulsetup);
  viewprintf(stdout, "    Initial misc. allocation time: %g\n", initalltime);
  viewprintf(stdout, "  Total iterative P*q = psi solve time: %g\n", ttlsolve);
  viewprintf(stdout, "    P*q product time, direct part: %g\n", dirtime);
  viewprintf(stdout, "    Total P*q time, multipole part: %g\n", multime);
  viewprintf(stdout, "      Upward pass time: %g\n", uptime);
  viewprintf(stdout, "      Downward pass time: %g\n", downtime);
  viewprintf(stdout, "      Evaluation pass time: %g\n", evaltime);
  viewprintf(stdout, "    Preconditioner solution time: %g\n", prectime);
  viewprintf(stdout, "    Iterative loop overhead time: %g\n", conjtime);

  if(DIRSOL == ON) {		/* if solution is done by Gaussian elim. */
    viewprintf(stdout,"\nTotal direct, full matrix LU factor time: %g\n",lutime);
    viewprintf(stdout,"Total direct, full matrix solve time: %g\n",fullsoltime);
    viewprintf(stdout, "Total direct operations: %d\n", fulldirops);
  }
  else if(EXPGCR == ON) {	/* if solution done iteratively w/o multis */
    viewprintf(stdout,"\nTotal A*q operations: %d (%d/iter)\n", 
	    fullPqops, fullPqops/ttliter);
  }
#endif

 if (indsys->opts->debug == ON) {
  viewprintf(stdout, "Total memory allocated: %d kilobytes ", memcount/1024);
  uallocEfcy(memcount);

  viewprintf(stdout, "  Q2M  matrix memory allocated: %7.d kilobytes\n",
	  memQ2M/1024);
  memcount = memQ2M;
  viewprintf(stdout, "  Q2L  matrix memory allocated: %7.d kilobytes\n",
	  memQ2L/1024);
  memcount += memQ2L;
  viewprintf(stdout, "  Q2P  matrix memory allocated: %7.d kilobytes\n",
	  memQ2P/1024);
  memcount += memQ2P;
  viewprintf(stdout, "  L2L  matrix memory allocated: %7.d kilobytes\n",
	  memL2L/1024);
  memcount += memL2L;
  viewprintf(stdout, "  M2M  matrix memory allocated: %7.d kilobytes\n",
	  memM2M/1024);
  memcount += memM2M;
  viewprintf(stdout, "  M2L  matrix memory allocated: %7.d kilobytes\n",
	  memM2L/1024);
  memcount += memM2L;
  viewprintf(stdout, "  M2P  matrix memory allocated: %7.d kilobytes\n",
	  memM2P/1024);
  memcount += memM2P;
  viewprintf(stdout, "  L2P  matrix memory allocated: %7.d kilobytes\n",
	  memL2P/1024);
  memcount += memL2P;
  viewprintf(stdout, "  Q2PD matrix memory allocated: %7.d kilobytes\n",
	  memQ2PD/1024);
  memcount += memQ2PD;
  viewprintf(stdout, "  Miscellaneous mem. allocated: %7.d kilobytes\n",
	  memMSC/1024);
  memcount += memMSC;
  viewprintf(stdout, "  Inductance mem allocated: %7.d kilobytes\n",
	  memIND/1024);
  memcount += memIND;
  viewprintf(stdout, "  Total memory (check w/above): %7.d kilobytes\n",
	  memcount/1024);
 }
#endif

  return sys;
}
