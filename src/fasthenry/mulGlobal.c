/* # ***** sort to /src/main
   # ***** */
#include "mulGlobal.h"

long memcount;	       	/* allocated memory counter */
long memQ2M;			/* allocated memory counters by function */
long memQ2L;
long memQ2P;
long memL2L;
long memM2M;
long memM2L;
long memM2P;
long memL2P;
long memQ2PD;
long memMSC;
long memIND;

#ifdef MATTDEBUG
long membins[1001];
#endif

/* 
  global timer and operation count accumulators
*/
double prectime;	       	/* time spent doing back solve for prec */
double prsetime;		/* time spent calculating preconditioner */
double conjtime;	        /* time spent doing everything but A*q */
double dirtime;			/* time for direct part of P*q */
double multime;			/* time for multipole part of P*q */
double uptime;			/* time in mulUp(), upward pass */
double downtime;		/* time in mulDown(), downward pass */
double evaltime;		/* time in mulEval(), evaluation pass */
int fulldirops;			/* total direct operations - DIRSOL=ON only */
double lutime;			/* factorization time DIRSOL=ON only */
double fullsoltime;		/* time for solves, DIRSOL=ON only */
int fullPqops;			/* total P*q ops using P on disk - EXPGCR=ON */

/*
  misc global
*/
NAME *start_name = NULL;	/* conductor name linked list head */
NAME *current_name;		/* conductor name linked list tail */
NAME *start_name_this_time;	/* cond name list for the current surface */
char *kill_name_list;		/* cond names whose columns are omitted */
ITER *kill_num_list;		/* cond numbers whose columns are omitted */
char *kinp_name_list;		/* cond names omitted from input */
ITER *kinp_num_list;		/* cond numbers omitted from input */
char *qpic_name_list;		/* cond column names that get q picture */
ITER *qpic_num_list;		/* cond column names that get q picture */
char *kq_name_list;		/* cond names removed from q picture */
ITER *kq_num_list;		/* cond numbers removed from q picture */

int num_dielec_panels;		/* number of dielectric interface panels */
int num_both_panels;		/* number of thin-cond-on-dielec-i/f panels */
int num_cond_panels;		/* number of thick conductor panels */
int up_size;			/* sum of above three (real panels) */
int num_dummy_panels;		/* number of off-panel eval pnt panels */
int eval_size;			/* sum of above two (total panel structs) */
double iter_tol;		/* iterative loop tolerence on ||r|| */

/*
  command line option variables - all have to do with ps file dumping
*/
#if CAPVEW == ON		/* eliminate messy globals if not needed */
char **argvals;			/* copy of argv */
int argcnt;			/* copy of argc */
int s_;				/* TRUE => insert showpage in .ps file(s) */
int n_;				/* TRUE => number faces with input ordering */
int g_;				/* TRUE => dump depth graph and quit */
int c_;				/* TRUE => print command line on .ps file(s) */
int x_;				/* TRUE => axes have been specified */
int k_;
int rc_;			/* TRUE => rm conductors in list from pic */
int rd_;			/* TRUE => remove all dielec i/fs from pic */
int rb_;			/* TRUE => rm BOTH-types in list from pic */
int q_;				/* TRUE => dump shaded plots of q_iter iters */
int rk_;			/* TRUE => rm chg den key in -q plots */
int m_;				/* TRUE => switch to plot gen mode */
int f_;				/* TRUE => don't fill faces (no hidden l rm) */
int dd_;			/* TRUE => dump ttl charges to .ps pictures */
double view[3];			/* absolute view point of 3D geometry */
double moffset[2];		/* image offset from lower left corner */
double elevation;		/* elevation of view rel to center of object */
double azimuth;			/* azimuth of view rel to center of object */
double rotation;		/* image rotation, degrees */
double distance;		/* relative distance from center (#radii-1) */
double linewd;			/* postscript line width */
double scale;			/* over all image scale factor */
double axeslen;			/* axes lengths in 3D distance */
double axes[10][2][3];		/* the 2d image of the coordinate axes */
int up_axis;			/* X,Y or Z => which axis is vertical in pic */
char *line_file;		/* pointer to .fig superimposed line file */
char *ps_file_base;		/* pointer to base name for .ps files */
#endif
