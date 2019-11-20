/* This is the main part of the code */


#include "FHWindow.h" // Enrico
#include "induct.h"

#include "Sparse/spMatrix.h"

// include to use path manipulation routines, Enrico
#include <stdlib.h>
// include for _chdir()
#include <direct.h>
// include for isspace()
#include <ctype.h>

/* these are missing in some math.h files */
// Not needed anymore, included in standard math.h
//extern double asinh();
//extern double atanh();

#define MAXPRINT 1
#define MAXCHARS 400
#define TIMESIZE 10

static int notblankline();

// Enrico, static vars moved to file scope to be initialized
static int colmax = 0;
static double *tcol = NULL;
static CX **trows = NULL; 

FILE *fp, *fp2, *fp3, *fptemp, *fb, *fROM;
int num_exact_mutual;
int num_fourfil;
int num_mutualfil;
int num_found;
int num_perp;
int forced = 0;  /* for debugging inside exact_mutual() */

char outfname[200];
char outfname2[200];

/* savemat_mod machine type */
#ifdef DEC
int machine = 0000;
#elif defined WIN32		// Enrico
int machine = 0000;
#else
int machine = 1000;
#endif

// function prototypes
int pick_subset(strlist *portlist, SYS *indsys);
void formMtrans(SYS *indsys);
void savemats(SYS *indsys);
void formMZMt(SYS *indsys);
void dump_M_to_text(FILE *fp, MELEMENT **Mlist, int num_mesh, int nnz);
void dump_M_to_matlab(FILE *fp, MELEMENT **Mlist, int num_mesh, int nnz, char *mname);
void dumpMat_totextfile(FILE *fp, double **A, int rows, int cols);
void dumpVec_totextfile(FILE *fp2, double *Vec, int size);
void cx_dumpMat_totextfile(FILE *fp, CX **Z, int rows, int cols);
void saveCarray(FILE *fp, char *fname, double **Arr, int rows, int cols);
void dump_to_Ycond(FILE *fp, int cond, SYS *indsys);
void fillMrow(MELEMENT **Mlist, int mesh, double *Mrow);
void fillA(SYS *indsys);
void fillZ(SYS *indsys);
int nnz_inM(MELEMENT **Mlist, int num_mesh);
void pick_ground_nodes(SYS *indsys);
void savecmplx2(FILE *fp, char *name, CX **Z, int rows, int cols);


int fastHenryMain(argc, argv)
int argc;
char *argv[];
{

  double freq;
  GROUNDPLANE *plane;                   /* CMS 7/7/92 */
  SEGMENT *seg;
  NODES *node;
  int i,j,k,m, last, err;
  char fname[80], tempstr[10];
  CX *vect, **pvect;         /* space needed by gmres */
  double tol = 1e-8;
  double ftimes[TIMESIZE];
  int maxiters, user_maxiters;
  CX *b, *x0;
  SYS *indsys;           /* holds all the big variables (inductance system) */
  CX **MtZM;
  int num_mesh, num_extern, num_fils,num_nodes, num_real_nodes, num_segs;
  int num_sub_extern;
  double totaltime;
  int dont_form_Z;       /* if fmin=0, don't from the L matrix */
  char *MRMt;            /* may be needed if ROM is requested */
  double **B, **C;
  int actual_order;

  int num_planes, nonp, planemeshes, tree_meshes;           /* CMS 6/7/92 */
  void fillgrids();                               /* function in addgroundplane.c */

  int *meshsect;
  double fmin, fmax, logofstep;
  double *R;
  EXTERNAL *ext;
  ind_opts *opts;
  extern long long memcount;

  char full[_MAX_PATH], drive[_MAX_DRIVE], dir[_MAX_DIR];

  charge *chglist, *chgend, chgdummy;
  ssystem *sys, *SetupMulti();

  // Enrico
  int iindex, jindex, kindex;

  // initialize all global and static variables in various modules, Enrico
  InitGlobAndStatVars();

  memcount = 0;

  for(i = 0; i < TIMESIZE; i++)
    ftimes[i] = 0;

  /**********************************/
  /* Read geometry                  */
  /**********************************/

  starttimer;

  indsys = (SYS *)MattAlloc(1, sizeof(SYS));
  indsys->externals = NULL;
  indsys->num_extern = 0;
  indsys->segment = NULL;
  indsys->nodes = NULL;
  indsys->planes = NULL;                               /* CMS 7/2/92 */
  indsys->logofstep = 0;

  indsys->opts = Parse_Command_Line(argc, argv);

  /*  indsys->r_height = indsys->r_width = indsys->opts->filratio;*/

  /* set the type equal to that specified on command line for now */
  indsys->precond_type = indsys->opts->precond;

  /* degugging counters for calls to functions used in mutual() */
  num_exact_mutual = 0;
  num_fourfil = 0;
  num_mutualfil = 0;
  num_found = 0;
  num_perp = 0;
  
  opts = indsys->opts;

  tol = opts->tol;
  user_maxiters = opts->maxiters;

  if (opts->fname == NULL) {
    viewprintf(stdout, "No input file given.  Nothing to do.\n"); // Enrico
    //fp = stdin;
	FHOnClosing(FH_NO_MTX_NORMAL_END);  // Enrico 
	return FH_NO_MTX_NORMAL_END;       // Enrico
  }
  else if (strcmp(opts->fname,"-") == 0) {
    viewprintf(stdout, "Reading from stdin...\n");
    fp = stdin;
  }
  else {
    fp = fopen(opts->fname, "r");
    if (fp == NULL) {
      viewprintf(stderr, "Couldn't open %s\n", opts->fname);
      FHExit(FH_GENERIC_ERROR);
    }
    viewprintf(stdout, "Reading from file: %s\n",opts->fname);
  }

  FHSetName(opts->fname); // Enrico 

  // change the working directory to the one containing
  // the input file, to solve output directory bug when
  // called from FastModel. Enrico.
  _fullpath(full, opts->fname, _MAX_PATH);
  _splitpath(full, drive, dir, NULL, NULL);
  _makepath(full, drive, dir, NULL, NULL);
  _chdir(full);

  /* read in geometry */
  err = readGeom(fp, indsys);
  if (err != 0) {
	FHOnClosing(FH_GENERIC_ERROR);  // Enrico 
	return FH_GENERIC_ERROR;
  }

  fclose(fp);

  /*  fprintf(stdout,"Number of nodes before breaking up: %d\n",
      indsys->num_nodes); */
  
  /* regurgitate input file to stdout */
  if (opts->regurgitate == TRUE) {
    regurgitate(indsys);
  }

  if ((opts->makeFastCapFile & HIERARCHY) 
      && !(opts->makeFastCapFile & (SIMPLE | REFINED))) {
    /* the hierarchy has been dumped already, and nothing more to do. so quit*/
	  FHOnClosing(FH_NO_MTX_NORMAL_END);  // Enrico 
	  return FH_NO_MTX_NORMAL_END;
  }

  /* make a file suitable for keith's postscript maker */

  if (opts->makeFastCapFile & SIMPLE) {
    viewprintf(stdout, "Making simple zbuf file...");
    fflush(stdout);
    concat4(outfname,"zbuffile",opts->suffix,"");
    concat4(outfname2,outfname,"_","shadings");
    writefastcap(outfname, outfname2, indsys);
    viewprintf(stdout, "Done\n");
    if (! (opts->makeFastCapFile & REFINED) ) {
      /* we are done, after making fastcap files for visualization */
		FHOnClosing(FH_NO_MTX_NORMAL_END);  // Enrico 
		return FH_NO_MTX_NORMAL_END;
	}
  }

  /* initialize Gaussian quadrature arrays for mutual(). */
  /* See dist_betw_fils.c */
  fill_Gquad();

  /* initialize Freeable allocation stuff for mutual terms lookup table */
  init_table();

  /* break each segment into filaments */
  num_fils = 0;
  chgend = &chgdummy;

  for(seg = indsys->segment; seg != NULL; seg = seg->next) 
    chgend = assignFil(seg, &num_fils, chgend);

  indsys->num_fils = num_fils;

  chgend->next = NULL;
  chglist = chgdummy.next;

  stoptimer;
  ftimes[0] = dtime;

  /**********************************/
  /* Multipole setup                */
  /**********************************/

  starttimer;

  /* set up multipole stuff */
  sys = SetupMulti(chglist, indsys);

  if (bFHContinue == FALSE) {
	FHOnClosing(FH_USER_BREAK);  // Enrico 
	return FH_USER_BREAK;
  }

#if 1 == 0
  if (indsys->opts->debug == ON && 1 == 0)  /* for debugging eval pass */
    dump_evalcnts(sys);
#endif

  stoptimer;
  ftimes[5] = dtime;

  if (indsys->opts->debug == ON) 
    viewprintf(stdout, "Time for Multipole Setup: %lg\n",dtime);

  /**********************************/
  /* Scanning graph                 */
  /**********************************/

  starttimer;

  viewprintf(stdout, "Scanning graph to find fundamental circuits...\n");
  /* find all the necessary meshes */
  make_trees(indsys);

  /* clear node->examined and seg->is_deleted */
  clear_marks(indsys);

  /* find meshes in groundplane due to holes and put them at the front
     of the list of trees so that they are created by fillM() first so
     that they are marked before regular ground plane meshes  */
  find_hole_meshes(indsys);

  /* determine if a subset of columns is to be computed and assign 
     col_Yindex.  This will happen if -x option is used. */
  num_sub_extern = indsys->num_sub_extern
                       = pick_subset(opts->portlist, indsys);

  stoptimer;
  ftimes[6] = dtime;

  /*                                       */
  /* Shifted below code to write to Zc.mat */
  /* Enrico                                */
  
  /* count nodes */
  indsys->num_real_nodes = 0;
  indsys->num_nodes = 0;
  last = -1;
  for(node = indsys->nodes, i = 0; node != NULL; node = node->next, i++) {
    indsys->num_nodes++;
    if (getrealnode(node) == node) {
      indsys->num_real_nodes++;
      if (last + 1 != i) {
      /* printf("Non equivalent nodes must be listed first??, 
	 no I take that back.\n"); */
      /* exit(1); */
      }
      last = i;
    }
  }
  
  /* CMS 7/3/92 ------------------------------------------------------------*/
  planemeshes = 0;
  nonp = 0;

  for(plane = indsys->planes; plane != NULL; plane = plane->next){
    planemeshes = planemeshes + plane->numesh;
    if(plane->external == 0)
      nonp++;
  }
  /*------------------------------------------------------------------------*/

/* moved lower for shading 
  if (opts->makeFastCapFile & REFINED) {
    fprintf(stdout, "Making refined zbuf file...");
    fflush(stdout);

    concat4(outfname,"zbuffile2",opts->suffix,"");
    

    concat4(outfname2,outfname,"_","shadings");
    writefastcap(outfname, outfname2, indsys);
    fprintf(stdout, "Done\n");
  }
*/

  num_fils = indsys->num_fils;      /* setup multi may alter this */
  num_segs = indsys->num_segs;
  num_nodes = indsys->num_nodes;
  num_real_nodes = indsys->num_real_nodes;
  num_planes = indsys->num_planes;                    /* CMS 7/2/92 */

  /* figure out number of meshes */
  tree_meshes = indsys->tree_meshes = count_tree_meshes(indsys->trees);

  indsys->extra_meshes = tree_meshes;

#if 1==0 
  unimplemented junk
  indsys->extra_meshes = estimate_extra_meshes(indsys->trees, FILS_PER_MESH);
#endif

  num_mesh = num_fils - num_segs + indsys->extra_meshes + planemeshes;
                                                          /* CMS 7/2/92 */
  indsys->num_mesh = num_mesh;
  num_extern = count_externals(indsys->externals);
  if (num_extern != indsys->num_extern) {
    viewprintf(stderr, "main:  discrepancy in num_extern and actual number\n");
    FHExit(FH_GENERIC_ERROR);
  }
  fmin = indsys->fmin;
  fmax = indsys->fmax;
  logofstep = indsys->logofstep;

  /*  dont_form_Z = fmin == 0 && opts->mat_vect_prod == DIRECT 
                        && opts->soln_technique == ITERATIVE; */

  dont_form_Z = fmin == 0 && opts->soln_technique == LUDECOMP;
  indsys->dont_form_Z = dont_form_Z;


  viewprintf(stdout, "Number of Groundplanes : %d \n",num_planes);        /* CMS 7/2/92 */
  viewprintf(stdout, "Number of filaments: %10d\n", num_fils);
  viewprintf(stdout, "Number of segments:  %10d\n", num_segs);
  viewprintf(stdout, "Number of nodes:     %10d\n", num_nodes);
  viewprintf(stdout, "Number of meshes:    %10d\n", num_mesh);
  viewprintf(stdout, "          ----from tree:                   %d \n",tree_meshes);
  viewprintf(stdout, "          ----from planes: (before holes)  %d \n",planemeshes);  
                                                            /* CMS 7/6/92 */
  viewprintf(stdout, "Number of conductors:%10d   (rows of matrix in Zc.mat) \n",
	 num_extern);
  viewprintf(stdout, "Number of columns:   %10d   (columns of matrix in Zc.mat) \n",
	 num_sub_extern);

  viewprintf(stdout, "Number of real nodes:%10d\n", num_real_nodes);  /* non equivalent */
  
/*  A = MatrixAlloc(num_real_nodes, num_fils, sizeof(double)); */
  if (opts->mat_vect_prod == DIRECT || opts->soln_technique == LUDECOMP) {
/*  indsys->M = MatrixAlloc(num_fils, num_mesh, sizeof(double)); */
    if (!dont_form_Z)
      indsys->Z = MatrixAlloc(num_fils, num_fils, sizeof(double));

    /* let's save space and allocate this later
    /* indsys->MtZM = (CX **)MatrixAlloc(num_mesh, num_mesh, sizeof(CX)); */
  }
  indsys->FinalY = (CX **)MatrixAlloc(num_extern, num_sub_extern, sizeof(CX));
  b = (CX *)MattAlloc(num_mesh,sizeof(CX));
  x0 = (CX *)MattAlloc(num_mesh,sizeof(CX));
  indsys->R = (double *)MattAlloc(num_fils,sizeof(double));
  indsys->meshsect = (int *)MattAlloc((num_extern+nonp+1),sizeof(int));
  vect = (CX *)MattAlloc(num_mesh,sizeof(CX));
  pvect = (CX **)MattAlloc(num_mesh,sizeof(CX *));
  indsys->Mlist = (MELEMENT **)MattAlloc(num_mesh,sizeof(MELEMENT *));
  indsys->Mtrans = (MELEMENT **)MattAlloc(num_fils,sizeof(MELEMENT *));
  indsys->m_info = (Minfo *)MattAlloc(num_mesh,sizeof(Minfo));
  indsys->Precond = (PRE_ELEMENT **)MattAlloc(num_mesh, sizeof(PRE_ELEMENT *));
  if (1 == 1) {  /* let's always dump the residual history */
    indsys->resids = MatrixAlloc(num_sub_extern,user_maxiters, sizeof(double));
    indsys->resid_real = MatrixAlloc(num_sub_extern, user_maxiters, 
				     sizeof(double));
    indsys->resid_imag = MatrixAlloc(num_sub_extern, user_maxiters, 
				     sizeof(double));
    indsys->niters = (double *)MattAlloc(num_sub_extern, sizeof(double));
  }

/*
  if (opts->mat_vect_prod == DIRECT || opts->soln_technique == LUDECOMP) {
    for(i = 0; i < num_fils; i++) {
      for(j = 0; j < num_mesh; j++)
	indsys->M[i][j] = 0;
    }
  }
*/

  meshsect = indsys->meshsect;
  R = indsys->R;

  /**********************************/
  /* Form A M and Z                 */
  /**********************************/

  starttimer;
/*
  printf("filling A...\n");
  fillA(segment, A, num_segs);
*/

  viewprintf(stdout, "filling M...\n");
  fillM(indsys);  /* expects hole meshes to be marked */

  if (opts->makeFastCapFile & REFINED) {
    viewprintf(stdout, "Making refined zbuf file...");
    fflush(stdout);
    concat4(outfname,"zbuffile2",opts->suffix,"");

    /* shadings file */
    concat4(outfname2,outfname,"_","shadings");

    /* output file to later turn into .ps. Note: this affects seg->is_deleted*/
    writefastcap(outfname, outfname2, indsys);
    viewprintf(stdout, "Done\n");
	FHOnClosing(FH_NO_MTX_NORMAL_END);  // Enrico 
	return FH_NO_MTX_NORMAL_END;  /* we are done after making fastcap files for visualization */;
  }


  if (num_mesh != indsys->num_mesh)
    viewprintf(stdout, "Exact number of meshes after ground plane holes removed: %d\n",
	   indsys->num_mesh);

  num_mesh = indsys->num_mesh;  /* actual number of meshes */

  /* let's save a little space and allocate MZMt now. */
  if (!dont_form_Z 
      && (opts->mat_vect_prod == DIRECT || opts->soln_technique == LUDECOMP))
    indsys->MtZM = (CX **)MatrixAlloc(num_mesh, num_mesh, sizeof(CX));

  MtZM = indsys->MtZM;

  if (opts->dumpMats & MESHES) {
    if (opts->kind & MATLAB)
      dump_mesh_coords(indsys);
    else
      dump_ascii_mesh_coords(indsys);
  }

  formMtrans(indsys); /* Form M transpose by row*/

  viewprintf(stdout, "filling R and L...\n");
  fillZ(indsys);

  stoptimer;
  ftimes[1] = dtime;

  if (indsys->opts->debug == ON) 
    viewprintf(stdout, "Time to Form M and Z: %lg\n",dtime);

  viewprintf(stdout, "Total Memory allocated: %d kilobytes\n",memcount/1024);

  if (indsys->opts->debug == ON) 
    viewprintf(stdout, "Memory used and freed by lookup table: %d kilobytes\n",
	   get_table_mem());

  /* free memory for lookup table */
  destroy_table();

  /* This 'starttimer' macro seems not matched by any 'stoptime'; it does not */
  /* produce any error, since the next 'starttimer' will reset the timer, however */
  /* it probably is a bug, Enrico */
  starttimer;
  choose_and_setup_precond(indsys);

  if (opts->dumpMats) {
    viewprintf(stdout, "saving some files disk...\n");
    savemats(indsys);
  }

  if (logofstep == 0.0) {
    viewprintf(stderr, "no frequency range data read!\n");
    FHExit(FH_GENERIC_ERROR);
  }

  if (fmin == 0) {
      viewprintf(stderr, "***First frequency is zero. Only the DC case will be run.***\n");
      if (!dont_form_Z)
        viewprintf(stderr, "Warning: First frequency is zero, but -sludecomp was not specified.\n\
      Use this setting to save time and memory.\n");
  }

  if (indsys->opts->debug == ON) {
    /* open Ycond.mat */
    concat4(outfname,"Ycond",opts->suffix,".mat");
//    fp = fopen(outfname, "w");
      fp = fopen(outfname, "wb");
    if (fp == NULL) {
      viewprintf(stderr, "couldn't open file %s\n",outfname);
      FHExit(FH_GENERIC_ERROR);
    }
  }

  /* open b.mat */
  if (opts->dumpMats != OFF) {
    concat4(outfname,"b",opts->suffix,".mat");
//    if ((fb = fopen(outfname,"w")) == NULL)
      if ((fb = fopen(outfname,"wb")) == NULL)
        viewprintf(stderr, "No open fb\n");
  }

  /* Model Order Reduction: create, compute and print the model if requested */
  if (opts->orderROM > 0) {

    /* create the matrix whose inverse we will need */
    createMRMt(&MRMt, indsys);
    /* now create the input and output matrices for the state */
    B = MatrixAlloc(num_mesh, num_extern, sizeof(double));
    C = MatrixAlloc(num_mesh, num_extern, sizeof(double));
    for(ext = indsys->externals, j = 0;
        ext != NULL; ext = ext->next, j++) {
      int_list *elem;
      /* create C first; C = N, where Vs = N Vterminal */
      for(elem = ext->indices; elem != NULL; elem = elem->next) {
        if (elem->index > num_mesh || ext->Yindex > num_extern) {
          viewprintf(stderr, "Indexing into matrix C out of bounds\n");
	      FHExit(FH_GENERIC_ERROR);
		}
        C[elem->index][ext->Yindex] = elem->sign;
        B[elem->index][ext->Yindex] = elem->sign;
      }
    }

    if (indsys->opts->mat_vect_prod != MULTIPOLE && !indsys->dont_form_Z) {
      viewprintf(stdout, "multiplying M*(L)*transpose(M)for model order reduction\n");
      /* this form of storage isn't the best */
      formMLMt(indsys);      /*form (M^t)*(L)*M and store in indys->MtZM*/
      

      actual_order = ArnoldiROM(B, C, (double **)NULL, MRMt, num_mesh,
                                num_extern, num_extern, opts->orderROM, 
                                realMatVect, indsys, sys, chglist);
    }
    else if (indsys->opts->mat_vect_prod == MULTIPOLE)
      actual_order = ArnoldiROM(B, C, (double **)NULL, MRMt, num_mesh, 
                                num_extern, num_extern, opts->orderROM, 
                                realComputePsi, indsys, sys, chglist); 

    if (indsys->opts->debug == ON) {
#if 1==0
      /* open orig.mat */
      concat4(outfname,"orig",opts->suffix,".mat");
//      if ((fROM = fopen(outfname,"w")) == NULL)
	  if ((fROM = fopen(outfname,"wb")) == NULL)	// Enrico
        viewprintf(stderr, "No open fROM\n");
      /* dump what we have of the original system */
      dumpROMbin(fROM, NULL, B, C, NULL,
                 num_mesh, num_extern, num_extern);
      fclose(fROM);
#endif

    /* open rom.m */
      concat4(outfname,"rom",opts->suffix,".m");
      if ((fROM = fopen(outfname,"w")) == NULL)
        viewprintf(stderr, "No open fROM\n");

    /* now dump the reduced order model */
      dumpROM(fROM, indsys->Ar, indsys->Br, indsys->Cr, indsys->Dr,
              actual_order * num_extern, num_extern, num_extern);
      /* close and save away the file */
      fclose(fROM);
    }

    /* generate equivalent circuit */
    concat4(outfname,"equiv_circuitROM",opts->suffix,".spice");
    if ((fROM = fopen(outfname,"w")) == NULL)
      viewprintf(stderr, "No open fROM\n");
    /* now dump the reduced order model */
    dumpROMequiv_circuit(fROM, indsys->Ar, indsys->Br, indsys->Cr, indsys->Dr,
            actual_order * num_extern, num_extern, num_extern, 
            indsys->title, opts->suffix, indsys);
    /* close and save away the file */
    fclose(fROM);
    /* end of Model Order Reduction generation */
  }

  /* NOTE: exit here if all we want is the reduced order model (ROM) */
  if (opts->orderROM > 0 && opts->onlyROM)
    FHExit(FH_NORMAL_END);

  /* The next section of code (row / col names) has been shifted here from above, Enrico */

  // save the number of rows / cols
  strctImpMatrix.m_lRowNum = num_extern;
  strctImpMatrix.m_lColNum = num_sub_extern;
  // allocate list of row / columns names, Enrico
  strctImpMatrix.m_sRowNames = (char **)MattAlloc(strctImpMatrix.m_lRowNum, sizeof(char*));
  strctImpMatrix.m_sColNames = (char **)MattAlloc(strctImpMatrix.m_lColNum, sizeof(char*));

  /* Write to Zc.mat only if we aren't only running for visualization */
  if ( !(opts->makeFastCapFile & (SIMPLE | REFINED)) 
       && !(opts->orderROM > 0 && opts->onlyROM)   ) {

	concat4(outfname,"Zc",opts->suffix,".mat");   /* put filnames together */

    fp3 = fopen(outfname, "w");
    if (fp3 == NULL) {
     viewprintf(stderr, "couldn't open file %s\n",outfname);
      FHExit(FH_GENERIC_ERROR);
    }

    for(ext = indsys->externals; ext != NULL; ext=ext->next) {
      /* printf("Row %d :  %s  to  %s\n",ext->Yindex,ext->source->node[0]->name,
         ext->source->node[1]->name); */

	  // alloc memory for row name, Enrico
      strctImpMatrix.m_sRowNames[ext->Yindex] = (char *)MattAlloc(FH_MAX_NAME_LEN, sizeof(char));

	  if (ext->portname[0] == '\0') {
        fprintf(fp3,"Row %d:  %s  to  %s\n",ext->Yindex+1,ext->name1,ext->name2);
		// store row name, Enrico
		sprintf(strctImpMatrix.m_sRowNames[ext->Yindex], "%s to %s", ext->name1, ext->name2);
	  }
      else {
        fprintf(fp3,"Row %d:  %s  to  %s, port name: %s\n",
                ext->Yindex+1,ext->name1,ext->name2,ext->portname);
		// store row name, Enrico
		sprintf(strctImpMatrix.m_sRowNames[ext->Yindex], "%s to %s, %s", ext->name1, ext->name2, ext->portname);
	  }
    }

    if (num_sub_extern != indsys->num_extern)
      for(ext = indsys->externals; ext != NULL; ext=ext->next) {
		if (ext->col_Yindex != -1) {
          fprintf(fp3,"Col %d: port name: %s\n", 
                  ext->col_Yindex+1,ext->portname);
	      // alloc memory for col name, Enrico
          strctImpMatrix.m_sColNames[ext->col_Yindex] = (char *)MattAlloc(FH_MAX_NAME_LEN, sizeof(char));
		  // store col name, Enrico
		  sprintf(strctImpMatrix.m_sColNames[ext->col_Yindex], "%s", ext->portname);
		}
      }
  }

  // Count the number of frequency points, Enrico
  // uses 'for' cycle to mimic FastHenry behaviour, otherwise calculation
  // is not trivial
  for(freq = fmin, strctImpMatrix.m_lFreqNum = 0; 
      (fmin != 0 && freq <= fmax*1.001) || (fmin == 0 && strctImpMatrix.m_lFreqNum == 0);
      strctImpMatrix.m_lFreqNum++, freq = (fmin != 0 ? pow(10.0,log10(fmin) + strctImpMatrix.m_lFreqNum*logofstep) : 0.0)) ;

  // Allocate memory to store the impedance matrices, Enrico
  // Use FastHenry memory allocation routines so in case of premature
  // termination (i.e. user stop or error) can free memory (no leakages)
  strctImpMatrix.m_daMatricesReal = (double ***)MattAlloc(strctImpMatrix.m_lFreqNum, sizeof(double**));
  strctImpMatrix.m_daMatricesImag = (double ***)MattAlloc(strctImpMatrix.m_lFreqNum, sizeof(double**));
  // alloc memory for frequencies
  strctImpMatrix.m_daFrequencies = (double *)MattAlloc(strctImpMatrix.m_lFreqNum, sizeof(double));


  for(freq = fmin, m = 0; 
      (fmin != 0 && freq <= fmax*1.001) || (fmin == 0 && m == 0);
      m++, freq = (fmin != 0 ? pow(10.0,log10(fmin) + m*logofstep) : 0.0)) {
    viewprintf(stdout, "Frequency = %lg\n",freq);

    /**********************************/
    /* Form M'ZM                      */
    /**********************************/

    starttimer;

    if (!dont_form_Z && (opts->mat_vect_prod == DIRECT || opts->soln_technique==LUDECOMP)) {

	  viewprintf(stdout, "multiplying M*(R + jL)*transpose(M)\n");
      formMZMt(indsys);      /*form transpose(M)*(R+jL)*M  (no w) */

      if (opts->dumpMats & MZMt) {
	    if (m == 0) {
	      if (opts->kind & MATLAB) {
            concat4(outfname,"MZMt",opts->suffix,".mat");
//	        if ( (fp2 = fopen(outfname,"w")) == NULL)  
	        if ( (fp2 = fopen(outfname,"wb")) == NULL) {
              viewprintf(stderr, "Couldn't open file\n");
	          FHExit(FH_GENERIC_ERROR);
			}
	        viewprintf(stdout, "Saving MZMt...\n");
	        savecmplx2(fp2,"MZMt",indsys->MtZM, indsys->num_mesh,indsys->num_mesh);
	        fclose(fp2);
		  }
	      if (opts->kind & TEXT) {
            concat4(outfname,"MZMt",opts->suffix,".dat");
	        if ( (fp2 = fopen(outfname,"w")) == NULL) {
	          viewprintf(stderr, "Couldn't open file\n");
	          FHExit(FH_GENERIC_ERROR);
			}
	        cx_dumpMat_totextfile(fp2, indsys->MtZM,
				  indsys->num_mesh,indsys->num_mesh );
	        fclose(fp2);
		  }
		}
      }

      viewprintf(stdout, "putting in frequency \n");
      
      /* put in frequency */
      for(i = 0; i < num_mesh; i++)
	    for(j = 0; j < num_mesh; j++)
	      MtZM[i][j].imag *= 2*PI*freq;
	}

    stoptimer;
    ftimes[2] += dtime;

    /**********************************/
    /* Form precond                   */
    /**********************************/

    starttimer;

    if (opts->soln_technique == ITERATIVE) {
      if (indsys->precond_type == LOC) {
	    viewprintf(stdout, "Forming local inversion preconditioner\n");
	
	    if (opts->mat_vect_prod == DIRECT)
	      indPrecond_direct(sys, indsys, 2*PI*freq);
	    else if (opts->mat_vect_prod == MULTIPOLE)
	      indPrecond(sys, indsys, 2*PI*freq);
	    else {
	      viewprintf(stderr, "Internal error: mat_vect_prod == %d\n",
		  opts->mat_vect_prod);
	      FHExit(FH_GENERIC_ERROR);
		}
      }
      else if (indsys->precond_type == SPARSE) {
	    viewprintf(stdout, "Forming sparse matrix preconditioner..\n");
	    // TBC warning Enrico: this is a lenghty operation, but maybe doesn't 
	    // need bFHContinue check
	    fill_spPre(sys, indsys, 2*PI*freq);

        if (bFHContinue == FALSE) {
	      FHOnClosing(FH_USER_BREAK);  // Enrico 
		  return FH_USER_BREAK;
		}

	    if (m == 0) {
	      if (indsys->opts->debug == ON)
	        viewprintf(stdout, "Number of nonzeros before factoring: %d\n",
		               spElementCount(indsys->sparMatrix));

	      /* Reorder and Factor the matrix */
	      // TBC warning Enrico: this is a lenghty operation, but bFHContinue
	      // check is skipped because it is precompiled in library Sparse;
	      // should modify and recompile
	      err = spOrderAndFactor(indsys->sparMatrix, NULL, 1e-3, 0.0, 1);

          if (bFHContinue == FALSE) {
	        FHOnClosing(FH_USER_BREAK);  // Enrico 
			return FH_USER_BREAK;
		  }

          if (indsys->opts->debug == ON)
	        viewprintf(stdout, "Number of fill-ins after factoring: %d\n",
		               spFillinCount(indsys->sparMatrix));
		}
	    else
	      err = spFactor(indsys->sparMatrix);

	    if (err != 0) {
	      viewprintf(stderr,"Error on factor: %d\n",err);
	      FHExit(FH_GENERIC_ERROR);
		}
	  }
	}
    else if (opts->soln_technique == LUDECOMP) {
      if (!dont_form_Z) {
        viewprintf(stdout, "Performing LU Decomposition...\n");
        cx_ludecomp(MtZM, num_mesh, FALSE);
        viewprintf(stdout, "Done.\n");
      }
      else {
        /* let's form a sparse version since fmin=0 and the L matrix = 0 */
        create_sparMatrix(indsys);
        viewprintf(stdout, "Filling sparse version M R Mt...");
        fill_diagR(indsys);

        /* dump matrix to disk if requested */
        if (indsys->opts->dumpMats & MZMt) {
          viewprintf(stdout, "saving sparse MZMt to MZMt.mat...\n");
          concat4(outfname,"MZMt",indsys->opts->suffix,".mat");
          if (spFileMatrix(indsys->sparMatrix, outfname, "MZMt", 0, 1, 1) == 0)
            viewprintf(stderr,"saving sparse matrix failed\n");
		}

        if (indsys->opts->debug == ON)
          viewprintf(stdout, "Number of nonzeros before factoring: %d\n",
                     spElementCount(indsys->sparMatrix));

        /* Reorder and Factor the matrix */
        /* since w = 0, this matrix is real but all complex computations 
           are done since that is how the sparse package was compiled.
           But changing it doesn't help that much.  Must be dominated by
           reordering and such.*/
        viewprintf(stdout, "Factoring the sparse matrix...\n");
        err = spOrderAndFactor(indsys->sparMatrix, NULL, 1e-3, 0.0, 1);

        viewprintf(stdout, "Done factoring.\n");
        if (indsys->opts->debug == ON)
          viewprintf(stdout, "Number of fill-ins after factoring: %d\n",
                     spFillinCount(indsys->sparMatrix));
	  }
	}
    else {
      viewprintf(stderr, "Internal error: Unknown type of soln_technique: %d\n",
	             opts->soln_technique);
      FHExit(FH_GENERIC_ERROR);
    }

    stoptimer;
    ftimes[3] += dtime;

    if (indsys->opts->debug == ON)
      viewprintf(stdout, "Time spent on forming Precond: %lg\n",dtime);

    /**********************************/
    /* GMRES                          */
    /**********************************/

    starttimer;
    for(ext = get_next_ext(indsys->externals), i=0; ext != NULL; 
	                       ext = get_next_ext(ext->next),i++) {
      viewprintf(stdout, "conductor %d from node %s\n",i, get_a_name(ext->source));

      /* initialize b */
      for(j = 0; j < num_mesh; j++)
	    b[j] = x0[j] = CXZERO;
      
      fill_b(ext, b);

      sprintf(fname, "b%d_%d",m,i);
      if (opts->dumpMats != OFF)
	    savecmplx(fb, fname, &b, 1, num_mesh);

      maxiters = MIN(user_maxiters, num_mesh+2);

#if OPCNT == ON
      maxiters = 2;
#endif

      if (opts->soln_technique == ITERATIVE) {
	    viewprintf(stdout, "Calling gmres...\n");
	    if (opts->mat_vect_prod == MULTIPOLE) 
	      gmres(MtZM, b, x0, inner, SetupComputePsi, num_mesh, maxiters, tol,
		        sys, chglist, 2*PI*freq, R, indsys, i);
	    else 
	      gmres(MtZM, b, x0, inner, directmatvec,    num_mesh, maxiters, tol,
		        sys, chglist, 2*PI*freq, R, indsys, i);

        if (bFHContinue == FALSE) {
	      FHOnClosing(FH_USER_BREAK);  // Enrico 
		  return FH_USER_BREAK;
		}

	    if (indsys->precond_type == LOC) {
	      multPrecond(indsys->Precond, x0, vect, num_mesh);
	      for (k=0; k < num_mesh; k++)
	        x0[k] = vect[k];
		}
	    else if (indsys->precond_type == SPARSE)
	      spSolve(indsys->sparMatrix, x0, x0);
	  
	  }
      else 
        if (!dont_form_Z)
          cx_lu_solve(MtZM, x0, b, num_mesh);
        else 
          spSolve(indsys->sparMatrix, b, x0);

      if (opts->dumpMats & GRIDS)
	    /* Do stuff to look at groundplane current distribution */
	    makegrids(indsys, x0, ext->Yindex, m);

      /* Extract appropriate elements of x0 that correspond to final Yc */
      extractYcol(indsys->FinalY, x0, ext, indsys->externals);

      if (indsys->opts->debug == ON) {
//	    fptemp = fopen("Ytemp.mat", "w");
	    fptemp = fopen("Ytemp.mat", "wb");
	    if (fptemp == NULL) {
	      viewprintf(stderr, "couldn't open file %s\n","Ytemp.mat");
	      FHExit(FH_GENERIC_ERROR);
		}
	    strcpy(fname,"Ycond");
	    sprintf(tempstr, "%d", m);
	    strcat(fname,tempstr);
	    savecmplx(fptemp, fname, indsys->FinalY, num_extern, num_sub_extern);
	    fclose(fptemp);
      }
	}
    stoptimer;
    ftimes[4] += dtime;

    if (i != indsys->num_sub_extern) {
      viewprintf(stderr, "Huh?  columns calculated = %d  and num extern = %d\n",
	             i, indsys->num_sub_extern);
      FHExit(FH_GENERIC_ERROR);
    }

    if (indsys->opts->debug == ON) {
      /* dump matrix to matlab file */
      dump_to_Ycond(fp, m, indsys);
    }

    /* dump matrix to text file */
    if (num_extern == num_sub_extern) {
      fprintf(fp3, "Impedance matrix for frequency = %lg %d x %d\n ", freq,
	          num_extern, num_extern);
	  /* in-place invert the matrix using Gauss-Jordan (comment by Enrico) */
      cx_invert(indsys->FinalY, num_extern);
    }
    else /* if -x option used, dump admittance matrix instead of impedance (comment by Enrico) */
      fprintf(fp3, "ADMITTANCE matrix for frequency = %lg %d x %d\n ", freq,
	      num_extern, num_sub_extern);
    
	/* actually dump matrix to file */
    cx_dumpMat_totextfile(fp3, indsys->FinalY, num_extern, num_sub_extern);

	// store frequency at which matrix was computed, Enrico
	strctImpMatrix.m_daFrequencies[m] = freq;
	// alloc memory for storing matrix, Enrico
    strctImpMatrix.m_daMatricesReal[m] = (double **)MatrixAlloc(strctImpMatrix.m_lRowNum,
		                                  strctImpMatrix.m_lColNum, sizeof(double));
    strctImpMatrix.m_daMatricesImag[m] = (double **)MatrixAlloc(strctImpMatrix.m_lRowNum,
		                                  strctImpMatrix.m_lColNum, sizeof(double));
	// store matrix, Enrico; note that the imaginary part is divided by w*f (2*PI*freq)
	// to directly yeld inductance values
	for(iindex=0; iindex<num_extern; iindex++) {
	  for(jindex=0; jindex<num_sub_extern; jindex++) {
        strctImpMatrix.m_daMatricesReal[m][iindex][jindex] = indsys->FinalY[iindex][jindex].real;
		if(freq != 0)
          strctImpMatrix.m_daMatricesImag[m][iindex][jindex] = indsys->FinalY[iindex][jindex].imag / (2*PI*freq);
        else
          strctImpMatrix.m_daMatricesImag[m][iindex][jindex] = 0;
	  }
	}

    fflush(fp3);
  }
  
  if (indsys->opts->debug == ON)
    fclose(fp);

  fclose(fp3);

  if (opts->dumpMats != OFF)
    fclose(fb);

  //
  // Dump impedance matrices on the screen, Enrico
  //

  if(indsys->opts->output_mat == TRUE) {
    if(strctImpMatrix.m_lRowNum == strctImpMatrix.m_lColNum) {
	  viewprintf(stdout, "\nComputed matrices (R+jL)\n");
	}
    else {
  	  viewprintf(stdout, "\nComputed admittance matrices\n");
	}
    // print row/col names
    for(iindex=0; iindex<strctImpMatrix.m_lRowNum; iindex++) {
      viewprintf(stdout," Row %d:  %s\n", iindex, strctImpMatrix.m_sRowNames[iindex]);
	}
    if(strctImpMatrix.m_lRowNum != strctImpMatrix.m_lColNum) {
      for(iindex=0; iindex<strctImpMatrix.m_lColNum; iindex++) {
        viewprintf(stdout," Col %d:  %s\n", iindex, strctImpMatrix.m_sColNames[iindex]);
	  }
	}
    // print matrices
    for(iindex=0; iindex<strctImpMatrix.m_lFreqNum; iindex++) {
      viewprintf(stdout," Freq = %lg\n", strctImpMatrix.m_daFrequencies[iindex]);
      for(jindex=0; jindex<strctImpMatrix.m_lRowNum; jindex++) {
        viewprintf(stdout, "  Row %d: ", jindex);
        for(kindex=0; kindex<strctImpMatrix.m_lColNum; kindex++) {
	      viewprintf(stdout, "%lg%+lgj ", strctImpMatrix.m_daMatricesReal[iindex][jindex][kindex],
                                          strctImpMatrix.m_daMatricesImag[iindex][jindex][kindex]);
		}
        viewprintf(stdout, "\n");
	  }
	}
  }

  viewprintf(stdout, "\nAll impedance matrices dumped to file Zc%s.mat\n\n",opts->suffix);
  totaltime = 0;
  for(i = 0; i < TIMESIZE; i++) 
    totaltime += ftimes[i];

  if (indsys->opts->debug == ON) {
    viewprintf(stdout, "Calls to exact_mutual: %15d\n",num_exact_mutual);
    viewprintf(stdout, "         fourfils:     %15d\n",num_fourfil);
    viewprintf(stdout, "         mutualfil:    %15d\n",num_mutualfil);
    viewprintf(stdout, "Number found in table: %15d\n",num_found);
    viewprintf(stdout, "Number perpendicular:  %15d\n",num_perp);
    viewprintf(stdout, "\n");
  }

  viewprintf(stdout, "Times:  Read geometry   %lg\n",ftimes[0]);
  viewprintf(stdout, "        Multipole setup %lg\n",ftimes[5]);
  viewprintf(stdout, "        Scanning graph  %lg\n",ftimes[6]);
  viewprintf(stdout, "        Form A M and Z  %lg\n",ftimes[1]);
  viewprintf(stdout, "        form M'ZM       %lg\n",ftimes[2]);
  viewprintf(stdout, "        Form precond    %lg\n",ftimes[3]);
  viewprintf(stdout, "        GMRES time      %lg\n",ftimes[4]);
  viewprintf(stdout, "   Total:               %lg\n",totaltime);

#ifdef MATTDEBUG
  /* print memory bins */
  for(i = 0; i < 1001; i++)
    viewprintf(stderr, "%d\n", membins[i]);
#endif

  FHOnClosing(FH_NORMAL_END);  // Enrico 
  return FH_NORMAL_END;
}

/* this will divide a rectangular segment into many filaments */
charge *assignFil(seg, num_fils, chgptr)
SEGMENT *seg;
int *num_fils;
charge *chgptr;
{

  int i,j,k;
  double delw, delh;
  double hx, hy, hz;
  int temp, counter;
  int Hinc, Winc;
  double Hdiv, Wdiv;
  double height, width;
  double h_from_edge, w_from_edge, min_height, min_width;
  FILAMENT *tempfil, *filptr;
  int countfils;
  double wx, wy, wz, mag;  /* direction of width */
  NODES *node0, *node1;
  charge *clast, *chg;
  charge *ctemp;  /* temporary place so can allocate all the charges at once*/
  surface *dummysurf;
  int indices[4], row, col;
  double r_width, r_height;
            /* ratio of element sizes (for geometrically increasing size)*/
 
  Hinc = seg->hinc;
  Winc = seg->winc;
  r_height = seg->r_height;
  r_width = seg->r_width;

  clast = chgptr;

  seg->num_fils = Winc*Hinc;
  seg->filaments = (FILAMENT *)MattAlloc(Hinc*Winc,sizeof(FILAMENT));

  ctemp = (charge *)MattAlloc(Hinc*Winc, sizeof(charge));

  /* clear pchg as a marker for checking later */
  for(i = 0; i < Hinc*Winc; i++)
    seg->filaments[i].pchg = NULL;

  dummysurf = (surface *)MattAlloc(1, sizeof(surface));
  dummysurf->type = CONDTR;
  dummysurf->next = NULL;
  dummysurf->prev = NULL;
  
  /*  To make the filaments have geometrically decreasing areas */
  /* Hdiv and Wdiv are 1/(smallest Width) */

  if (fabs(1.0 - r_height) < EPS)
    Hdiv = Hinc;
  else {
    temp = Hinc/2;
    Hdiv = 2*(1.0 - pow(r_height, (double)temp))/(1.0 - r_height);
    if (Hinc%2 == 1) Hdiv += pow(r_height, (double)(temp));
  }

  if (fabs(1.0 - r_width) < EPS)
    Wdiv = Winc;
  else {
    temp = Winc/2;
    Wdiv = 2*(1.0 - pow(r_width, (double)temp))/(1.0 - r_width);
    if (Winc%2 == 1) Wdiv += pow(r_width, (double)(temp));
  }

  node0 = seg->node[0];
  node1 = seg->node[1];
  /* determine direction of width */
  if (seg->widthdir != NULL) {
    wx = seg->widthdir[XX];
    wy = seg->widthdir[YY];
    wz = seg->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(node1->y - node0->y)*1.0;
    wy = (node1->x - node0->x)*1.0;
    wz = 0;
    if ( fabs(wx/seg->length) < EPS && fabs(wy/seg->length) < EPS) {
      /* if all of x-y is perpendic to length, then choose x direction */
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }
  /* height direction perpendicular to both */
  hx = -wy*(node1->z - node0->z) + (node1->y - node0->y)*wz;
  hy = -wz*(node1->x - node0->x) + (node1->z - node0->z)*wx;
  hz = -wx*(node1->y - node0->y) + (node1->x - node0->x)*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;

  filptr = seg->filaments;
  counter = 0;

  /* this will fill the 'filament' array. It uses symmetry wrt z and y */
  /* it generates the four corner fils first, then the next one in...  */
  /* 6/25/93 - added stuff to place fils in filament array so that adjacent */
  /*           fils are adjacent in the array.  This will make the meshes   */
  /*           in M consist of fils that are near each other.               */
  h_from_edge = 0.0;
  min_height = 1.0/Hdiv;  /* fil of smallest height */
  for(i = 0; i < Hinc/2 || (Hinc%2 == 1 && i == Hinc/2); i++) {

    /* height of the filaments for this row */
    height = (seg->height/Hdiv)*pow(r_height, (double)i);
    /* delh = (seg->height)*(0.5 - 1.0/Hdiv*
			  ( (1 - pow(r_height, (double)(i+1)))/(1-r_height)
			   - 0.5*pow(r_height, (double) i) )
			 );
			 */
    if (i == 0)
      h_from_edge += min_height/2;
    else
      h_from_edge += min_height/2*pow(r_height,(double)(i-1))*(1+r_height);

    delh = (seg->height)*(0.5 - h_from_edge);

    if (delh < 0.0 && fabs(delh/seg->height) > EPS) 
     viewprintf(stderr, "uh oh, delh < 0. delh/height = %lg\n", delh/seg->height);
    
    w_from_edge = 0;
    min_width = 1.0/Wdiv;
    for(j = 0; j < Winc/2 || (Winc%2 == 1 && j == Winc/2); j++) {
      width = (seg->width/Wdiv)*pow(r_width, (double)j );
      /*delw = (seg->width)*
	(0.5 - 1.0/Wdiv*
	 ( (1 - pow(r_width, (double)(j+1)))/(1 - r_width)
	  - 0.5*pow(r_width, (double) j) )
	 );
	 */

      if (j == 0)
	w_from_edge += min_width/2;
      else
	w_from_edge += min_width/2*pow(r_width,(double)(j-1))*(1+r_width);
      
      delw = (seg->width)*(0.5 - w_from_edge);

      if (delw < 0.0 && fabs(delw/seg->width) > EPS) 
	printf("uh oh, delw < 0. delw/width = %lg\n", delw/seg->width);
/*    tempfil = filptr; */
      countfils = 0;
      
      row = i; 
      col = j;
      if (row%2 == 1)
	col = (Winc - 1) - col;
      indices[countfils] = col + Winc*row;
      filptr = &(seg->filaments[indices[countfils]]);

      filptr->x[0] = node0->x + hx*delh + wx*delw;
      filptr->x[1] = node1->x + hx*delh + wx*delw;
      filptr->y[0] = node0->y + hy*delh + wy*delw;
      filptr->y[1] = node1->y + hy*delh + wy*delw;
      filptr->z[0] = node0->z + hz*delh + wz*delw;
      filptr->z[1] = node1->z + hz*delh + wz*delw;
      filptr->pchg = &ctemp[counter++];
/*    filptr++; */
      countfils++;
      
      /* do symmetric element wrt y */
      if(j != Winc/2) {
	row = i; 
	col = (Winc - 1) - j;
	if (row%2 == 1)
	  col = (Winc - 1) - col;
	indices[countfils] = col + Winc*row;
	filptr = &(seg->filaments[indices[countfils]]);

	filptr->x[0] = node0->x + hx*delh - wx*delw;
	filptr->x[1] = node1->x + hx*delh - wx*delw;
	filptr->y[0] = node0->y + hy*delh - wy*delw;
	filptr->y[1] = node1->y + hy*delh - wy*delw;
	filptr->z[0] = node0->z + hz*delh - wz*delw;
	filptr->z[1] = node1->z + hz*delh - wz*delw;
	filptr->pchg = &ctemp[counter++];
/*	filptr++; */
	countfils++;
      }
      
      /* wrt z */
      if(i != Hinc/2) {
        row = (Hinc - 1) - i; 
	col = j;
	if (row%2 == 1)
	  col = (Winc - 1) - col;
	indices[countfils] = col + Winc*row;
	filptr = &(seg->filaments[indices[countfils]]);

	filptr->x[0] = node0->x - hx*delh + wx*delw;
	filptr->x[1] = node1->x - hx*delh + wx*delw;
	filptr->y[0] = node0->y - hy*delh + wy*delw;
	filptr->y[1] = node1->y - hy*delh + wy*delw;
	filptr->z[0] = node0->z - hz*delh + wz*delw;
	filptr->z[1] = node1->z - hz*delh + wz*delw;
	filptr->pchg = &ctemp[counter++];
	filptr++;
	countfils++;
      }
      
      /* wrt z and y */
      if( i != Hinc/2 && j != Winc/2) {
	row = (Hinc - 1) - i; 
	col = (Winc - 1) - j;
	if (row%2 == 1)
	  col = (Winc - 1) - col;
	indices[countfils] = col + Winc*row;
	filptr = &(seg->filaments[indices[countfils]]);

	filptr->x[0] = node0->x - hx*delh - wx*delw;
	filptr->x[1] = node1->x - hx*delh - wx*delw;
	filptr->y[0] = node0->y - hy*delh - wy*delw;
	filptr->y[1] = node1->y - hy*delh - wy*delw;
	filptr->z[0] = node0->z - hz*delh - wz*delw;
	filptr->z[1] = node1->z - hz*delh - wz*delw;
	filptr->pchg = &ctemp[counter++];
/*	filptr++; */
	countfils++;
      }
      
      for(k = 0; k < countfils; k++) {
	tempfil = &(seg->filaments[indices[k]]);
	tempfil->length = seg->length;
	tempfil->area = width*height;
	tempfil->width = width;
	tempfil->height = height;
	tempfil->filnumber = (*num_fils)++;
	tempfil->segm = seg;

	tempfil->lenvect[XX] = tempfil->x[1] - tempfil->x[0];
	tempfil->lenvect[YY] = tempfil->y[1] - tempfil->y[0];
	tempfil->lenvect[ZZ] = tempfil->z[1] - tempfil->z[0];
	/* do stuff for multipole */
	/*   make linked list entry */
	chg = tempfil->pchg;
	clast->next = chg;
	clast = chg;
	/* fill charge structure */
	chg->max_diag = chg->min_diag = tempfil->length;
	chg->x = (tempfil->x[0] + tempfil->x[1])/2.0;
	chg->y = (tempfil->y[0] + tempfil->y[1])/2.0;
	chg->z = (tempfil->z[0] + tempfil->z[1])/2.0;
	chg->surf = dummysurf;
	chg->dummy = FALSE;
	chg->fil = tempfil;

/*	tempfil++; */
      }
      
    }
  }

  i = 0;
  while(i < Hinc*Winc && seg->filaments[i].pchg != NULL)
    i++;

  if (i != Hinc*Winc) {
    viewprintf(stderr, "Hey, not all filaments created in assignfil()! \n");
    FHExit(FH_GENERIC_ERROR);
  }

  return clast;
}


double **MatrixAlloc(rows, cols, size)
int rows, cols, size;
{

  double **temp;
  int i;

  temp = (double **)MattAlloc(rows,sizeof(double *));
  if (temp == NULL) {
       viewprintf(stderr, "not enough space for matrix allocation\n");
        FHExit(FH_GENERIC_ERROR);
      }

  for(i = 0; i < rows; i++) 
    temp[i] = (double *)MattAlloc(cols,size);

  if (temp[rows - 1] == NULL) {
       viewprintf(stderr, "not enough space for matrix allocation\n");
        FHExit(FH_GENERIC_ERROR);
      }
  return(temp);
}

void fillA(indsys)
SYS *indsys;
{
  SEGMENT *seg;
  NODES *node1, *node2;
  MELEMENT **Alist;
  int i, counter;
  FILAMENT *fil;

  indsys->Alist = (MELEMENT **)MattAlloc(indsys->num_real_nodes,
					 sizeof(MELEMENT *));

  pick_ground_nodes(indsys);

  Alist = indsys->Alist;
  counter = 1;  /* ground is chosen already */

  for(seg = indsys->segment; seg != NULL; seg = seg->next) {
    node1 = getrealnode(seg->node[0]);
    node2 = getrealnode(seg->node[1]);
    if (node1->index == -1) {
      node1->index = counter++;
      Alist[node1->index] = NULL;
    }
    if (node2->index == -1) {
      node2->index = counter++;
      Alist[node2->index] = NULL;
    }
    for(i = 0; i < seg->num_fils; i++) {
      fil = &seg->filaments[i];
      Alist[node1->index] = insert_in_list(make_melement(fil->filnumber, 
							 fil, 1),
					   Alist[node1->index]);
      Alist[node2->index] = insert_in_list(make_melement(fil->filnumber, 
							 fil, -1),
					   Alist[node2->index]);	   
    }
  }

  if (counter != indsys->num_real_nodes - indsys->num_trees + 1) {
    viewprintf(stderr,"Internal error when forming A: counter %d != num_real_nodes %d\n",
	    counter, indsys->num_real_nodes);
  }

}

/* this fills the kircoff's voltage law matrix (Mesh matrix) */
/* it maps a matrix of mesh currents to branch currents */
/* it might actually be what some think of as the transpose of M */
/* Here, M*Im = Ib  where Im are the mesh currents, and Ib the branch */
/* 6/92 I added Mlist which is a vector of linked lists to meshes. 
   This replaces M.  But I keep M around for checking things in matlab. */
void old_fillM(indsys)
SYS *indsys;
{
}
    
void fillZ(indsys)
SYS *indsys;
{
  int j, m;
  FILAMENT *fil_j, *fil_m;
  int filnum_j, filnum_m;
  SEGMENT *seg1, *seg2;
  double **Z, *R;

  Z = indsys->Z;
  R = indsys->R;

  for(seg1 = indsys->segment; seg1 != NULL; seg1 = seg1->next) {
    for(j = 0; j < seg1->num_fils; j++) {
      fil_j = &(seg1->filaments[j]);
      filnum_j = fil_j->filnumber;
      R[filnum_j] = resistance(fil_j, seg1->sigma);
      if (indsys->opts->mat_vect_prod != MULTIPOLE 
          && !indsys->dont_form_Z) {
	for(seg2 = indsys->segment; seg2 != NULL; seg2 = seg2->next) {
	  for(m = 0; m < seg2->num_fils; m++) {
	    fil_m = &(seg2->filaments[m]);
	    filnum_m = fil_m->filnumber;
	    if (filnum_m == filnum_j)
	      Z[filnum_m][filnum_m] = selfterm(fil_m); /* do self-inductance */
	    else
	      if (filnum_m > filnum_j) /*we haven't done it yet */
		
		Z[filnum_m][filnum_j] 
		  = Z[filnum_j][filnum_m] = mutual(fil_j, fil_m);
	  }
	}
      }  /* end if (multipole or dont_form_Z) */
    }
  }
}

/* calculates resistance of filament */
double resistance(fil, sigma)
FILAMENT *fil;
double sigma;  /* conductivitiy */
{
  return  fil->length/(sigma * fil->area);
}

/* mutual inductance functions moved to mutual.c */

int matherr(exc)
struct exception *exc;
{
 viewprintf(stderr, "Err in math\n");
  return(0);
}

/* This counts the nonblank lines of the file  fp (unused) */
int countlines(fp)
FILE *fp;
{
  int count;
  char temp[MAXCHARS];

  count = 0;
  while( fgets(temp,MAXCHARS, fp) != NULL)
    if ( notblankline(temp) ) count++;

  return count;
}

/* returns 1 if string contains a nonspace character */
static int notblankline(string)
char *string;
{
   while( *string!='\0' && isspace(*string))
     string++;

   if (*string == '\0') return 0;
     else return 1;
}

/* This saves various matrices to files and optionally calls fillA() if
   the incidence matrix, A, is requested */
void savemats(indsys)
SYS *indsys;
{
  int i, j;
  FILE *fp, *fp2;
  FILAMENT *tmpf;
  SEGMENT *tmps;
  double *dumb, *dumbx, *dumby, *dumbz;
  int num_mesh, num_fils, num_real_nodes;
  double *Mrow;   /* one row of M */
  double **Z, *R;
  int machine, nnz, nnz0;
  ind_opts *opts;
  MELEMENT **Mlist = indsys->Mlist;
  MELEMENT **Alist;

  num_mesh = indsys->num_mesh;
  num_fils = indsys->num_fils;
  num_real_nodes = indsys->num_real_nodes;

  Mrow = (double *)MattAlloc(num_fils,sizeof(double));
  R = indsys->R;
  Z = indsys->Z;
  opts = indsys->opts;

  if (opts->dumpMats & DUMP_M) {
    /* count non-zero entries */
    nnz = nnz_inM(indsys->Mlist, num_mesh);

    if (opts->kind & TEXT) {
      concat4(outfname,"M",opts->suffix,".dat");
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	FHExit(FH_GENERIC_ERROR);
      }
      dump_M_to_text(fp2, Mlist, num_mesh, nnz);
      fclose(fp2);
    }

    if (opts->kind & MATLAB) {
      concat4(outfname,"M",opts->suffix,".mat");
//      if ( (fp = fopen(outfname,"w")) == NULL) {
      if ( (fp = fopen(outfname,"wb")) == NULL) {
	printf("Couldn't open file\n");
	FHExit(FH_GENERIC_ERROR);
      }
      dump_M_to_matlab(fp, Mlist, num_mesh, nnz, "M");
      fclose(fp);
    }
  }

  if (opts->dumpMats & DUMP_A) {

    /* fill the A matrix */
    fillA(indsys);
    Alist = indsys->Alist;

    /* count non-zero entries */
    nnz = nnz_inM(indsys->Alist, num_real_nodes);

    /* how many non-zeros in ground node 0? */
    nnz0 = nnz_inM(indsys->Alist, 1);  

    /* now dump it to a file without ground node */

    if (opts->kind & TEXT) {
      concat4(outfname,"A",opts->suffix,".dat");
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	FHExit(FH_GENERIC_ERROR);
      }
      dump_M_to_text(fp2, &Alist[1], num_real_nodes - 1, nnz - nnz0);
      fclose(fp2);
    }

    if (opts->kind & MATLAB) {
      concat4(outfname,"A",opts->suffix,".mat");
//      if ( (fp = fopen(outfname,"w")) == NULL) {
      if ( (fp = fopen(outfname,"wb")) == NULL) {
	printf("Couldn't open file\n");
	FHExit(FH_GENERIC_ERROR);
      }
      dump_M_to_matlab(fp, &Alist[1], num_real_nodes - 1, nnz - nnz0,"A");
      fclose(fp);
    }
  }

  /* save text files */
  if (opts->kind & TEXT && opts->dumpMats & DUMP_RL) {
    /*
    if ( (fp2 = fopen("M.dat","w")) == NULL) {
      printf("Couldn't open file\n");
      exit(1);
    }

    for(i = 0; i < num_fils; i++)
      Mrow[i] = 0;

    for(i = 0; i < num_mesh; i++) {
      fillMrow(indsys->Mlist, i, Mrow);
      dumpVec_totextfile(fp2, Mrow, num_fils);
    }
    
    fclose(fp2);
    */

    concat4(outfname,"R",opts->suffix,".dat");
    if ( (fp2 = fopen(outfname,"w")) == NULL) {
     viewprintf(stderr, "Couldn't open file\n");
      FHExit(FH_GENERIC_ERROR);
    }
    
    dumpVec_totextfile(fp2, R, num_fils);
    
    fclose(fp2);
    
    if (!indsys->dont_form_Z && indsys->opts->mat_vect_prod == DIRECT) {
      /* do imaginary part of Z */
      
      concat4(outfname,"L",opts->suffix,".dat");
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	FHExit(FH_GENERIC_ERROR);
      }
      
      dumpMat_totextfile(fp2, Z, num_fils, num_fils);
      
      fclose(fp2);
      
    }
    else {
      if (indsys->dont_form_Z)
       viewprintf(stderr, "L matrix not dumped since L = 0 since fmin = 0\n");
      else
       viewprintf(stderr, "L matrix not dumped.  Run with mat_vect_prod = DIRECT if this is desired\n");
    }
  }

  /* Save Matlab files */
  
  if (opts->kind & MATLAB && opts->dumpMats & DUMP_RL) {
    concat4(outfname,"RL",opts->suffix,".mat");
//    if ( (fp = fopen(outfname,"w")) == NULL) {
    if ( (fp = fopen(outfname,"wb")) == NULL) {
     viewprintf(stderr, "Couldn't open file\n");
      FHExit(FH_GENERIC_ERROR);
    }

    machine = 1000;
#ifdef DEC
    machine = 0000;
#endif
#ifdef WIN32		// Enrico
    machine = 0000;
#endif

    dumbx =  (double *)DMALCORE(num_fils*sizeof(double));	// Enrico, MALCORE instead of malloc
    dumby =  (double *)DMALCORE(num_fils*sizeof(double));	
    dumbz =  (double *)DMALCORE(num_fils*sizeof(double));	
    dumb =  (double *)DMALCORE(num_fils*sizeof(double));	
    if (dumb == NULL) {
      viewprintf(stderr, "no space for R\n");

	  DFREECORE(dumb);     // added to free memory (and close file) if exits beforehand
      DFREECORE(dumbx);
      DFREECORE(dumby);
      DFREECORE(dumbz);
      fclose(fp);

      FHExit(FH_GENERIC_ERROR);
    }
    
    /* save and print matrices */

    /*
    for(i = 0; i < num_fils; i++)
      Mrow[i] = 0;

    for(i = 0; i < num_mesh; i++) {
      fillMrow(indsys->Mlist, i, Mrow);
      savemat_mod(fp, machine+100, "M", num_mesh, num_fils, 0, Mrow, 
		  (double *)NULL, i, num_fils);
    }
    */

    /* this only saves the real part (savemat_mod.c modified) */
    savemat_mod(fp, machine+100, "R", 1, num_fils, 0, R, (double *)NULL, 
		0, num_fils);
    
    if (!indsys->dont_form_Z && indsys->opts->mat_vect_prod == DIRECT) {
      /* do imaginary part of Z */
      for(i = 0; i < num_fils; i++) {
	savemat_mod(fp, machine+100, "L", num_fils, num_fils, 0, Z[i], 
		    (double *)NULL, i, num_fils);
      }
    }
    else {
      if (indsys->dont_form_Z)
       viewprintf(stderr, "L matrix not dumped since L = 0 since fmin = 0\n");
      else
       viewprintf(stderr, "L matrix not dumped.  Run with mat_vect_prod = DIRECT if this is desired\n");
    }
    
    /* save vector of filament areas */
    for(tmps = indsys->segment; tmps != NULL; tmps = tmps->next) {
      for(j = 0; j < tmps->num_fils; j++) {
	tmpf = &(tmps->filaments[j]);
	dumb[tmpf->filnumber] = tmpf->area;
	dumbx[tmpf->filnumber] = tmpf->x[0];
	dumby[tmpf->filnumber] = tmpf->y[0];
	dumbz[tmpf->filnumber] = tmpf->z[0];
      }
    }
    savemat_mod(fp, machine+0, "areas", num_fils, 1, 0, dumb, (double *)NULL,0,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumbx, (double *)NULL, 0,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumby, (double *)NULL, 1,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumbz, (double *)NULL, 1,
		num_fils);
    
    DFREECORE(dumb);
    DFREECORE(dumbx);
    DFREECORE(dumby);
    DFREECORE(dumbz);
    
    fclose(fp);
  }
}

void savecmplx(fp, name, Z, rows, cols)
FILE *fp;
char *name;
CX **Z;
int rows, cols;
{
  int i,j;
  int machine;

  machine = 1000;
#ifdef DEC
  machine = 0000;
#endif
#ifdef WIN32		// Enrico
  machine = 0000;
#endif

  /* this only saves the real part (savemat_mod.c modified) one byte per call*/
  for(i = 0; i < rows; i++) 
    for(j = 0; j < cols; j++) 
      savemat_mod(fp, machine+100, name, rows, cols, 1, &Z[i][j].real, 
		  (double *)NULL, i+j, 1);

  /* do imaginary part of Z */
  for(i = 0; i < rows; i++) 
    for(j = 0; j < cols; j++)
      savemat_mod(fp, machine+100, name, rows, cols, 0, &Z[i][j].imag, 
		  (double *)NULL, 1, 1);
}

/* saves a complex matrix more efficiently? */
void savecmplx2(fp, name, Z, rows, cols)
FILE *fp;
char *name;
CX **Z;
int rows, cols;
{
  int i,j;
  int machine;
  static double *temp;
  // static int colmax = 0; // Enrico, moved at file scope

  if (colmax < cols) {
    temp = (double *)MALCORE(cols*sizeof(double));	// Enrico, MALCORE instead of malloc
    colmax = cols;
  }

  machine = 1000;
#ifdef DEC
  machine = 0000;
#endif
#ifdef WIN32		// Enrico
  machine = 0000;
#endif

  /* this only saves the real part (savemat_mod.c modified) one byte per call*/
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      temp[j] = Z[i][j].real;
    savemat_mod(fp, machine+100, name, rows, cols, 1, temp, 
		(double *)NULL, i, cols);
  }

  /* do imaginary part of Z */
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      temp[j] = Z[i][j].imag;
    savemat_mod(fp, machine+100, name, rows, cols, 0, temp, 
		(double *)NULL, 1, cols);
  }
}

/* This computes the product M*Z*Mt in a better way than oldformMZMt */
void formMZMt(indsys)
SYS *indsys;
{
  double tempsum, tempR;
  // Enrico, moved to file scope
  // static double *tcol = NULL;   /* temporary storage for extra rows */
  int i, j, mesh, mesh2;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  MELEMENT *mt, *mt2;    /* elements of M transpose */
  double **M, **L, *R;
  CX **Zm;
  int rows, cols, num_mesh;
  MELEMENT **Mlist, **Mt;

  Zm = indsys->MtZM;
  L = indsys->Z;
  R = indsys->R;
  M = indsys->M;
  nfils = rows = indsys->num_fils;
  nmesh = cols = num_mesh = indsys->num_mesh;
  Mlist = indsys->Mlist;
  Mt = indsys->Mtrans;

  if (nmesh > nfils) {
    viewprintf(stderr, "Uh oh, more meshes than filaments, I'm confused\n");
    FHExit(FH_GENERIC_ERROR);
  }

  /* allocate a temporary column for Z*Mt */
  if (tcol == NULL)
    tcol = (double *)MattAlloc(rows, sizeof(double));

  /* this does L*(Mt)_i where (Mt)_i is a single column (row) of Mt (M) and
     saves it in the temp space, tcol */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j< nfils; j++)
      tcol[j] = 0;
    /* note, this next set of nested loops could be reversed, but I think
       this might be more efficient for pipelining, ? */
    for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
      for(j = 0; j < nfils; j++)
	tcol[j] += L[j][melem->filindex]*melem->sign;
    for(mesh2 = 0; mesh2 < num_mesh; mesh2++) {
      tempsum = 0;
      for(melem2 = Mlist[mesh2]; melem2 != NULL; melem2 = melem2->mnext)
	tempsum += melem2->sign*tcol[melem2->filindex];
      Zm[mesh2][mesh].imag = tempsum;
    }
  }

  /* R is diagonal, so M*R*Mt is easy */
  for(i = 0; i < num_mesh; i++)
    for(j = 0; j < num_mesh; j++)
      Zm[i][j].real = 0;

  for(i = 0; i < nfils; i++)
    for(mt = Mt[i]; mt != NULL; mt = mt->mnext) {
      tempR = R[i]*mt->sign;
      for(mt2 = Mt[i]; mt2 != NULL; mt2 = mt2->mnext)
	Zm[mt2->filindex][mt->filindex].real += mt2->sign*tempR;
    }
}

oldformMZMt(indsys)
SYS *indsys;
{
  double tempsum;
  // Enrico, moved to file scope
  //static CX **trows = NULL;   /* temporary storage for extra rows */
  int i, j, mesh;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  double **M, **L, *R;
  CX **Zm;
  int rows, cols, num_mesh;
  MELEMENT **Mlist;

  Zm = indsys->MtZM;
  L = indsys->Z;
  R = indsys->R;
  M = indsys->M;
  nfils = rows = indsys->num_fils;
  nmesh = cols = num_mesh = indsys->num_mesh;
  Mlist = indsys->Mlist;

  if (nmesh > nfils) {
    viewprintf(stderr, "Uh oh, more meshes than filaments, I'm confused\n");
    FHExit(FH_GENERIC_ERROR);
  }

  /* allocate some extra rows so we can use Zm as temp space */
  if (rows > cols && trows == NULL)
    trows = (CX **)MatrixAlloc(rows - cols, cols, sizeof(CX));

  /* this does L*Mt and saves it in Zm.real.  This could be done better i
     think (use the fact that MZMt is symmetric)
     Also, could only track through each mesh list once, and store
     temporary values in Zm.real as we go along. (someday) */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j < nfils; j++) {
      tempsum = 0;
      for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext) 
	tempsum += L[j][melem->filindex]*melem->sign;
      if (j < nmesh)
	Zm[j][mesh].real = tempsum;
      else
	trows[j - nmesh][mesh].real = tempsum;
    }
  }

/*      savecmplx(fp, "step1", Zm, num_mesh, num_mesh); */

  /* this does M*(Zm.real) where Zm.real is hopefully L*Mt */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j < nmesh; j++) {
      tempsum = 0;
      for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext) {
	if (melem->filindex < nmesh)
	  tempsum += Zm[melem->filindex][j].real*melem->sign;
	else
	  tempsum += trows[melem->filindex - nmesh][j].real*melem->sign;
      }
      Zm[mesh][j].imag = tempsum;
    }
  }
/*      savecmplx(fp, "step2", Zm, num_mesh, num_mesh); */
  
  /* R is diagonal, so M*R*Mt is easy */
  for(i = 0; i < num_mesh; i++) {
    for(j = i; j < num_mesh; j++) {
      tempsum = 0;
      /* brute force search for elements */
      for(melem = Mlist[i]; melem != NULL; melem = melem->mnext) {
	for(melem2 = Mlist[j]; melem2 != NULL; melem2 = melem2->mnext) {
	  if (melem2->filindex == melem->filindex) 
	    tempsum += melem->sign*R[melem->filindex]*melem2->sign;
	}
      }
      Zm[i][j].real = Zm[j][i].real = tempsum;
    }
  }
  /* don't free the temp space */    
  /*
  if (rows > cols) {
    for(i = 0; i < rows - cols; i++)
      free(trows[i]);
    free(trows);
  }
  */
}

char* MattAlloc(number, size)
int number, size;
{

  char *blah;

/*  blah = (char *)malloc(number*size); */
  CALLOC(blah, number*size, char, OFF, IND);

  if (blah == NULL) {
    viewprintf(stderr, "MattAlloc: Couldn't get space. Needed %d\n",number*size);
    FHExit(FH_GENERIC_ERROR);
  }

  return blah;
}

/* forms the transpose of M.  Makes a linked list of each row.  Mlist is 
    a linked list of its rows. */
/* Note: This uses the same struct as Mlist but in reality, each linked list
   is composed of mesh indices, not fil indices. (filindex is a mesh index) */
void formMtrans(indsys)
SYS *indsys;
{
  int i, j;
  MELEMENT *m, *mt;
  MELEMENT **Mlist, **Mtrans;
  int meshes, fils;
  int last_filindex;
  MELEMENT *create_melem();

  fils = indsys->num_fils;
  meshes = indsys->num_mesh;
  Mlist = indsys->Mlist;
  Mtrans = indsys->Mtrans;

  /* clear it (should be garbage only) */
  for(i = 0; i < fils; i++)
    Mtrans[i] = NULL;

  for(j = 0; j < meshes; j++) {
    last_filindex = -1;
    for(m = Mlist[j]; m != NULL; m = m->mnext) {
      mt = make_melement(j, NULL, m->sign);
      Mtrans[m->filindex] = insert_in_list(mt,Mtrans[m->filindex]);
      if (last_filindex == m->filindex)
	printf("Internal Warning: Mesh %d contains the same filament multiple times\n",j);
      last_filindex = m->filindex;
    }
  }

#if 1==0
  /* this was the old inefficient way to from Mt */

  /* allocated for comparison */
  Mtrans_check = (MELEMENT **)MattAlloc(fils,sizeof(MELEMENT *));

  for(i = 0; i < fils; i++) {
    mt = &mtdum;
    /* scan through mesh list j looking for a filament i. (j,i) in M */
    for(j = 0; j < meshes; j++) {
      for(m = Mlist[j]; m->filindex < i && m->mnext != NULL; m = m->mnext)
	;
      if (m->filindex == i) {
	mt->mnext = make_melement(j, NULL, m->sign);
	mt = mt->mnext;
	if (m->mnext != NULL && m->mnext->filindex == i) {
	  for(count = 1; m->mnext->filindex == i; count++) {
	    m = m->mnext;
	    mt->mnext = make_melement(j, NULL, m->sign);
	    mt = mt->mnext;
	  }
	 viewprintf(stderr, "Internal Warning: Mesh %d contains the same filament %d times\n",j,count);
	} 
      }
    }
    mt->mnext = NULL;
    Mtrans_check[i] = mtdum.mnext;
  }

  /* check old and new ways */
  for(i = 0; i < fils; i++) {
    for(mt = Mtrans[i], mt2 = Mtrans_check[i]; mt != NULL && mt2 != NULL;
                  mt = mt->mnext, mt2 = mt2->mnext)
      if (mt->filindex != mt2->filindex || mt->sign != mt2->sign)
       viewprintf(stderr, "different: %d  %d %d\n",i,mt->filindex, mt2->filindex);
    if (mt != mt2) 
     viewprintf(stderr, "both not NULL:  mt: %u  mt2: %u\n", mt, mt2);
  }
#endif
      
}

compare_meshes(mesh1, mesh2)
MELEMENT *mesh1, *mesh2;
{

  while(mesh1 != NULL && mesh2 != NULL && mesh1->filindex == mesh2->filindex && mesh1->sign == mesh2->sign) {
    mesh1 = mesh1->mnext;
    mesh2 = mesh2->mnext;
  }

  if (mesh1 != NULL || mesh2 != NULL) {
    viewprintf(stderr, "meshes don't match!\n");
    FHExit(FH_GENERIC_ERROR);
  }
}

void cx_dumpMat_totextfile(fp, Z, rows, cols)
FILE *fp;
CX **Z;
int rows, cols;
{
  int i, j;

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) 
      fprintf(fp, "%13.6lg %+13.6lgj ", Z[i][j].real, Z[i][j].imag);
    fprintf(fp, "\n");
  }
  return;
}

void dumpMat_totextfile(fp, A, rows, cols)
FILE *fp;
double **A;
int rows, cols;
{
  int i, j;

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++) 
      fprintf(fp, "%13.6lg ", A[i][j]);
    fprintf(fp, "\n");
  }
  return;
}

void dumpVec_totextfile(fp2, Vec, size)
FILE *fp2;
double *Vec;
int size;
{
  dumpMat_totextfile(fp2, &Vec, 1, size);
}

void fillMrow(Mlist, mesh, Mrow)
MELEMENT **Mlist;
int mesh;
double *Mrow;
{
  MELEMENT *melem;

  if (mesh != 0)
    /* remove non-zeros from last call */
    for(melem = Mlist[mesh - 1]; melem != NULL; melem=melem->mnext)
      Mrow[melem->filindex] = 0;

  for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
    Mrow[melem->filindex] = melem->sign;
}

void dump_to_Ycond(fp, cond, indsys)
FILE *fp;
int cond;
SYS *indsys;
{
  static char fname[20], tempstr[5];
  int maxiters = indsys->opts->maxiters;

  sprintf(tempstr, "%d", cond);

  strcpy(fname,"Ycond");
  strcat(fname,tempstr);
  savecmplx(fp, fname, indsys->FinalY, indsys->num_extern, 
	    indsys->num_sub_extern);

  strcpy(fname,"resids");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resids, indsys->num_sub_extern, maxiters);

  strcpy(fname,"resid_real");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resid_real, indsys->num_sub_extern, maxiters);

  strcpy(fname,"resid_imag");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resid_imag, indsys->num_sub_extern, maxiters);

  strcpy(fname,"niters");
  strcat(fname,tempstr);
  saveCarray(fp, fname, &(indsys->niters), 1, indsys->num_sub_extern);

}

void saveCarray(fp, fname, Arr, rows, cols)
FILE *fp;
char *fname;
double **Arr;
int rows, cols;
{
  int i;
  int machine;

    machine = 1000;
#ifdef DEC
    machine = 0000;
#endif
#ifdef WIN32		// Enrico
    machine = 0000;
#endif

  for(i = 0; i < rows; i++) {
    savemat_mod(fp, machine+100, fname, rows, cols, 0, Arr[i], (double *)NULL,
		i, cols);
  }
}

int nnz_inM(Mlist, num_mesh)
MELEMENT **Mlist;
int num_mesh;
{
  int counter, i;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++)
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext)
      counter++;

  return counter;
}

void dump_M_to_text(fp, Mlist, num_mesh, nnz)
FILE *fp;
MELEMENT **Mlist;
int num_mesh, nnz;
{
  int counter, i;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++)
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext) {
      fprintf(fp, "%d\t%d\t%d\n", i+1, mesh->filindex+1, mesh->sign);
      counter++;
    }

  if (counter != nnz)
    viewprintf(stderr,"Internal Warning: nnz %d != counter %d\n",nnz,counter);
  
}

void dump_M_to_matlab(fp, Mlist, num_mesh, nnz, mname)
FILE *fp;
MELEMENT **Mlist;
int num_mesh, nnz;
char *mname;
{

  double onerow[3];
  int i, counter;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++) {
    onerow[0] = i + 1;
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext) {
      onerow[1] = mesh->filindex + 1;
      onerow[2] = mesh->sign;
      savemat_mod(fp, machine+100, mname, nnz, 3, 0, onerow,
		  (double *)NULL, counter++, 3);
    }
  }

  if (counter != nnz)
    viewprintf(stderr,"Internal Warning: nnz %d != counter %d\n",nnz,counter);
}  

/* this picks one node in each tree to be a ground node */
void pick_ground_nodes(indsys)
SYS *indsys;
{
  TREE *atree;
  char type;
  NODES *node;

  for(atree = indsys->trees; atree != NULL; atree = atree->next) {
    type = atree->loops->path->seg.type;
    if (type == NORMAL)
      node = getrealnode(((SEGMENT *)atree->loops->path->seg.segp)->node[0]);
    else
      node = getrealnode(((PSEUDO_SEG *)atree->loops->path->seg.segp)->node[0]);
      
    if (node->index != -1) {
      viewprintf(stderr,"huh? index == %d in pick_ground_node\n",node->index);
      FHExit(FH_GENERIC_ERROR);
    }
    node->index = 0;
  }
}

int pick_subset(portlist, indsys)
strlist *portlist;
SYS *indsys;
{
  strlist *oneport;
  EXTERNAL *ext;
  int counter;

  if (portlist == NULL) {
    for(ext = indsys->externals; ext != NULL; ext=ext->next)
      ext->col_Yindex = ext->Yindex;
    return indsys->num_extern;
  }
    
  for(ext = indsys->externals; ext != NULL; ext=ext->next)
    ext->col_Yindex = -1;

  counter = 0;
  for(oneport = portlist; oneport != NULL; oneport = oneport->next) {
    ext = get_external_from_portname(oneport->str,indsys);
    if (ext == NULL) {
      viewprintf(stderr,"Error: unknown portname %s\n",oneport->str);
      FHExit(FH_GENERIC_ERROR);
    }

    ext->col_Yindex = counter;
    counter++;
  }

  return counter;
}

/* concatenates so that s1 = s1 + s2 + s3 + s4 */
void concat4(s1,s2,s3,s4)
char *s1, *s2, *s3, *s4;
{
  s1[0] = '\0';
  strcat(s1,s2);
  strcat(s1,s3);
  strcat(s1,s4);
}

// Enrico, added explicit initialization of global and static vars
// to be able to re-use FastHenry as thread on multiple calls
void InitGlobAndStatVars(void)
{
	// initialize this module's global and static vars
	colmax = 0;
	tcol = NULL;
	trows = NULL; 

	forced = 0;
#ifdef DEC
	machine = 0000;
#elif defined WIN32
	machine = 0000;
#else
	machine = 1000;
#endif

	// initialize global and static vars in other modules 
	InitAddGndPlaneVars();
	InitBarnoldiVars();
	InitCalcpVars();
	InitDirectVars();
	InitFillMVars();
	InitFindPathsVars();
	InitGMResVars();
	InitJoelselfVars();
	InitMulDoVars();
	InitMulGlobalVars();
	InitMulMatsVars();
	InitMutualVars();
	InitPrecondVars();
	InitReadTreeVars();
	InitReadGeomVars();
	InitSetupComputePsiVars();
	InitUglierAllocVars();

}






