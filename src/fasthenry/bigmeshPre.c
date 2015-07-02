/* This preconditions one row. It is mostly a duplicate of indPrecond()
 code from olmulPrcond() */

#include "induct.h"
#define PARTMESH OFF

/* This near picks up only the hamming distance one cubes. */    
#define HNEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) + ABS((nbr)->k - (nk)) + ABS((nbr)->l - (nl))) <= 1)

/* This near picks up all 27 neighboring cubes. */
#define NEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) <= 1) && \
 (ABS((nbr)->k - (nk)) <= 1) && \
 (ABS((nbr)->l - (nl)) <= 1))

/* This near picks only the diagonal, for testing. */
#define DNEAR(nbr, nj, nk, nl) \
(((nbr)->j == (nj)) && \
 ((nbr)->k == (nk)) && \
 ((nbr)->l == (nl)) )

FILE *fp;


bigmeshPre(sys, indsys, w)
ssystem *sys;
SYS *indsys;
double w;
{

  cube *getcube(), *main_nc;
  static cube **nc_list;   /* an array of pointers to cubes */
  int msh, p;
  int numblocks, blocknum;
  static int *startrows;

  cube *nc, *nnbr, *nnnbr;
  static double **mat, **nmat;
  int i, j, k, l, m;
  int maxsize = 0, nsize, nnsize, nnnsize, *reorder;
  int nj, nk, nl, offset, noffset;
  int dindex, *nc_dummy, *nnbr_dummy, *nnnbr_dummy;
  static int *is_dummy;		/* local dummy flag vector, stays around */
  static int big_mat_size = 0;	/* size of previous mat */
  charge **nnnbr_pc, **nnbr_pc, **nc_pc, **mpc, *dp;
  surface *surf;
  double factor;

  /* FastHenry stuff */
  static CX **meshmat = NULL;
  static int meshmax = 0;       /* size of previous meshmat */
  static int *filcount = NULL;  /* number of fils per mesh for this cube and nbrs */
  static int *maxfilcount = NULL; /* max fils per mesh for any cube and nbrs */
  static int *indx = NULL;      /* index of real mesh number in meshmat */
  static int *meshnum = NULL;   /* mesh number corresponding to a row in meshmat */
         /* meshnum and indx should be inverses of each other (sort of). */
         /* i.e.  indx[meshnum[i]] == i. */
         /* meshnum[indx[i]] == i if i is one of the meshes in this cube */
  static int *fillist;      /* list of the filament numbers to which rows */
                            /* and cols of mat correspond */
  static int *findx;        /* For every filament, -1 if not in fillist, row */
                            /* number in mat if in fillist.  */
                            /* fillist and findx are inverses of each other */

  MELEMENT *getnext();

  int num_mesh = indsys->num_mesh;
  int num_fils = indsys->num_fils;
  int filnum, count;
  MELEMENT *mtran, *mtranj, *mtrani, *melem;
  MELEMENT **Mtrans = indsys->Mtrans;
  MELEMENT **Mlist = indsys->Mlist;
  PRE_ELEMENT **Precond = indsys->Precond;
  PRE_ELEMENT *pre, *prelast;
  double *R = indsys->R;
  int meshsize, realmrow;
  int counter, mrow, mcol;
  static int *is_in_nc;
  static DUPS *is_dup;
  static int *is_partial;
  int debug = 0;

  int xi, yi, zi;
  CX tempsum;
  charge *filchg;
  double length = sys->length;
  double minx = sys->minx, miny = sys->miny, minz = sys->minz;

  if (filcount == NULL) {
    CALLOC(filcount, num_mesh, int, ON, IND);
    CALLOC(is_partial, num_mesh, int, ON, IND);
    CALLOC(maxfilcount, num_mesh, int, ON, IND);
    CALLOC(indx, num_mesh, int, ON, IND);
    CALLOC(findx, num_fils, int, ON, IND);
    CALLOC(is_in_nc, num_mesh, int, ON, IND);
    CALLOC(nc_list, num_fils, cube*, ON, IND);
    CALLOC(startrows, num_fils + 1, int, ON, IND);

  }
  for (i = 0; i < num_mesh; i++)
    maxfilcount[i] = 0;


/* Figure out the max number of elements in any set of near cubes. */
  for(i = maxsize = 0; i < num_mesh; i++) 
    if (Precond[i] == NULL) {
      nsize = 0;
      for(j = 0; j < num_fils; j++)
	nc_list[j] = NULL;
      numblocks = 0;
      startrows[0] = 0;
      for(melem = Mlist[i]; melem != NULL; melem = melem->mnext) {
	filchg = melem->fil->pchg;
	main_nc = getcube(filchg, sys);
	if (main_nc == NULL) {
	  fprintf(stderr, "Hey, how come nc isn't in a cube?\n");
	  exit(1);
	}
	for(p = -1; p < main_nc->numnbrs; p++) {
	  if (p == -1)
	    nc = main_nc;
	  else
	    nc = main_nc->nbrs[p];
	  if (is_in_list(nc, nc_list, &blocknum) == FALSE) {  
	    /* nc isn't already in the list */
	    addtolist(nc, nc_list, &numblocks, startrows);
	    nsize += nc->directnumeles[0];
/*
	    nj = nc->j;
	    nk = nc->k;
	    nl = nc->l;
	    for(i=0; i < nc->numnbrs; i++) {
	      nnbr = nc->nbrs[i];
	      if(NEAR(nnbr, nj, nk, nl)) nsize += nnbr->directnumeles[0];
	    }
*/
	  }
	}
      }
      maxsize = MAX(nsize, maxsize);
    } /* end if Precond */

/* Allocate a matrix big enough for any set of 7. */
#if JACDBG == ON
  printf("max direct size =%d\n", maxsize);
#endif

    CALLOC(is_dummy, maxsize, int, ON, AMSC);
    MALLOC(mat, maxsize, double*, ON, AMSC);
    MALLOC(fillist, maxsize, int, ON, IND);   /* filament numbers  IND stuff */
    for(i=0; i < maxsize; i++) {
      MALLOC(mat[i], maxsize, double, ON, AMSC);
    }

      CALLOC(meshnum, maxsize, int, ON, IND);
      CALLOC(meshmat, maxsize, CX*, ON, IND);
      for(i = 0; i < maxsize; i++) 
	CALLOC(meshmat[i], maxsize, CX, ON, IND);
      CALLOC(is_dup, maxsize, DUPS, ON, IND);

/* Now go fill-in a matrix. */
  for(msh = 0; msh < num_mesh; msh++) {
    if (Precond[msh] != NULL) 
      continue;

    printf("mesh %d being worked on\n",msh);

    for(j = 0; j < num_fils; j++)
      nc_list[j] = NULL;
    numblocks = 0;
    startrows[0] = 0;

    for(i = 0; i < maxsize; i++)
      for(j = 0; j < maxsize; j++)
	mat[i][j] = 0;

    offset = 0;
    for(i = 0; i < num_mesh; i++) {
      filcount[i] = 0;
      is_in_nc[i] = 0;
      is_partial[i] = 0;
    }
    for(i = 0; i < num_fils; i++) 
      findx[i] = -1;
	
    for(melem = Mlist[msh]; melem != NULL; melem = melem->mnext) {
      filchg = melem->fil->pchg;
      main_nc = getcube(filchg, sys);
      if (main_nc == NULL) {
	fprintf(stderr, "Hey, how come filchg isn't in a cube?\n");
	exit(1);
      }
      for(p = -1; p == -1 /*p < main_nc->numnbrs */; p++) {
	if (p == -1) 
	  nc = main_nc;
	else
	  nc = main_nc->nbrs[p];
	if (is_in_list(nc, nc_list, &blocknum) == FALSE) {  
	  /* nc isn't already in the list */
	  addtolist(nc, nc_list, &numblocks, startrows);
	  nsize = nc->directnumeles[0];
	  nc_pc = nc->chgs;
	  nj = nc->j;
	  nk = nc->k;
	  nl = nc->l;
	  for(i = nsize - 1; i >= 0; i--) {
	    
	    filnum = fillist[offset + i] = nc_pc[i]->fil->filnumber;   /* IND stuff. 8/92 */
	    findx[nc_pc[i]->fil->filnumber] = offset + i;
	    /* find all the meshes that this filament is contained within */
	    for(mtran = indsys->Mtrans[filnum]; mtran != NULL; mtran=mtran->mnext) {
	      filcount[mtran->filindex]++;
	      is_in_nc[mtran->filindex] = 1;
	    }
	    
	    for(j = nsize - 1; j >= 0; j--) {
	      mat[offset + i][offset + j] = nc->directmats[0][i][j];
	    }
	  }
	  /* see if any of it's neighbors are in the list */
	  for(j = 0; j < nc->numnbrs; j++) {
	    nnbr = nc->nbrs[j];
	    if (is_in_list(nnbr, nc_list, &blocknum) == TRUE) {
	      nmat = nc->directmats[j+1];
	      nnsize = nc->directnumeles[j+1];
	      noffset = startrows[blocknum];
	      if (nnsize != startrows[blocknum +1] - startrows[blocknum]) {
		fprintf(stderr, "Hey, nnsize != start row difference \n");
		exit(1);
	      }
	      for(k = 0; k < nsize; k++)
		for(l = 0; l < nnsize; l++) 
		  mat[offset + k][noffset + l] = mat[noffset + l][offset + k] 
		    = nmat[k][l];
	    }
	  }
	  offset += nsize;
	}
      }
    }

    /* FastHenry stuff */

#if PARTMESH == OFF
    /* check to see if a mesh is only partly in the cube + neighbors */
    for(i = 0; i < num_mesh; i++)
      if (filcount[i] > 0) {
	count = 0;
	j = 0;
	for(melem = Mlist[i]; melem != NULL; melem = melem->mnext) {
	  count++;
	  j += melem->sign;
	}
	if (count != filcount[i]) {
	  if (count <= 2 || 1==1) {
	    filcount[i] = -1;   /* this is a partial mesh */
	    /*	  fprintf(stdout,"removed partial mesh #%d\n",i); */
	  }
	  else
	    is_partial[i] = TRUE;
	}
	else if (count > 2) 
	  printf("mesh number %d being worked on?\n", i);
      }
#endif

    /* count total number of meshes */
    meshsize = 0;
    for(i = 0; i < num_mesh; i++)
      if (filcount[i] > 0) meshsize++;

    if (meshsize > maxsize) {
      fprintf(stderr, "More meshes than filaments?????\n");
      exit(1);
    }

    /* fill indx and meshnum vectors */
    counter = 0;
    for(i = 0; i < num_mesh; i++) {
      if (filcount[i] > 0) {
	indx[i] = counter;
	meshnum[counter++] = i;
      }
      else {
	indx[i] = -1;
      }
    }
    if (counter != meshsize) {
      fprintf(stderr, "Hey, counter should equal meshsize\n");
      exit(1);
    }

    for(i = 0; i < meshsize; i++)
      for(j = 0; j < meshsize; j++)
	meshmat[i][j] = CXZERO;

    /* for each element in mat, determine it's contribution to */
    /* meshmat = M*mat*Mtrans */
    /* there may be a more efficient way to do this with some more */
    /* temporary storage. Like a temp matrix for mat*Mtran */

    for(i = 0; i < offset; i++)
      for(j = 0; j < offset; j++)
	for(mtranj=Mtrans[fillist[j]]; mtranj != NULL; mtranj = mtranj->mnext)
	  {
	    if (filcount[mtranj->filindex] > 0) {  /* mesh not partial */
	      mcol = indx[mtranj->filindex];
	      for(mtrani=Mtrans[fillist[i]]; mtrani!=NULL;mtrani=mtrani->mnext)
		{
		  if (filcount[mtrani->filindex] > 0) { /*mesh not partial */
		    mrow = indx[mtrani->filindex];
		    meshmat[mrow][mcol].imag 
		      += w*mtrani->sign*mat[i][j]*mtranj->sign;
		    if (i == j)
		      meshmat[mrow][mcol].real 
			+= mtrani->sign*R[fillist[i]]*mtranj->sign;
		  }
		}
	    }
	  }

    /* check if duplicate partial meshes (fills is_dup) */
    mark_dup_mesh(Mlist, meshnum, meshsize, is_dup, findx);

    if (debug == 1) {
      fp = fopen("chkinv.mat","w");
      if (fp == NULL) {printf("no open\n"); exit(1); }
      savecmplx(fp, "before", meshmat, meshsize, meshsize);
    }

    if (indsys->opts->debug == ON)
      fprintf(stdout, "Inverting a %d x %d matrix\n",meshsize,meshsize);

    /* now invert meshmat and skip duplicate rows and cols */
    cx_invert_dup(meshmat, meshsize, is_dup);

    if (debug == 1) {
      savecmplx(fp, "after", meshmat, meshsize, meshsize);
      fclose(fp);
    }

    /* add the rows to the preconditioner */

    i = indx[msh];
    if (i != -1) {
      realmrow = meshnum[i];
      if (filcount[ realmrow ] > maxfilcount[ realmrow ]) {
	maxfilcount[realmrow] = filcount[realmrow];
	if (Precond[realmrow] == NULL) {
	  CALLOC(Precond[realmrow], 1, PRE_ELEMENT, ON, IND);
	  Precond[realmrow]->next = NULL;
	}
	prelast = NULL;
	for(j = 0, pre = Precond[realmrow]; j < meshsize; j++) {
	  if (pre == NULL) {
	    CALLOC(pre, 1, PRE_ELEMENT, ON, IND);
	    pre->next = NULL;
	    if (prelast == NULL) {
	      fprintf(stderr, "Hey, prelast is null!\n");
	      exit(1);
	    }
	    prelast->next = pre;
	  }

	  pre->meshcol = meshnum[j];
	  if (is_dup[j].sign == 0)
	    pre->value = meshmat[i][j];
	  else
	    /* it's a duplicate, so use the duplicates inverse value. */
	    /* this effectively 'adds' the mesh currents of all duplicates */
	    cx_scalar_mult(pre->value, 
			   is_dup[j].sign, meshmat[i][is_dup[j].dup]);
	  prelast = pre;
	  pre = pre->next;
	}
      }
    }
    else {
      printf("Huh? how come msh isn't in matrix?\n");
      exit(1);
    }
  }
  

#if defined(DUMPPRE)
  dumpPrecond(Precond, num_mesh);
#endif


}


/* returns the cube that contains the charge */
cube *getcube(chg, sys)
charge *chg;
ssystem *sys;
{

  int xi, yi, zi;

  xi =  (chg->x - sys->minx) / sys->length;
  yi =  (chg->y - sys->miny) / sys->length;
  zi =  (chg->z - sys->minz) / sys->length;
  return sys->cubes[sys->depth][xi][yi][zi];

}

is_in_list(nc, list, blocknum)
cube *nc, **list;
int *blocknum;
{
  int i;

  i = 0;
  while(list[i] != NULL && list[i] != nc)
    i++;

  if (list[i] == nc) {
    *blocknum = i;
    return TRUE;
  }
  else {
    return FALSE;
  }
}

/* numblocks is the total number of blocks (cubes) in the matrix */
/* startrows is the row number where the block starts */

addtolist(nc, list, numblocks, startrows)
cube *nc, **list;
int *numblocks, *startrows;
{

  list[*numblocks] = nc;
  startrows[*numblocks + 1] = startrows[*numblocks] + nc->directnumeles[0];
  (*numblocks)++;

}
    
