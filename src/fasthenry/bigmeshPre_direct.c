/* this preconditions the 
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

bigmesh_direct(sys, indsys, w)
ssystem *sys;
SYS *indsys;
double w;
{
  cube *nc, *nnbr, *nnnbr;
  int nsize, nnsize;
  charge **nc_pc, **nnbr_pc;
  int meshmax = 0, *meshnum;
  CX **meshmat = NULL;
  int *filcount = NULL, *is_in_nc, *maxfilcount, *indx;
  int num_mesh = indsys->num_mesh;
  int filnum, i, j, k, nj, nk, nl;
  MELEMENT *mtran;
  MELEMENT **Mtrans = indsys->Mtrans;
  MELEMENT **Mlist = indsys->Mlist;
  PRE_ELEMENT **Precond = indsys->Precond;
  PRE_ELEMENT *pre, *prelast;
  int meshsize, realmrow;
  int counter;
  int debug = 0;
  CX **MtZM = indsys->MtZM;

  meshsize = count_tree_meshes(indsys->trees);

  if (meshsize > meshmax) {
    CALLOC(meshmat, meshsize + 10, CX*, ON, IND);
    for(i = 0; i < meshsize + 10; i++) 
      CALLOC(meshmat[i], meshsize + 10, CX, ON, IND);
    meshmax = meshsize + 10;
  }

  /* conveniently, all the bigmeshes (tree meshes) are at the front of Mlist */
  for(i = 0; i < meshsize; i++)
    for(j = 0; j < meshsize; j++)
      meshmat[i][j] = MtZM[i][j];

  if (indsys->opts->debug == ON) {
    fprintf(stdout, "For big meshes:\n");
    fprintf(stdout, "Inverting a %d x %d matrix\n",meshsize,meshsize);
  }

  /* now invert meshmat and skip duplicate rows and cols */
  cx_invert(meshmat, meshsize);

  if (debug == 1) {
    savecmplx(fp, "after", meshmat, meshsize, meshsize);
    fclose(fp);
  }

#if 1==0 /* the stupid way (adding to what is there) */
    /* add the rows (don't overwrite what's there) to the preconditioner */
    for(i = 0; i < meshsize; i++) {

      prelast = NULL;
      for(j = 0, pre = Precond[i]; j < meshsize; j++) {
	if (is_in_Precond(Precond[i], j, &prelast) == 0) {
	  CALLOC(pre, 1, PRE_ELEMENT, ON, IND);
	  pre->meshcol = j;
	  pre->value = meshmat[i][j];
	  if (prelast == NULL) {
	    pre->next = Precond[i];
	    Precond[i] = pre;
	  }
	  else {
	    pre->next = prelast->next;
	    prelast->next = pre;
	  }
	}
      }	  
    }
#endif

    /* add the rows to the preconditioner */
    for(i = 0; i < meshsize; i++) {

      Precond[i] = NULL;  /* a quick dumb solution */

      if (Precond[i] == NULL) {
        CALLOC(Precond[i], 1, PRE_ELEMENT, ON, IND);
        Precond[i]->next = NULL;
        }
      prelast = NULL;
      for(j = 0, pre = Precond[i]; j < meshsize; j++) {
        if (pre == NULL) {
          CALLOC(pre, 1, PRE_ELEMENT, ON, IND);
          pre->next = NULL;
          if (prelast == NULL) {
            fprintf(stderr, "Hey, prelast is null!\n");
            exit(1);
          }
          prelast->next = pre;
        }

        pre->meshcol = j;
        pre->value = meshmat[i][j];
        prelast = pre;
          pre = pre->next;
      }
    }

}  

is_in_Precond(prelist, col, last)
PRE_ELEMENT *prelist, **last;
int col;
{
  if (prelist == NULL) {
    *last = NULL;
    return 0;
  }
  
  if (prelist->meshcol == col) 
    return 1;
  else if (prelist->meshcol > col) {
    *last = NULL;
    return 0;
  }
  else {
    while(prelist->next != NULL && prelist->next->meshcol < col)
      prelist = prelist->next;

    if (prelist->next->meshcol == col)
      return 1;
    else {
      *last = prelist;
      return 0;
    }
  }
}



