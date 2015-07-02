/******************************************************************************
 *
 * File Name: barnoldi.c
 *	      (c) 1994 LMS, MATT
 * Author:    Luis Miguel Silveira
 * Revision:  23 Jun 95
 *
 * NAME
 *     barnoldi - creation of a reduced order model via the block arnoldi
 *
 * SYNOPSIS
 *       #include "induct.h"
 *
 * DESCRIPTION
 *       This is an implementation of the block-Arnoldi algorithm aimed at
 *       generating a reduced order state-space model of the system.
 *       Uses coordinate transform to make H symmetric.  All computation
 *       are done in the p-norm, defined by the matrix P.
 *       It fills the indsys->Ar,Br,Cr matrices with the model and
 *       returns the order of the model computed.
 * 
 * DIAGNOSTICS
 *
 * Modifications:
 *    9/9/96 - Added functions to create equivalent circuit - Mattan Kamon
 *
 *****************************************************************************/


#include "induct.h"
#include "sparse/spMatrix.h"

ArnoldiROM(B, C, D, P, size, numinp, numout, q_orig,
           matvec, indsys, sys, chglist)
  double **B, **C, **D;
  char **P;  /* really in sparse matrix form */
  int size, numinp, numout, q_orig;
  int (*matvec)();
  SYS *indsys;
  ssystem *sys;
  charge *chglist;
{
  static double **H = NULL, **V, **R1, **W, **X, **Z;
  int i, j, k, iter, bter;
  int s, q;
  
  /* s is the number of right-hand sides for Arnoldi */
  s = numinp;
  q = q_orig;
  
  if (H != NULL)
    fprintf(stderr, "ArnoldiROM called more than once, you should fix memory \
to make this efficient\n");

  /* allocate space for the matrices returned by Arnoldi */
  H = MatrixAlloc((q+1) * s, q * s, sizeof(double));
  V = MatrixAlloc(size, (q+1) * s, sizeof(double));
  R1 = MatrixAlloc(s, s, sizeof(double));
  /* for intermediate matrices here are a few tricks to save space */
  W = MatrixAlloc(size, s, sizeof(double));
  Z = MatrixAlloc(size, (q+1) * s, sizeof(double));
  X = MatrixAlloc(s, s, sizeof(double));

  ZeroMatrix(H, (q+1) * s, q*s);
  ZeroMatrix(R1, s, s);
  ZeroMatrix(X, s, s);

  if (size <= (q)*s) {
    q = size/s; 
    fprintf(stderr, "\n**Warning: Reduced order model would have higher order\
 (%d x %d = %d) than\n** the original system (%d)!  Will compute a %dth order\
 model instead. **\n\n",s,q_orig,s*q_orig,size, q*s);
  }

  /* QR-decompose the input matrix: [V(:, 1:s), R] = qr(B,0); */
  qr_P(B, V, R1, Z, size, s, 0, P);

  printf("Computing %d x %d matrix vector products for Reduced order model:\n",
         q, s);

  /* compute the Arnoldi iterates */
  for (iter = 0; iter < q; iter++) {
    printf("%d: ",iter+1);

    /* compute  W = A * Z(:, jms:js); use multipole expansions */
    matvec(W, sys, Z, chglist, indsys, size, s, iter * s);

    /* new basis block */
    for ( bter = (iter ? iter - 1 : 0); bter <= iter; bter++) {
      
      /* compute H(ims:is, jms:js) = V(:, ims:is)^T * W; */
      for (i = 0; i < s; i++) {
        for (j = 0; j < s; j++) {
          for (k = 0; k < size; k++) {
            /* note that V is transposed */
            H[bter * s + i][iter * s + j] +=
              Z[k][bter * s + i] * W[k][j];
          }
        }
      }

      /* compute W = W - V(:, ims:is) * H(ims:is, jms:js); */
      for (i = 0; i < size; i++) {
        for (j = 0; j < s; j++) {
          for (k = 0; k < s; k++) {
            W[i][j] -= V[i][bter * s + k] *
              H[bter * s + k][iter * s + j];
          }
        }
      }
    }

    /* compute [V(:, js1:j1s), X] = qr(W,0); */
    qr_P(W, V, X, Z, size, s, iter+1, P);
    for (i = 0; i < s; i++) {
      for (j = 0; j < s; j++) {
        H[(iter+1) * s + i][iter * s + j] = X[i][j];
      }
    }
  } 
  
  /* allocate and compute the reduced order model matrices */
  indsys->Ar = H;
  indsys->Br = MatrixAlloc(q * s, s, sizeof(double));
  /* compute Br = E1 * R1; top s by s block is R1, rest is zero */
  for (i = 0; i < s; i++) {
    for (j = 0; j < s; j++) {
      indsys->Br[i][j] = R1[i][j];
    }
  }
  /* Cr is trickier; well not really */
  indsys->Cr = MatrixAlloc(q * s, numout, sizeof(double));
  for (i = 0; i < (q * s); i++) {
    for (j = 0; j < numout; j++) {
      for (k = 0; k < size; k++) {
        /* note that Z is transposed */
        indsys->Cr[i][j] += Z[k][i] * C[k][j];
      }
    }
  }

  /* Dr is really hard ;-) */
  indsys->Dr = D;

  /* could clean some stuff but that's not the way FastHenry works... */
  
  return(q);
}


/******************************************************************************
 * qr()
 *
 * Arguments: Bmat - the matrix to be decomposed
 *            Qmat, Rmat - the Q and R blocks
 *            numlin, numcol - dimensions of Bmat
 *            block - Qmat is a row of matrices; this argument indicates
 *                    which block element to use
 * Returns: (nothing really)
 * Side-Effects: the Q and R blocks are filled in; the matrix is altered
 *
 * Description: QR skinny factorization using modified Gram-Schmidt
 *                       Bmat = Qmat Rmat
 *       Bmat has numlin rows and numcol columns
 *       Qmat has numlin rows but has K*numcol columns; block indicates what
 *       piece will be used to store the orthogonal matrix
 *
 *****************************************************************************/

qr(Bmat, Qmat, Rmat, numlin, numcol, block)
  double **Bmat, **Qmat, **Rmat;
  int numlin, numcol, block;
{
  int k, j, i;
  double normsq;
  
  for (k = 0; k < numcol; k++) {
    /* get the 2-norm of the k-th column */
    for (normsq = 0.0, i = 0; i < numlin; i++)
      normsq += (Bmat[i][k] * Bmat[i][k]);
    Rmat[k][k] = sqrt(normsq);
    /* normalize current vector */
    for (i = 0; i < numlin; i++)
      Qmat[i][block * numcol + k] = Bmat[i][k] / Rmat[k][k];
    /* orthogonalize */
    for (j = k+1; j < numcol; j++) {
      /* compute and store cross-product */
      Rmat[k][j] = 0.0;
      for (i = 0; i < numlin; i++)
        Rmat[k][j] += Qmat[i][block * numcol + k] * Bmat[i][j];
      for (i = 0; i < numlin; i++)
        Bmat[i][j] -= Qmat[i][block * numcol + k] * Rmat[k][j];
    }
  }
  
  /* yeah, the lint thing */
  return (1);
}

/******************************************************************************
 * qr_P()
 *
 * Arguments: Bmat - the matrix to be decomposed
 *            Qmat, Rmat - the Q and R blocks
 *            numlin, numcol - dimensions of Bmat
 *            block - Qmat is a row of matrices; this argument indicates
 *                    which block element to use
 *            P     - matrix to define P-norm
 * Returns: (nothing really)
 * Side-Effects: the Q and R blocks are filled in; the matrix is altered
 *
 * Description: QR skinny factorization using modified Gram-Schmidt
 *                       Bmat = Qmat Rmat
 *       Bmat has numlin rows and numcol columns
 *       Qmat has numlin rows but has K*numcol columns; block indicates what
 *       piece will be used to store the orthogonal matrix
 *
 *
 *       Q is orthogonal in the P-norm!  
 *
 *
 *****************************************************************************/

qr_P(Bmat, Qmat, Rmat, Z, numlin, numcol, block, P)
  double **Bmat, **Qmat, **Rmat, **Z;
  char *P;
  int numlin, numcol, block;
{
  int k, j, i;
  double normsq;
  static double *tcol = NULL;

  if (tcol == NULL)
    tcol = (double *) MattAlloc(numlin, sizeof(double));
  
  for (k = 0; k < numcol; k++) {
    /* copy column, tcol = w */
    for (i = 0; i < numlin; i++)
      tcol[i] = Bmat[i][k];
    /* Pw = z */
    spSolve(P, tcol, tcol);
    /* get the P-norm of the k-th column (wk' * zk) */
    for (normsq = 0.0, i = 0; i < numlin; i++)
      normsq += (tcol[i] * Bmat[i][k]);
    Rmat[k][k] = sqrt(normsq);
    /* normalize current vector, and copy tcol into Z */
    for (i = 0; i < numlin; i++) {
      Qmat[i][block * numcol + k] = Bmat[i][k] / Rmat[k][k];
      Z[i][block * numcol + k] = tcol[i]/Rmat[k][k];
    }
    /* orthogonalize */
    for (j = k+1; j < numcol; j++) {
      /* compute and store cross-product */
      Rmat[k][j] = 0.0;
      for (i = 0; i < numlin; i++)
        Rmat[k][j] += Z[i][block * numcol + k] * Bmat[i][j];
      for (i = 0; i < numlin; i++)
        Bmat[i][j] -= Qmat[i][block * numcol + k] * Rmat[k][j];
    }
  }
  
  /* yeah, the lint thing */
  return (1);
}


/******************************************************************************
 * dumpROM()
 *
 * Arguments: Ar, Br, Cr, Dr - the state-space representation of the
 *              reduced order model matrixes
 *            size - size of A (and the model really)
 *            numinp. numout - number of inputs and outputs
 * Returns: (nothing really)
 * Side-Effects: none
 *
 * Description: dumps the reduced-order model obatined via the block Arnoldi
 *      process to a file.  The file is written in textual matlab.  This can
 *      be improved...
 *****************************************************************************/

dumpROM(fp, Ar, Br, Cr, Dr, size, numinp, numout)
  FILE *fp;
  double **Ar, **Br, **Cr, **Dr;
  int size, numinp, numout;
{
  int i, j;

  /* put in a small header */
  fprintf(fp, "%%%% Reduced-order model state-space representation\n");
  fprintf(fp, "%%%% (read into matlab)\n\n");

  /* dump Ar */
  if (Ar != (double **)NULL) {  
    fprintf(fp, "Ar = [");
    for (i = 0; i < size; i++) {
      for (j = 0; j < size; j++) {
        fprintf(fp, "%13.6lg  ", Ar[i][j]);
        if (j%5 == 4)
          fprintf(fp, "...\n");
      }
      if (i < (size - 1))
        fprintf(fp, "\n");
      else
        fprintf(fp, "];\n\n");
    }
  }

  /* dump Br */
  if (Br != (double **)NULL) {
    fprintf(fp, "Br = [");
    for (i = 0; i < size; i++) {
      for (j = 0; j < numinp; j++) {
        fprintf(fp, "%13.6lg  ", Br[i][j]);
        if (j%5 == 4)
          fprintf(fp, "...\n");
      }
      if (i < (size - 1))
        fprintf(fp, "\n");
      else
        fprintf(fp, "];\n\n");
    }
  }

  /* dump Cr */
  if (Cr != (double **)NULL) {
    fprintf(fp, "Cr = [");
    for (i = 0; i < size; i++) {
      for (j = 0; j < numout; j++) {
        fprintf(fp, "%13.6lg  ", Cr[i][j]);
        if (j%5 == 4)
          fprintf(fp, "...\n");
      }
      if (i < (size - 1))
        fprintf(fp, "\n");
      else
        fprintf(fp, "];\n\n");
    }
  }

  /* dump Dr; check to see if there is really a direct term */
  if (Dr != (double **)NULL) {
    fprintf(fp, "Dr = [");
    for (i = 0; i < numout; i++) {
      for (j = 0; j < numinp; j++) {
        fprintf(fp, "%13.6lg  ", Dr[i][j]);
        if (j%5 == 4)
          fprintf(fp, "...\n");
      }
      if (i < (numout - 1))
        fprintf(fp, "\n");
      else
        fprintf(fp, "];\n\n");
    }
  }

  /* add extra stuff to compute residue form */
  if (Ar != (double **)NULL &&
      Br != (double **)NULL && 
      Cr != (double **)NULL) {
    fprintf(fp, "\n%%%% Matlab code to obtain pole-residue form\n\n");
    fprintf(fp, "%%%% first some miscellaneous stuff\n");
    fprintf(fp, "[n, n2] = size(Ar);\n");
    fprintf(fp, "[n, m]  = size(Br);\n");
    fprintf(fp, "[n, p]  = size(Cr);\n");
    fprintf(fp, "q = n / m;\n");
    fprintf(fp, "Rpoles = zeros(n, 1);\n");
    fprintf(fp, "Rresid = zeros(n, p*m);\n");
    fprintf(fp, "%%%% now the poles\n");
    fprintf(fp, "[Sr, EVr] = eig(Ar);\n");
    fprintf(fp, "Rpoles = 1 ./ diag(EVr);\n");
    fprintf(fp, "%%%% now for the zeros\n");
    fprintf(fp, "mu = Cr.' * Sr * inv(sparse(EVr));\n");
    fprintf(fp, "nu = inv(Sr) * Br;\n");
    fprintf(fp, "for (ii = 1:p)   %%%%  for each row\n");
    fprintf(fp, "  for (jj = 1:m)  %%%% for each column\n");
    fprintf(fp, "     Rresid(:, (jj-1)*m + ii) = -mu(ii,:).' .* nu(:,jj);\n");
    fprintf(fp, "  end\n");
    fprintf(fp, "end\n\n\n");
  }

  /* 'nuff said */
  return(1);
}

/******************************************************************************
 * dumpROMequiv_circuit()
 *
 * Arguments: Ar, Br, Cr, Dr - the state-space representation of the
 *              reduced order model matrixes
 *            size - size of A (and the model really)
 *            numinp. numout - number of inputs and outputs
 * Returns: (nothing really)
 * Side-Effects: none
 *
 * Description: creates a spice equivalent circuit composed of resistors,
 * capacitors and Voltage controlled current sources to represent the 
 * reduced-order model obatined via the block Arnoldi process.
 *
 * The assumed state space model is supposed to be of the form
 *    d/dt Ar x = x + Br v
 *            i = Cr^t x
 *
 * which is implemented by creating a subcircuit ROMequiv which uses
 * capacitors to represent Ar.  The state vector x will be the capacitor
 * nodes. As KCL equations:
 *
 *  - d/dt Ar x + x + Br v = 0
 * 
 * Thus we also need to represesnt x (=  I*x) which will be resistors
 *  of 1 ohm to ground and also Br*v which are VCCS's.
 *
 * The i = Cr^t x  terms are also VCCS's to be added to the current equations
 *   for the port nodes.
 *
 * I think the Br in FastHenry is missing a sign (should have used -inv(R)) so
 * I actually will swap its sign as noted in the comments below.
 *
 * Also note the minus sign in front of Ar
 *
 *   -Mattan Kamon  9/96
 *
 *****************************************************************************/

dumpROMequiv_circuit(fp, Ar, Br, Cr, Dr, size, numinp, numout, title, suffix,
                     indsys)
  FILE *fp;
  double **Ar, **Br, **Cr, **Dr;
  int size, numinp, numout;
  char *title;
  char *suffix;
  SYS *indsys;
{
  
  int i,j;
  EXTERNAL *ext;

  if (numinp != numout) {
    fprintf(stderr,"Error: Inputs and outputs assumed to be coincident. No equiv cicuited dumped\n");
    return;
  }
  
  fprintf(fp, "* Equivalent circuit for state-space model from analysis of FastHenry file:\n");
  fprintf(fp, "* %s\n",title);
  fprintf(fp, "* Fasthenry generates  A dx/dt = x + B Vin,  Iout = C^t x + D Vin\n");
  fprintf(fp, "* Note that regardless of what the original system modeled (say inductance),\n* this file will use only R's, C's, and VCCS's to represent the system\n");
  fprintf(fp, "* This is a %d-port model with the port nodes listed pairwise\n", numinp);
  fprintf(fp, "* Correspondence to FastHenry nodes:\n");

  for(ext = indsys->externals, j = 0;
      ext != NULL; ext = ext->next, j++) {
    if (ext->portname[0] == '\0')
      fprintf(fp,"*   Port %d:  %s  to  %s\n",ext->Yindex,ext->name1,ext->name2);
    else
      fprintf(fp,"*   Port %d:  %s  to  %s, port name: %s\n",
              ext->Yindex,ext->name1,ext->name2,ext->portname);
  }

  fprintf(fp, "* pn mn is a plus/minus port pair of nodes, where n is the port number\n");
  fprintf(fp, "\n.subckt ROMequiv%s ", suffix);
  
  for(i = 0; i < numinp; i++)
    fprintf(fp, " p%d m%d",i,i);
  fprintf(fp, "\n");

  fprintf(fp, "* The coefficient of x is the identity == 1 ohm resistors to ground.\n");
  for(i = 0; i < size; i++)
    fprintf(fp, "RH_%d_%d state%d 0 1\n", i,i,i);
  fprintf(fp, "\n");

  /* do B */
  if (Br != (double **)NULL) {
    fprintf(fp, "* The B matrix. VCCS dependent on port voltages\n");
    for (i = 0; i < size; i++) {
      for (j = 0; j < numinp; j++) {
        if (Br[i][j] != 0)
          /* i think there is a sign missing in the original definition
             of B, so this normally would be B, but instead i'll put -B */
          fprintf(fp, "GB_%d_%d  state%d 0 p%d m%d   %13.6lg\n", i, j, i, j, j,
                  -Br[i][j]);
      }
    }
  }

  /* do Cr */
  if (Cr != (double **)NULL) {
    fprintf(fp, "* The C matrix. VCCS dependent on state variable\n");
    for (i = 0; i < size; i++) {
      for (j = 0; j < numout; j++) {
        if (Cr[i][j] != 0.0)
          fprintf(fp, "GC_%d_%d   p%d m%d state%d 0   %13.6lg\n", i, j,j,j, i,
                  Cr[i][j]);
      }
    }
  }

  /* dump Dr; check to see if there is really a direct term */
  if (Dr != (double **)NULL) {
    fprintf(stderr, "D matrix not supported! (but would just be a bunch of VCCS's)\n");
    return;
  }

  /* do Ar */
  if (Ar != (double **)NULL) {  
    double sum;
    fprintf(fp, "* The A matrix.  Coefficient of dx/dt.  constructed from capacitors. A = A^t\n");
    for (i = 0; i < size; i++) {
      /* Assuming symmetric, we need only do above the diagonal */
      for (j = 0; j < i; j++) {
        if (Ar[i][j] != 0) {
          /* must choose value so stamp of capacitance matrix
             produces -Ar[i][j].  Off diag terms are -1*capacitance so -(-Ar)*/
          fprintf(fp, "CA_%d_%d state%d state%d  %13.6lg ic=0\n", i,j,i,j,
                  Ar[i][j]);
        }
      }

      sum = 0.0;
      /* sum the off diags */
      for(j = 0; j < size; j++)
        if (i != j)
          sum += Ar[i][j];

      /* do the diagonal.  Add in the off-diag parts and set the sum to 
         ground */
      fprintf(fp, "CA_%d_%d state%d 0  %13.6lg ic=0\n", i,i,i,-(Ar[i][i]+sum));

    }
  }

  fprintf(fp, ".ends ROMequiv%s\n", suffix);
  
}
/******************************************************************************
 * dumpROMbin()
 *
 * Arguments: A, B, C, D - the state-space representation of the system
 *            size - size of A (and the model really)
 *            numinp. numout - number of inputs and outputs
 * Returns: (nothing really)
 * Side-Effects: none
 *
 * Description: dumps the original system representation to a file
 *****************************************************************************/

dumpROMbin(fp, A, B, C, D, size, numinp, numout)
  FILE *fp;
  double **A, **B, **C, **D;
  int size, numinp, numout;
{
  int i, j, realsize;
  static double *temp = (double *) NULL;
  int machine;

  machine = 1000;
#ifdef DEC
  machine = 0000;
#endif

  if (temp == (double *) NULL) {
    realsize = MAX(size, numinp);
    realsize = MAX(realsize, numout);
    temp = (double *) calloc(realsize, sizeof(double));
  }

  if (A != (double **)NULL) {
    for(j = 0; j < size; j++) {
      for(i = 0; i < size; i++)
        temp[i] = A[i][j];
      savemat_mod(fp, machine+100, "Asys", size, size, 0, temp,
                  (double *)NULL, j, size);
    }
  }
  if (B != (double **)NULL) {
    for(j = 0; j < numinp; j++) {
      for(i = 0; i < size; i++)
        temp[i] = B[i][j];
      savemat_mod(fp, machine+100, "Bsys", size, numinp, 0, temp, 
                  (double *)NULL, j, size);
    }
  }
  if (C != (double **)NULL) {
    for(j = 0; j < numout; j++) {
      for(i = 0; i < size; i++)
        temp[i] = C[i][j];
      savemat_mod(fp, machine+100, "Csys", size, numout, 0, temp, 
                  (double *)NULL, j, size);
    }
  }
  if (D != (double **)NULL) {
    for(j = 0; j < numinp; j++) {
      for(i = 0; i < numout; i++) 
        temp[i] = D[i][j];
      savemat_mod(fp, machine+100, "Dsys", numout, numinp, 0, temp, 
                  (double *)NULL, j, numout);
    }
  }
  return(1);
}


/******************************************************************************
 * createMRMt()
 *
 * Arguments: MRMt_Ptr - pointer to the M * R * Mt matrix
 *            indsys - fasthenry sys structure
 * Returns: (nothing really)
 * Side-Effects: MRMt gets filled in
 *
 * Description: Creates and fills the sparse matrix MRMt used in the model
 *      generation
 *****************************************************************************/

createMRMt(MRMt_Ptr, indsys)
  char **MRMt_Ptr;
  SYS *indsys;
{
  SEGMENT *seg;
  double valR;
  double *R = indsys->R;
  int i, num_fils;
  MELEMENT *mtranj, *mtrani;
  MELEMENT **Mtrans = indsys->Mtrans;
  double *elem;
  int filnum, err;
  char *MRMt;

  MRMt = (char *) spCreate(indsys->num_mesh, 0, &err);
  if (err != 0) {
    fprintf(stderr,"Couldn't create sparse matrix, err %d\n", err);
    exit(1);
  }

  for (seg = indsys->segment; seg != NULL; seg = seg->next) {
    num_fils = seg->num_fils;
    for (i = 0; i < num_fils; i++) {
      filnum = seg->filaments[i].filnumber;
      valR = R[filnum];
      for (mtranj=Mtrans[filnum]; mtranj != NULL; mtranj = mtranj->mnext) {
	  for (mtrani=Mtrans[filnum]; mtrani != NULL; mtrani = mtrani->mnext) {
            if (mtrani->filindex + 1 > indsys->num_mesh ||
                mtranj->filindex + 1 > indsys->num_mesh) {
              fprintf(stderr, "Indexing into matrix MRMt out of bounds\n");
              exit(1);
            }
	    *(elem = (double *) spGetElement(MRMt, mtrani->filindex+1,
                                             mtranj->filindex+1))
              += mtrani->sign * valR * mtranj->sign;
	  }
      }
    }
  }
  err = spFactor(MRMt);
  if (err != 0) {
    fprintf(stderr,"Error on factor: %d\n",err);
    exit(1);
  }

  /* return the value */
  *MRMt_Ptr = MRMt;
  
  return(1);
}


/******************************************************************************
 * createMRMtinvMLMt()
 *
 * Arguments: MRMtinvMLMt_Ptr - pointer to (M * R * Mt)^-1 (M * L * Mt) matrix
 *            indsys - fasthenry sys structure
 *            MRMt - the M * R * Mt matrix
 * Returns: (nothing really)
 * Side-Effects: MRMtinvMLMt gets filled in with
 *
 * Description: Creates and fills the matrix MRMtinvMLMt used in the model
 *      generation.  
 *****************************************************************************/

createMRMtinvMLMt(MRMtinvMLMt_Ptr, indsys, MRMt)
  double ***MRMtinvMLMt_Ptr;
  SYS *indsys;
  char *MRMt;
{
  double **MRMtinvMLMt;
  int m, n, p;
  double tempsum, tempR, tempsumR;
  static double *tcol = NULL, *temp = NULL;  /* temp storage for extra rows */
  int i, j, k, mesh, mesh2, nodeindx;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  MELEMENT *mt, *mt2;    /* elements of M transpose */
  double **M, **L;
  int rows, cols, num_mesh;
  MELEMENT **Mlist, **Mt;

  MRMtinvMLMt = (double **)
    MatrixAlloc(indsys->num_mesh, indsys->num_mesh, sizeof(double));

  L = indsys->Z;
  M = indsys->M;
  nfils = rows = indsys->num_fils;
  nmesh = cols = num_mesh = indsys->num_mesh;
  Mlist = indsys->Mlist;
  Mt = indsys->Mtrans;

  /* allocate a temporary column for L*Mt */
  tcol = (double *) MattAlloc(rows, sizeof(double));
  temp = (double *) MattAlloc(nmesh, sizeof(double));

  /* this does L*(Mt)_i where (Mt)_i is a single column (row) of Mt (M) and
     saves it in the temp space, tcol */
  for (mesh = 0; mesh < num_mesh; mesh++) {
    for (j = 0; j< nfils; j++)
      tcol[j] = 0;
    /* note, this next set of nested loops could be reversed, but I think
       this might be more efficient for pipelining, ? */
    for (melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
      for (j = 0; j < nfils; j++)
	tcol[j] += L[j][melem->filindex] * melem->sign;
    for (mesh2 = 0; mesh2 < num_mesh; mesh2++) {
      tempsum = 0;
      for (melem2 = Mlist[mesh2]; melem2 != NULL; melem2 = melem2->mnext)
	tempsum += melem2->sign * tcol[melem2->filindex];
      temp[mesh2] = tempsum;
    }
    /* multiply on the left by -(M*R*Mt)^-1 and
     * put back the vector in the returning matrix; so actually solve
     *      -(M*R*Mr) X = Vs
     */
    spSolve(MRMt, temp, temp);
    for (mesh2 = 0; mesh2 < num_mesh; mesh2++) {
      MRMtinvMLMt[mesh2][mesh] = -temp[mesh2];
    }
  }
  *MRMtinvMLMt_Ptr = MRMtinvMLMt;

  return(1);
}





/******************************************************************************
 * realComputePsi()
 *
 * Arguments: Prod - the result, i.e. the product of two "matrices"
 *            sys - fastcap sys structure
 *            B - the second matrix
 *            chglist - info on how to setup the charges
 *            indsys - fasthenry sys structure
 *            size, numRhs - dimensions for C
 *            initcol - initial column of B to use in products
 * Returns: (nothing)
 * Side-Effects: none other than Prod gets filled in
 *
 * Description: Matrix - matrix multiply
 *                Prod = -L * B
 *       however A is not really formed and the product A * B is done
 *       one column of B at a time using the multipole algorithm for the
 *       matrix-vector product;  this is done by calling Fastcap's
 *       computePsi which is a matrix vector product; therefore numRHS calls
 *       will be needed, in fact 3 times that since it will be called for
 *       each coordinate direction.
 *       Note that the inversion shown is not actually performed.  Instead
 *       multiple solves are done on the matrix MRMt which BTW is already
 *       factored
 *****************************************************************************/

realComputePsi(Prod, sys, B, chglist, indsys, size, numRHS, initcol)
  double **Prod, **B;
  ssystem *sys;
  charge *chglist;
  SYS *indsys;
  int size, numRHS, initcol;
{
  static double *Ib = NULL, *Vb = NULL;
  static double *Im = NULL, *Vs = NULL;
  int branches;
  MELEMENT **Mtrans, **Mlist;
  MELEMENT *mtemp;
  int currRHS;
  int ind_opcnt_mult = 0, ind_opcnt_real = 0;
  double *q, *p;
  extern double dirtime;
  charge *chg;
  int i, j;

  branches = indsys->num_fils;
  Mtrans = indsys->Mtrans;
  Mlist = indsys->Mlist;

  if (Ib == NULL) {
    Ib = (double *) MattAlloc(branches, sizeof(double));
    Im = (double *) MattAlloc(size, sizeof(double));
    Vb = (double *) MattAlloc(branches, sizeof(double));
    Vs = (double *) MattAlloc(size, sizeof(double));
  }

  /* pick up one column of B */
  for (currRHS = 0; currRHS < numRHS; currRHS++) {

    /* fill it in in Im */
    for (i = 0; i < size; i++)
      Im[i] = B[i][initcol + currRHS];
    for (i = 0; i < branches; i++)
      Vb[i] = 0.0;

    q = sys->q;
    p = sys->p;
    ASSERT(size == indsys->num_mesh);

    /* Do all of the non-direct parts first; start with Ib = Mtrans*Im */
    for (i = 0; i < branches; i++) {
      Ib[i] = 0.0;
      for (mtemp = Mtrans[i]; mtemp != NULL; mtemp = mtemp->mnext) {
        if (mtemp->sign == 1) 
          Ib[i] += Im[mtemp->filindex];
        else
          Ib[i] -= Im[mtemp->filindex];
      }
    }

    /* Evaluate (M*)L*Mt*Im = (M*)L*Ib using the multipole algorithm */
    sys->DirectEval = FALSE;
    for (i = 0; i < 3; i++) {  /* for each of the coordinate directions */
      for (chg = chglist; chg != NULL; chg = chg->next) {
        /* fill the pseudo-charge vector */
        q[chg->index] = Ib[chg->fil->filnumber] * chg->fil->lenvect[i];
#if OPCNT == ON
        ind_opcnt_mult++;
#endif
      }
      computePsi(sys, q, p, branches, chglist);
      for (chg = chglist; chg != NULL; chg = chg->next) {
        /* add potential due to i direction */
        Vb[chg->fil->filnumber] +=
          p[chg->index] * chg->fil->lenvect[i] * MUOVER4PI;
#if OPCNT == ON
        ind_opcnt_mult++;
#endif
      }
    }
  
    /* do the final direct parts M*Vb = M*(L*Ib) = M*(L*Mt*Im) */
    sys->DirectEval = TRUE;
    
    for (i = 1; i <= branches; i++)
      p[i] = 0;
    for (chg = chglist; chg != NULL; chg = chg->next) 
      /* fill the pseudo-charge vector */
      q[chg->index] = Ib[chg->fil->filnumber];
    
    /* starttimer; */
    mulDirect(sys);
    mulEval(sys);
    /* stoptimer; */
    dirtime += dtime;
    
    for (chg = chglist; chg != NULL; chg = chg->next) {
      /* add potential due to i direction */
      Vb[chg->fil->filnumber] += p[chg->index];
    }
    
    /* do Vs = M * Vb */
    for (i = 0; i < size; i++) {
      Vs[i] = 0.0;
      for (mtemp = Mlist[i]; mtemp != NULL; mtemp = mtemp->mnext)
        if (mtemp->sign == 1) 
          Vs[i] += Vb[mtemp->filindex];
        else
          Vs[i] -= Vb[mtemp->filindex];
    }

    /* multiply by -1 and store in result */
    for (i = 0; i < size; i++) {
      Prod[i][currRHS] = -Vs[i];
    }

    printf("%d ",currRHS);
    fflush(stdout);
  }
  printf("\n");

  return(1);
}


/******************************************************************************
 * realMatVect()
 *
 * Arguments: Prod - the result, i.e. the product of two "matrices"
 *            sys - fastcap sys structure
 *            B - the second matrix
 *            chglist - info on how to setup the charges
 *            indsys - fasthenry sys structure
 *            size, numRhs - dimensions for C
 *            initcol - initial column of B to use in products
 * Returns: (nothing)
 * Side-Effects: none other than Prod gets filled in
 *
 * Description: Matrix - matrix multiply
 *                Prod = -MLMt * B
 *        in this situation the matrix has indeed been compuded explicitely
 *****************************************************************************/

realMatVect(Prod, sys, B, chglist, indsys, size, numRHS, initcol)
  double **Prod, **B;
  ssystem *sys;
  charge *chglist;
  SYS *indsys;
  int size, numRHS, initcol;
{
  int i, j, currRHS;
  double temp;
  CX **Zm = indsys->MtZM;

  /* pick up one column of B */
  for (currRHS = 0; currRHS < numRHS; currRHS++) {
    printf("%d ",currRHS);
    fflush(stdout);
    /* now do the multiplication */
    for (i = 0; i < size; i++) {
      Prod[i][currRHS] = 0.0;
      for (j = 0; j < size; j++) {
        temp = -Zm[i][j].imag * B[j][initcol + currRHS];
        Prod[i][currRHS] += temp;
      }
    }
  }

  printf("\n");
  return(1);
}


/******************************************************************************
 * printRowCol()
 *
 * Arguments: mat - matrix to be debugged
 *            rowcol - row or column or diagonal
 *            rownum, colnum - number of row or column to print (or diagonal)
 *            size - length of vector to print
 * Returns: (nothing really)
 * Side-Effects: none
 *
 * Description: Use for debugging purposes; given a matrix and a flag prints
 *         either a row or column of that matrix to file "DUMPS".  Each
 *         subsequent call overwrites the file.
 *****************************************************************************/

printRowCol(mat, rowcol, rownum, colnum, size)
  double **mat;
  int rowcol, rownum, colnum;
{
  int i;
  FILE *localfp;

  localfp = fopen("DUMPS", "w");
  if (localfp == NULL) {
    printf("Couldn't open DUMPS file\n");
    exit(1);
  }

  switch(rowcol) {
  default:
  case 0:
    /* print column colnum */
    for (i = 0; i < size; i++)
      fprintf(localfp, "%13.6lg  ", mat[i][colnum]);
    break;
  case 1:
    /* print row rownum */
    for (i = 0; i < size; i++)
      fprintf(localfp, "%13.6lg  ", mat[rownum][i]);
    break;
  case 2:
    /* print diagonal rownum */
    for (i = 0; i < size; i++)
      fprintf(localfp, "%13.6lg  ", mat[i][i + rownum]);
    break;

  }
  fclose(localfp);

  return(1);
}
  
/* This computes the product M*L*Mt and stores it in imag(indsys->MtZM) */
/*  It is just formMZMt() without the real part */
formMLMt(indsys)
SYS *indsys;
{
  int m,n,p;
  double tempsum, tempR, tempsumR;
  static double *tcol = NULL;   /* temporary storage for extra rows */
  int i,j, k, mesh, mesh2, nodeindx;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  MELEMENT *mt, *mt2;    /* elements of M transpose */
  double **M, **L, *R;
  CX **Zm, *tempZ;
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
    fprintf(stderr, "Uh oh, more meshes than filaments, I'm confused\n");
    exit(1);
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

}

/* sets everything to zero in a matrix */
ZeroMatrix(A, rows, cols)
  double **A;
  int rows, cols;
{
  int i, j;

  for(i = 0; i < rows; i++)
    for(j = 0; j < cols; j++)
      A[i][j] = 0.0;
}
