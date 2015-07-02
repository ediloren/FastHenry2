/*
 * savemat - C language routine to save a matrix in a MAT-file.
 *
 * We recommend that you use this routine, and its companion loadmat.c,
 * for all writing and reading of MAT-files.  These routines implement
 * "access methods" for MAT-files.  By using these routines, instead
 * of writing your own code to directly access the MAT-file format,
 * you will be unaffected by any changes that may be made to the MAT-file
 * structure at a future date.
 *
 * Here is an example that uses 'savemat' to save two matrices to disk,
 * the second of which is complex:
 *
 *	FILE *fp;
 *	double xyz[1000], ar[1000], ai[1000];
 *	fp = fopen("foo.mat","wb");
 *	savemat(fp, 2000, "xyz", 2, 3, 0, xyz, (double *)0);
 *	savemat(fp, 2000, "a", 5, 5, 1, ar, ai);
 *      fclose(fp);
 *
 * Author J.N. Little 11-3-86
 * Revised 7-23-91 to support ANSI-C
 */
#include <stdio.h>
// Enrico
#include <string.h>

#ifdef ALPHA
typedef struct {
     int type;   /* type */
     int mrows;  /* row dimension */
     int ncols;  /* column dimension */
     int imagf;  /* flag indicating imag part */
     int namlen; /* name length (including NULL) */
} Fmatrix;
#else
typedef struct {
     long type;   /* type */
     long mrows;  /* row dimension */
     long ncols;  /* column dimension */
     long imagf;  /* flag indicating imag part */
     long namlen; /* name length (including NULL) */
} Fmatrix;
#endif

#ifdef __STDC__
void savemat(FILE *fp, int type, char *pname, int mrows, int ncols,
             int imagf, double *preal, double *pimag)
#else
void savemat(fp, type, pname, mrows, ncols, imagf, preal, pimag)
FILE *fp;       /* File pointer */
int type;       /* Type flag: Normally 0 for PC, 1000 for Sun, Mac, */
		/* Apollo, and other Motorola format, */
		/* 2000 for VAX D-float, 3000 for VAX G-float, and */
		/* 4000 for CRAY */
		/* Add 1 for text variables, 2 for sparse matrices */
		/* See LOAD in reference section of guide for more info.*/
char *pname;    /* pointer to matrix name */
int mrows;      /* row dimension */
int ncols;      /* column dimension */
int imagf;	/* imaginary flag */
double *preal;  /* pointer to real data */
double *pimag;  /* pointer to imag data */
#endif
{
	Fmatrix x;
	int mn;

	x.type = type;
	x.mrows = mrows;
	x.ncols = ncols;
	x.imagf = imagf;
	x.namlen = strlen(pname) + 1;
	mn = x.mrows * x.ncols;

	fwrite(&x, sizeof(Fmatrix), 1, fp);
	fwrite(pname, sizeof(char), (int)x.namlen, fp);
	fwrite(preal, sizeof(double), mn, fp);
	if (imagf) {
	     fwrite(pimag, sizeof(double), mn, fp);
	}
}

/*
  MODIFIED version of above: added wr_flag to allow multiple writes
  to same matrix
  wr_flag = 0 => open, print header (like old matlab setup)
  wr_flag = 1 => update, print without header
  Doesn't work for complex.  For complex, full real part must precede
  full imaginary.
*/
#ifdef __STDC__
void savemat_mod(FILE *fp, int type, char *pname, int mrows, int ncols,
		 int imagf, double *preal, double *pimag, int wr_flag, int mn)
#else
void savemat_mod(fp, type, pname, mrows, ncols, imagf, preal, pimag,
		 wr_flag, mn)
FILE *fp;       /* File pointer */
int type;       /* Type flag: Normally 0 for PC, 1000 for Sun, Mac, */
		/* Apollo, and other Motorola format, */
		/* 2000 for VAX D-float, 3000 for VAX G-float, and */
		/* 4000 for CRAY */
		/* Add 1 for text variables, 2 for sparse matrices */
		/* See LOAD in reference section of guide for more info.*/
char *pname;    /* pointer to matrix name */
int mrows;      /* row dimension */
int ncols;      /* column dimension */
int imagf;	/* imaginary flag */
double *preal;  /* pointer to real data */
double *pimag;  /* pointer to imag data */
int wr_flag;			/* 0 for open, 1 to add to matrix */
int mn;				/* real #entries, this dump only */
#endif
{
	Fmatrix x;

	if(wr_flag == 0) {
	  x.type = type;
	  x.mrows = mrows;
	  x.ncols = ncols;
	  x.imagf = imagf;
	  x.namlen = strlen(pname) + 1;

	  fwrite(&x, sizeof(Fmatrix), 1, fp);
	  fwrite(pname, sizeof(char), (int)x.namlen, fp);
	}
	fwrite(preal, sizeof(double), mn, fp);
/*	if (imagf) {
	     fwrite(pimag, sizeof(double), mn, fp);
	}
*/
}
