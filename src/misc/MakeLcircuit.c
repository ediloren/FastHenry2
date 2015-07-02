/* This makes an inductor circuit out of a single impedance matrix */

/* the function ReadZc returns a pointer to a linked list called 'Zlist'
   which contains all relevant information about a matrix.  
*/
/* main() is an example of how to call ReadZc and access what it returns */
/* To include this in other code, remove main() and also #include "cmplx.h"
   and the typedef for Zlist */

/* Also, ReadZc was written quickly, so it doesn't do some vital
   things like check to make sure that malloc() returns nonNULL.  I
   suggest tailoring it to your needs */


#include <stdio.h>
#include <math.h>
#include "cmplx.h"
// Enrico
#include <stdlib.h>

#define LINEMAX 100000
#define PI 3.141592654

/* a linked list of the matrices in the given file */
typedef struct _Zlist {
  CX **Z;
  double freq;
  int size1, size2;
  struct _Zlist *next;
} Zlist;

Zlist *ReadZc();

main(argc, argv)
int argc;
char *argv[];
{
  FILE *fp;
  Zlist *Zc, *oneZ;
  int i,j, size, intnode, nodecount;
  double freq;
  CX **Z;
  double eps = 1e-6;

  if (argc != 2) {
    fprintf(stderr,"error: give one argument, and make it a filename to read\n");
    exit(1);
  }

  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    fprintf(stderr, "No open\n");
    exit(1);
  }

  Zc = ReadZc(fp);

  if (Zc->next != NULL)
    fprintf(stderr,"** WARNING: More than one frequency in input. Using freq=%lg\n", 
            Zc->freq);

  Z = Zc->Z;

  size = Zc->size1;
  freq = Zc->freq;

  /* use the first Zc */
  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++)
      /* divide out freq */
      Z[i][j].imag = Z[i][j].imag/(2*PI*freq);

  /* symmetrize it */
  for(i = 1; i < size; i++)
    for(j = 0; j < i; j++) {
      Z[i][j].real = Z[j][i].real = (Z[i][j].real + Z[j][i].real)/2.0;
      Z[i][j].imag = Z[j][i].imag = (Z[i][j].imag + Z[j][i].imag)/2.0;
    }

  intnode = 2;
  /* self terms */
  for(i = 0; i < size; i++) {
/*    printf("rZ_%d %d int_node%d %lg\n",i,2*i+1,intnode,Z[i][i].real);*/
    printf("LZ_%d %d int_node%d_%d %lg\n",i,i*2+1,i,intnode,
	   Z[i][i].imag);
/*    printf("Vam_%d int_node%d_%d %d dc=0v\n",i,i,size+intnode,i*2+2);*/
  }
  
  /* mutual inductance */
  for(i = 1; i < size; i++)
    for(j = 0; j < i; j++)
      printf("KZ_%d_%d LZ_%d LZ_%d %lg\n",i,j,i,j,
	     Z[i][j].imag/sqrt(Z[i][i].imag * Z[j][j].imag));
      
  intnode = 2;
  /* mutual resistance */
  for(i = 0; i < size; i++) {
    nodecount = 0;
    for(j = 0; j < size; j++)
      if (fabs(Z[i][j].real/Z[i][i].real) > eps 
	  || fabs(Z[i][j].real/Z[j][j].real) > eps) {
	if (i == j) 
	  printf("RZ_%d_%d int_node%d_%d int_node%d_%d %lg\n",i,j,i,
		 intnode+nodecount,i,intnode+nodecount+1, Z[i][j].real);
	else
	  printf("HZ_%d_%d int_node%d_%d int_node%d_%d Vam_%d %lg\n",i,j,i,
		 intnode+nodecount,i,intnode+nodecount+1, j, Z[i][j].real);
	nodecount++;
      }
    printf("Vam_%d int_node%d_%d %d dc=0v\n",i,i,intnode+nodecount,i*2+2);
  }

}

Zlist *ReadZc(fp)
FILE *fp;
{
  static char line[LINEMAX];
  char *ptr;
  int i, j, size1, size2, skip;
  double freq;
  Zlist *Zc, *oneZ;
  double real, imag;

  Zc = NULL;

  while(fgets(line, LINEMAX, fp) != NULL) {
    if (sscanf(line, "Impedance matrix for frequency = %lg %d x %d",
	       &freq, &size1, &size2) == 3) {
      oneZ = (Zlist *)malloc(sizeof(Zlist));
      fprintf(stderr,"Reading Frequency %lg\n", freq);
      oneZ->freq = freq;
      oneZ->size1 = size1;
      oneZ->size2 = size2;
      oneZ->next = Zc;  /* tack new matrix to front of list */
      Zc = oneZ;
      oneZ->Z = (CX **)malloc(size1*sizeof(CX *));
      for(i = 0; i < size1; i++)
	oneZ->Z[i] = (CX *)malloc(size2*sizeof(CX));

      for(i = 0; i < size1; i++) {
	if (fgets(line, LINEMAX, fp) == NULL) {
	  fprintf(stderr,"Unexpected end of file\n");
	  printerr(freq,i,0);
	}
	ptr = line;
	for(j = 0; j < size2; j++) {
	  if (sscanf(ptr, "%lf%n", &real, &skip) != 1)
	    printerr(freq,size1, size2);
	  else
	    oneZ->Z[i][j].real = real;
	  ptr += skip;

	  if (sscanf(ptr, "%lf%n", &imag, &skip) != 1)
	    printerr(freq, i, j);
	  else
	    oneZ->Z[i][j].imag = imag;
	  ptr += skip;

	  if (ptr[0] != 'j') {
	    fprintf(stderr,"Couldn't read j off of imaginary part\n");
	    printerr(freq, i, j);
	  }
	  ptr += 1;
	}
      }
    }
    else
      fprintf(stderr,"Not part of any matrix: %s", line);
  }

  return Zc;
}

printerr(freq, row, col)
double freq;
int row, col;
{
  fprintf(stderr, "Error reading value in matrix for frequency %lg at row %d, col %d\n", freq, row, col);
  exit(1);
}
