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
 
*/
/* The following function is an example of how to read in the list of matrices
   contained in Zc.mat */
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
#include "cmplx.h"
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
  int i,j;

  if (argc != 2) {
    printf("error: give one argument, and make it a filename to read\n");
    exit(1);
  }

  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    fprintf(stderr, "No open\n");
    exit(1);
  }

  Zc = ReadZc(fp);

  for(oneZ = Zc; oneZ != NULL; oneZ = oneZ->next) {
    printf("freq = %lg\n",oneZ->freq);
    for(i = 0; i < oneZ->size1; i++) {
      printf("Row %d: ",i);
      for(j = 0; j < oneZ->size2; j++)
	printf("%lg%+lgj ",oneZ->Z[i][j].real,
	       oneZ->Z[i][j].imag/(2*PI*oneZ->freq));
      printf("\n");
    }
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
      printf("Reading Frequency %lg\n", freq);
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
	  printf("Unexpected end of file\n");
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
	    printf("Couldn't read j off of imaginary part\n");
	    printerr(freq, i, j);
	  }
	  ptr += 1;
	}
      }
    }
    else
      printf("Not part of any matrix: %s", line);
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
