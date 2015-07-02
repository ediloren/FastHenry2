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
 
*//* Produces a fasthenry input file of filaments on a UNIFORM
   plane discretization */

#include <stdio.h>
#include <math.h>

FILE *f_uni, *f_nonuni;
char *MattAlloc();

main(argc, argv)
int argc;
char *argv[];
{
  int i,j;
  double min_len, x_len, y_len, thickness, z, x_cur;
  double *y, *x_horiz, *x_vert;
  double *vert_width;
  double horiz_len, horiz_width;
  double extra_x, extra_y;
  double r0, ratio;
  int numequiv;

  x_len = 21;
  y_len = 21;
  thickness = 0.01;
  z = 0;
  r0 = 0.05;
  ratio = 1.0/3.0;

  if (sscanf(argv[1], "%lf", &r0) != 1) {
    printf("need one integer argument:  The radius of contacts\n");
    exit(1);
  }

  numequiv = 200;

  /* f_nonuni = fopen("fh_nonuni.inp","w"); */
  f_nonuni = stdout;

  fprintf(f_nonuni,"* Nonuniform test for r=%lg\n",r0);
  fprintf(f_nonuni,".units m\n.default sigma=5.8e7\n\n");

  extra_x = (x_len - 1)/2.0;
  extra_y = (y_len - 1)/2.0;

  fprintf(f_nonuni,"g1 x1=%lg y1=%lg z1=%lg x2=%lg y2=%lg z2=%lg x3=%lg y3=%lg z3=%lg\n",
	  -extra_x, -extra_y, z, 1 + extra_x, -extra_y, z, 1 + extra_x, 
	  1 + extra_y, z);

  fprintf(f_nonuni,"+ file=NONE thick=%lg \n",
	  thickness);
  fprintf(f_nonuni,"+ contact decay_rect (0,0, 0, %lg, %lg, %lg, %lg, 3, 3)\n",
	  r0, r0, r0*ratio, r0*ratio);
  fprintf(f_nonuni,"+ contact decay_rect (1,1, 0, %lg, %lg, %lg, %lg, 3, 3)\n",
	  r0, r0, r0*ratio, r0*ratio);


  /* equiv points in a circle */
  for(i = 0; i < numequiv; i++) {
    fprintf(f_nonuni,"+ n_in_%d (%lg,%lg,%lg)\n",i, r0*cos(i*2*M_PI/numequiv), 
	    r0*sin(i*2*M_PI/numequiv),z);
    fprintf(f_nonuni,"+ n_out_%d (%lg,%lg,%lg)\n",
	    i, 1 + r0*cos(i*2*M_PI/numequiv), 
	    1 + r0*sin(i*2*M_PI/numequiv),z);
  }

  fprintf(f_nonuni,"\n.equiv n_in ");
  for(i = 0; i < numequiv; i++) {
    fprintf(f_nonuni, "n_in_%d ",i);
    if (i%10 == 0)
      fprintf(f_nonuni, "\n+ ");
  }

  fprintf(f_nonuni,"\n.equiv n_out ");
  for(i = 0; i < numequiv; i++) {
    fprintf(f_nonuni, "n_out_%d ",i);
    if (i%10 == 0)
      fprintf(f_nonuni, "\n+ ");
  }


  fprintf(f_nonuni,"\n\n");
  fprintf(f_nonuni,".external n_in n_out\n");

  fprintf(f_nonuni,".freq fmin=0 fmax=1e9 ndec=1\n");
  fprintf(f_nonuni,".end\n");

  fclose(f_nonuni);
}

double **MatrixAlloc(rows, cols, size)
int rows, cols, size;
{

  double **temp;
  int i;

  temp = (double **)MattAlloc(rows,sizeof(double *));
  if (temp == NULL) {
        printf("not enough space for matrix allocation\n");
        exit(1);
      }

  for(i = 0; i < rows; i++) 
    temp[i] = (double *)MattAlloc(cols,size);

  if (temp[rows - 1] == NULL) {
        printf("not enough space for matrix allocation\n");
        exit(1);
      }
  return(temp);
}

char* MattAlloc(number, size)
int number, size;
{

  char *blah;

  blah = (char *)malloc(number*size);

  if (blah == NULL) {
    fprintf(stderr, "MattAlloc: Couldn't get space. Needed %d\n",number*size);
    exit(1);
  }

  return blah;
}
