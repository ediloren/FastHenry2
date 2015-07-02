/* this is a function from Joel Phillips for generating matlab matrices
   for visualization of a fastcap structure */
/*  It produces the file panels.mat */

#include "mulGlobal.h"

dump_struct(chglist,qv)
charge *chglist; 
double *qv; 
{

  charge *cp; 
  double *x, *y, *z, *q; 
  int index; 
  FILE *fp, *fopen(); 
  int type; 
  int size, i; 
  extern char *title, *ps_file_base, *in_file_name;
  char *fname, *concat3();
  
  fname = concat3(ps_file_base,".","mat");

  fprintf(stdout,"Dumping structure to Matlab file %s\n", fname);

  size = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    size++; 
  }
  if ( (fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "Can't open %s\n",fname); 
    exit(1);
  }

  x = (double *) calloc( 4*size, sizeof(double)); 
  y = (double *) calloc( 4*size, sizeof(double)); 
  z = (double *) calloc( 4*size, sizeof(double)); 
  q = (double *) calloc( 4*size, sizeof(double)); 

  /* do the triangles */ 
  index = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    if (cp->shape == 3) {
      for (i=0; i<cp->shape; i++) {
        x[3*index+i] = cp->corner[i][0] * cp->X[0] + cp->corner[i][1]*cp->Y[0] + 
          cp->corner[i][2] * cp->Z[0] + cp->x; 

        y[3*index+i] = cp->corner[i][0] * cp->X[1] + cp->corner[i][1]*cp->Y[1] + 
          cp->corner[i][2] * cp->Z[1] + cp->y; 

        z[3*index+i] = cp->corner[i][0] * cp->X[2] + cp->corner[i][1]*cp->Y[2] + 
          cp->corner[i][2] * cp->Z[2] + cp->z; 
        
      }
      if (qv != NULL)
        q[index] = qv[cp->index]; 
      else 
        q[index] = 0.0; 
      index++; 
    }
  }
  
  
#ifdef sun 
  type = 1000; 
#else 
  type = 0000; 
#endif 
  if (index > 0) {
    savemat(fp, type, "xt", 3, index, 0, x, NULL);
    savemat(fp, type, "yt", 3, index, 0, y, NULL);
    savemat(fp, type, "zt", 3, index, 0, z, NULL);
    savemat(fp, type, "qt", 3, index, 0, q, NULL);
  }
  
  /* now the quads */ 

  index = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    if (cp->shape == 4) {
      for (i=0; i<cp->shape; i++) {
        x[4*index+i] = cp->corner[i][0] * cp->X[0] + cp->corner[i][1]*cp->Y[0] + 
          cp->corner[i][2] * cp->Z[0] + cp->x; 

        y[4*index+i] = cp->corner[i][0] * cp->X[1] + cp->corner[i][1]*cp->Y[1] + 
          cp->corner[i][2] * cp->Z[1] + cp->y; 

        z[4*index+i] = cp->corner[i][0] * cp->X[2] + cp->corner[i][1]*cp->Y[2] + 
          cp->corner[i][2] * cp->Z[2] + cp->z; 
        
      }
      if (qv != NULL) 
        q[index] = qv[cp->index]; 
      else 
        q[index] = 0.0; 
      index++; 
    }
  }
  if (index > 0) {
    savemat(fp, type, "xq", 4, index, 0, x, NULL);
    savemat(fp, type, "yq", 4, index, 0, y, NULL);
    savemat(fp, type, "zq", 4, index, 0, z, NULL);
    savemat(fp, type, "qq", 4, index, 0, q, NULL);
  }

  fclose(fp); 

}

