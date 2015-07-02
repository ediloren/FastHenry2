#include <stdio.h>
#include <math.h>
#include "induct.h"


/* this writes a file of faces of the segments suitable to be
   read in by keith's program for generating postscript images */

writefastcap(fname, shading_name, indsys)
char *fname, *shading_name;
SYS *indsys;
{
  SEGMENT *seg;
  NODES *node0, *node1;
  double x,y ,z, xbotr, ybotr, zbotr;
  int i, j, k;
  double hx,hy,hz,mag,wx,wy,wz, height, width;
  FILE *fp, *shade_fp;
  SEGMENT *segment = indsys->segment;
  GROUNDPLANE *gp;
  int thin;

  fp = fopen(fname,"w");
  if (fp == NULL) {
    printf("can't open zbuffile to write \n");
    exit(1);
  }

  shade_fp = fopen(shading_name,"w");
  if (fp == NULL) {
    printf("can't open shading file to write \n");
    exit(1);
  }

  /* choose shades for segments and save in seg.is_deleted */
  assign_shades(indsys);

  fprintf(fp, "0 This is an inductance picture \n");
  for(i = 0, seg = segment; seg != NULL; seg = seg->next, i++) {
    if (seg->type == GPTYPE && indsys->opts->gp_draw == OFF)
      continue;   /* skip gp segs */
    node0 = seg->node[0];
    node1 = seg->node[1];
    height = seg->height;
    width = seg->width;
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
    hx = -wy*(node1->z - node0->z) + (node1->y - node0->y)*wz;
    hy = -wz*(node1->x - node0->x) + (node1->z - node0->z)*wx;
    hz = -wx*(node1->y - node0->y) + (node1->x - node0->x)*wy;
    mag = sqrt(hx*hx + hy*hy + hz*hz);
    hx = hx/mag;
    hy = hy/mag;
    hz = hz/mag;

    if (seg->type == GPTYPE && indsys->opts->gp_draw == THIN)
      thin = TRUE;
    else
      thin = FALSE;


#if 1==0  
    /* first node end */
    fprintf(fp, "Q %d ",i);
    x = node0->x - 0.5 * wx * width + 0.5 * hx * height;
    y = node0->y - 0.5 * wy * width + 0.5 * hy * height;
    z = node0->z - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);

    x = node0->x - 0.5 * wx * width - 0.5 * hx * height;
    y = node0->y - 0.5 * wy * width - 0.5 * hy * height;
    z = node0->z - 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    x = node0->x + 0.5 * wx * width - 0.5 * hx * height;
    y = node0->y + 0.5 * wy * width - 0.5 * hy * height;
    z = node0->z + 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = node0->x + 0.5 * wx * width + 0.5 * hx * height;
    y = node0->y + 0.5 * wy * width + 0.5 * hy * height;
    z = node0->z + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    fprintf(fp,"\n");
    /* write out shade of grey to color */
    fprintf(shade_fp,"%d\n",seg->is_deleted);

    /* end at other node */
    fprintf(fp, "Q %d ",i);
    x = node1->x - 0.5 * wx * width + 0.5 * hx * height;
    y = node1->y - 0.5 * wy * width + 0.5 * hy * height;
    z = node1->z - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);
  
    x = node1->x - 0.5 * wx * width - 0.5 * hx * height;
    y = node1->y - 0.5 * wy * width - 0.5 * hy * height;
    z = node1->z - 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 

    x = node1->x + 0.5 * wx * width - 0.5 * hx * height;
    y = node1->y + 0.5 * wy * width - 0.5 * hy * height;
    z = node1->z + 0.5 * wz * width - 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = node1->x + 0.5 * wx * width + 0.5 * hx * height;
    y = node1->y + 0.5 * wy * width + 0.5 * hy * height;
    z = node1->z + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
    fprintf(fp,"\n");
    /* write out shade of grey to color */
    fprintf(shade_fp,"%d\n",seg->is_deleted);
#endif

    if (thin == FALSE) {
      /* left side */
      fprintf(fp, "Q %d ",i);
      x = node0->x - 0.5 * wx * width + 0.5 * hx * height;
      y = node0->y - 0.5 * wy * width + 0.5 * hy * height;
      z = node0->z - 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z);

      x = node0->x - 0.5 * wx * width - 0.5 * hx * height;
      y = node0->y - 0.5 * wy * width - 0.5 * hy * height;
      z = node0->z - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 

      x = node1->x - 0.5 * wx * width - 0.5 * hx * height;
      y = node1->y - 0.5 * wy * width - 0.5 * hy * height;
      z = node1->z - 0.5 * wz * width - 0.5 * hz * height;
  
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = node1->x - 0.5 * wx * width + 0.5 * hx * height;
      y = node1->y - 0.5 * wy * width + 0.5 * hy * height;
      z = node1->z - 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z);

      fprintf(fp,"\n");
      /* write out shade of grey to color */
      fprintf(shade_fp,"%d\n",seg->is_deleted);

    /* right side */
      fprintf(fp, "Q %d ",i);
      x = node0->x + 0.5 * wx * width - 0.5 * hx * height;
      y = node0->y + 0.5 * wy * width - 0.5 * hy * height;
      z = node0->z + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
      x = node0->x + 0.5 * wx * width + 0.5 * hx * height;
      y = node0->y + 0.5 * wy * width + 0.5 * hy * height;
      z = node0->z + 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
      x = node1->x + 0.5 * wx * width + 0.5 * hx * height;
      y = node1->y + 0.5 * wy * width + 0.5 * hy * height;
      z = node1->z + 0.5 * wz * width + 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
      x = node1->x + 0.5 * wx * width - 0.5 * hx * height;
      y = node1->y + 0.5 * wy * width - 0.5 * hy * height;
      z = node1->z + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      fprintf(fp, "\n");
      /* write out shade of grey to color */
      fprintf(shade_fp,"%d\n",seg->is_deleted);
    }
  
    /* top */
    fprintf(fp, "Q %d ",i);
    x = node0->x - 0.5 * wx * width + 0.5 * hx * height;
    y = node0->y - 0.5 * wy * width + 0.5 * hy * height;
    z = node0->z - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);
   
    x = node0->x + 0.5 * wx * width + 0.5 * hx * height;
    y = node0->y + 0.5 * wy * width + 0.5 * hy * height;
    z = node0->z + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
   
    x = node1->x + 0.5 * wx * width + 0.5 * hx * height;
    y = node1->y + 0.5 * wy * width + 0.5 * hy * height;
    z = node1->z + 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z); 
    x = node1->x - 0.5 * wx * width + 0.5 * hx * height;
    y = node1->y - 0.5 * wy * width + 0.5 * hy * height;
    z = node1->z - 0.5 * wz * width + 0.5 * hz * height;
    fprintf(fp,"%lg %lg %lg ",x,y,z);

    fprintf(fp,"\n");
    /* write out shade of grey to color */
    fprintf(shade_fp,"%d\n",seg->is_deleted);

    if (thin == FALSE) {
      /* bottom */
      fprintf(fp, "Q %d ",i);

      x = node0->x - 0.5 * wx * width - 0.5 * hx * height;
      y = node0->y - 0.5 * wy * width - 0.5 * hy * height;
      z = node0->z - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 

      x = node0->x + 0.5 * wx * width - 0.5 * hx * height;
      y = node0->y + 0.5 * wy * width - 0.5 * hy * height;
      z = node0->z + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = node1->x + 0.5 * wx * width - 0.5 * hx * height;
      y = node1->y + 0.5 * wy * width - 0.5 * hy * height;
      z = node1->z + 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
      x = node1->x - 0.5 * wx * width - 0.5 * hx * height;
      y = node1->y - 0.5 * wy * width - 0.5 * hy * height;
      z = node1->z - 0.5 * wz * width - 0.5 * hz * height;
      fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
      fprintf(fp,"\n");
      /* write out shade of grey to color */
      fprintf(shade_fp,"%d\n",seg->is_deleted);
    }
  }
  
  if (indsys->opts->gp_draw == OFF) 
   for(gp = indsys->planes; gp != NULL; gp = gp->next) {
      
    unit_cross_product(gp->x[0]-gp->x[1], gp->y[0]-gp->y[1], gp->z[0]-gp->z[1],
		       gp->x[2]-gp->x[1], gp->y[2]-gp->y[1], gp->z[2]-gp->z[1],
		       &wx, &wy, &wz);
    if (is_nonuni_gp(gp))
      width = gp->thick;
    else
      width = gp->segs1[0][0]->height;
    
    for(j = 0; j < 3; j++) {
      k = j+1;
      fprintf(fp, "Q %d ",i + j);
      make_between(fp, gp,wx,wy,wz,width,j,k);
      /* shade it white */
      fprintf(shade_fp,"%d\n",0);
    }
    k = 0;
    j = 3;
    fprintf(fp, "Q %d ",i + j);
    make_between(fp, gp, wx,wy,wz,width,j,k);
    fprintf(shade_fp,"%d\n",0);
    i += 4;
   }    

  /* clear is_deleted */
  clear_marks(indsys);

  fclose(fp);
  fclose(shade_fp);
}

make_between(fp, gp, wx,wy,wz,width,j,k)
GROUNDPLANE *gp;
double wx,wy,wz,width;
int j,k;
FILE *fp;
{
  double x,y,z;

  x = gp->x[j] - 0.5 * wx * width;
  y = gp->y[j] - 0.5 * wy * width;
  z = gp->z[j] - 0.5 * wz * width;
  fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
  x = gp->x[j] + 0.5 * wx * width;
  y = gp->y[j] + 0.5 * wy * width;
  z = gp->z[j] + 0.5 * wz * width;
  fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
  x = gp->x[k] + 0.5 * wx * width;
  y = gp->y[k] + 0.5 * wy * width;
  z = gp->z[k] + 0.5 * wz * width;
  fprintf(fp,"%lg %lg %lg ",x,y,z); 

  x = gp->x[k] - 0.5 * wx * width;
  y = gp->y[k] - 0.5 * wy * width;
  z = gp->z[k] - 0.5 * wz * width;
  fprintf(fp,"%lg %lg %lg ",x,y,z); 
  
  fprintf(fp,"\n");
}

unit_cross_product(x1,y1,z1,x2,y2,z2,cx, cy,cz)
double x1,y1,z1,x2,y2,z2,*cx, *cy,*cz;
{
  double magc;

  *cx = (y1*z2 - y2*z1);
  *cy = (z1*x2 - z2*x1);
  *cz = (x1*y2 - x2*y1);

  magc = mag(*cx, *cy, *cz);

  if (magc/mag(x1,y1,z1) > 1e-13) {
    /* normalize */
    *cx = *cx/magc;
    *cy = *cy/magc;
    *cz = *cz/magc;
  }
}
    
/* attempts to choose shades for each segment to make seeing the condutors 
   easier.  To see shades, run with the zbuf -q options */

assign_shades(indsys)
SYS *indsys;
{
  int_list *elem;
  EXTERNAL *ext;
  MELEMENT **Mlist = indsys->Mlist, *m;
  int greylev, nodes1, nodes2, i, j;
  GROUNDPLANE *gp;

  /* clear is_deleted */
  clear_marks(indsys);

  greylev = 0;
  /* color each set of external meshes a different shade of grey */
  for(ext = indsys->externals; ext != NULL; ext = ext->next) {
    greylev++;

    /* for each mesh in this set of ext (usually only one) */
    for(elem = ext->indices; elem != NULL; elem = elem->next)

      /* for each filament in this mesh */
      for(m = Mlist[elem->index]; m != NULL; m = m->mnext)
        
        m->fil->segm->is_deleted = greylev;
  }

  greylev += 10;

  /* now do ground planes */
  for(gp = indsys->planes; gp != NULL; gp = gp->next) {
    greylev++;

    nodes2 = gp->num_nodes2;
    nodes1 = gp->num_nodes1;
    
    for(i = 0; i < nodes2; i++)
      for(j = 0; j < (nodes1 - 1); j++)
        if (gp->segs1[j][i] != NULL)  /* it's not part of a hole */
          gp->segs1[j][i]->is_deleted = greylev;

    for(i = 0; i < (nodes2 - 1); i++)
      for(j = 0; j < nodes1; j++)
        if (gp->segs2[j][i] != NULL)  /* it's not part of a hole */
          gp->segs2[j][i]->is_deleted = greylev;
  }
}
