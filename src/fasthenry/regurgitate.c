/* This regurgitates the input file to the output with some 
   translations and replications if desired */
/* It outputs in SI units only */

#include "induct.h"

regurgitate(indsys)
SYS *indsys;
{
  /* functions to shift coordinates */
  void translate(), reflect_x(), reflect_y(), reflect_origin(), do_nothing();

  /* spit out geometries and .externals */
  spit(indsys, do_nothing, "");

/*
  spit(indsys, reflect_x, "_X");
  spit(indsys, reflect_y, "_Y");
  spit(indsys, reflect_origin, "_XY");
  */

  do_end_stuff(indsys);
}

do_end_stuff(indsys)
SYS *indsys;
{
  fprintf(stdout, ".freq fmin=%lg fmax=%lg ndec=%lg\n",indsys->fmin, 
          indsys->fmax, 1.0/indsys->logofstep);

  fprintf(stdout, ".end\n");

}



static double delta_x, delta_y, delta_z;

/* set amount of translation for translate() */
set_translate(x, y, z)
double x, y, z;
{
  delta_x = x;
  delta_y = y;
  delta_z = z;
}

void translate(x,y,z,new_x, new_y, new_z)
double x,y,z,*new_x, *new_y, *new_z;
{
  *new_x = x + delta_x;
  *new_y = y + delta_y;
  *new_z = z + delta_z;
}

/* reflect about x axis (negate y) */
void reflect_x(x,y,z,new_x, new_y, new_z)
double x,y,z,*new_x, *new_y, *new_z;
{
  *new_x = x;
  *new_y = -y;
  *new_z = z;
}

/* reflect about y axis (negate x) */

void reflect_y(x,y,z,new_x, new_y, new_z)
double x,y,z,*new_x, *new_y, *new_z;
{
  *new_x = -x;
  *new_y = y;
  *new_z = z;
}

/* reflect about origin (negate x and y) */

void reflect_origin(x,y,z,new_x, new_y, new_z)
double x,y,z,*new_x, *new_y, *new_z;
{
  *new_x = -x;
  *new_y = -y;
  *new_z = z;
}

void do_nothing(x,y,z,new_x, new_y, new_z)
double x,y,z,*new_x, *new_y, *new_z;
{
  *new_x = x;
  *new_y = y;
  *new_z = z;
}


/* spit out nodes and segments and .externals */

spit(indsys, new_coords, suffix)
SYS *indsys;
void (*new_coords)();
char *suffix;
{
  
  NODES *node;
  NODELIST *nodel;
  SEGMENT *seg;
  GROUNDPLANE *gp;
  EXTERNAL *ext;
  int i;
  double x,y,z;

  /* do nodes */
  fprintf(stdout, "* NODES for %s\n",suffix);
  for(node  = indsys->nodes; node != NULL; node = node->next) {
    if (node->type == NORMAL) {
      new_coords(node->x, node->y, node->z, &x, &y, &z);
      fprintf(stdout, "%s%s x=%lg y=%lg z=%lg\n", node->name, suffix, x, y, z);
    }
  }

  /*do segments */
  fprintf(stdout, "* Segments for %s\n",suffix);
  for(seg = indsys->segment; seg != NULL; seg = seg->next) {
    if (seg->type == NORMAL) {
      fprintf(stdout, "%s%s %s%s %s%s w=%lg h=%lg nhinc=%d nwinc=%d rw=%lg rh=%lg sigma=%lg",
              seg->name, suffix, seg->node[0]->name, suffix, 
              seg->node[1]->name, suffix, seg->width, 
              seg->height, seg->hinc, 
              seg->winc, seg->r_width, seg->r_height, seg->sigma);
      if (seg->widthdir != NULL)
        fprintf(stdout, " wx=%lg wy=%lg wz=%lg \n", seg->widthdir[XX],
                seg->widthdir[YY], seg->widthdir[ZZ]);
      else
        fprintf(stdout, "\n");
    }
  }

  /*do planes */
  fprintf(stdout, "* Planes for %s\n",suffix);
  for(gp = indsys->planes; gp != NULL; gp = gp->next) {

    if (is_nonuni_gp(gp)) {
      fprintf(stdout, "Nonuniform planes not supported for regurgitation at this time!\n");
      continue;
    }

    fprintf(stdout, "%s%s \n",gp->name, suffix);
    i = 0;
    new_coords(gp->x[i], gp->y[i], gp->z[i], &x, &y, &z);
    fprintf(stdout, "+ x1=%lg y1=%lg z1=%lg\n",x, y, z);

    i = 1;
    new_coords(gp->x[i], gp->y[i], gp->z[i], &x, &y, &z);
    fprintf(stdout, "+ x2=%lg y2=%lg z2=%lg\n",x, y, z);

    i = 2;
    new_coords(gp->x[i], gp->y[i], gp->z[i], &x, &y, &z);
    fprintf(stdout, "+ x3=%lg y3=%lg z3=%lg\n",x, y, z);

    fprintf(stdout, "+ seg1=%d seg2=%d\n", gp->seg1, gp->seg2);
    fprintf(stdout, "+ thick=%lg ", gp->thick);
    fprintf(stdout, "segwid1=%lg segwid2=%lg ", gp->segwid1, gp->segwid2);
    fprintf(stdout, "sigma=%lg\n", gp->sigma);

    for(nodel = gp->usernode_coords; nodel != NULL; nodel=nodel->next) {
      fprintf(stdout, "+ %s (%lg, %lg, %lg)\n", nodel->name, nodel->x, 
              nodel->y, nodel->z);
    }

    if (gp->list_of_holes != NULL)
      fprintf(stdout, "Holes cannot be regurgitated at this time!\n");
      
  }

  /* do .equivs */
  fprintf(stdout, "* .equivs for %s\n",suffix);
  for(node  = indsys->nodes; node != NULL; node = node->next) {
    if (node->type == NORMAL && node != getrealnode(node)) {
      fprintf(stdout, ".equiv %s%s  %s%s\n", node->name, suffix, 
              getrealnode(node)->name, suffix);
    }
  }

  /* do .externals */
  
  for(ext = indsys->externals; ext != NULL; ext = ext->next) {
    fprintf(stdout, ".external %s%s %s%s ",
            getrealnode(ext->source->node[0])->name, suffix,
            getrealnode(ext->source->node[1])->name, suffix);
    if (strcmp(ext->portname,"") != 0)
      fprintf(stdout, "%s%s\n",ext->portname, suffix);
    else
      fprintf(stdout, "\n");
  }
}
    
                
      
  
