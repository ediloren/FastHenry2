/* Calculates the minimum distance between two filaments */
/* Among other things, like Gaussian Quadrature weights */
#include "induct.h"

/* this is where the Gaussian Quadrature are defined */
double **Gweight, **Gpoint;

double dist_between();
double min_endpt_sep();
double dist_betw_pt_and_fil();

/* The line defined by fil1 is r = s1 + t1*D1 where s1 is the vector
   (fil1->x[0], y[0], z[0]) and D1 is (x[1] - x[0], y[1] - y[0], z[1] - z[0]).
    t1 is a scalar from 0 to 1.  Similarly for fil2 */

double dist_betw_fils(fil1, fil2, parallel)
FILAMENT *fil1, *fil2;
int *parallel;
{
  double k1,k2,k3,k4,c1,c2,c3,c4,a1,a2,a3,b1,b2,b3;
  double x1,y1,z1,x2,y2,z2;
  double D1[3], D2[3], s1[3], s2[3], e[3], *D, t1, t2, s1ms2[3], s1me[3];
  double D1D1, D1D2, D2D2, D1s1s2, D2s1s2;
  double tmp1, tmp2;

  s1[XX] = fil1->x[0];
  s1[YY] = fil1->y[0];
  s1[ZZ] = fil1->z[0];
  s2[XX] = fil2->x[0];
  s2[YY] = fil2->y[0];
  s2[ZZ] = fil2->z[0];

  s1ms2[XX] = s1[XX] - s2[XX];
  s1ms2[YY] = s1[YY] - s2[YY];
  s1ms2[ZZ] = s1[ZZ] - s2[ZZ];

  getD(fil1, D1);
  getD(fil2, D2);

  // square of the length of fil1 (square of the modulus of the vector joining the endpoints)
  D1D1 = vdotp(D1,D1);
  // cosine of the angle between fil1 and fil2 (times both modulus)
  D1D2 = vdotp(D1,D2);
  D2D2 = vdotp(D2,D2);
  D1s1s2 = vdotp(D1,s1ms2);
  D2s1s2 = vdotp(D2,s1ms2);

  // Compiled with gcc under Fedora 16, this code triggered an error (t = 1
  // in dist_betw_pt_and_fil() )
  if(fabs(D1D1) < EPS) {
    printf("Internal error in dist_betw_fils: D1D1 too small!\n");
  }
  if(fabs(D2D2) < EPS) {
    printf("Internal error in dist_betw_fils: D2D2 too small!\n");
  }

  tmp1 = D1D2*D1D2/D1D1 - D2D2;

  if(fabs(D2D2) < EPS) {
    printf("Internal error in dist_betw_fils: D2D2 too small!\n");
  }

  if (fabs(tmp1/D2D2) < EPS) {
    /* fils are parallel */
    *parallel = 1;
    return min_endpt_sep(fil1,fil2);
  }
  else
    *parallel = 0;

  t2 = (D1D2*D1s1s2/D1D1 - D2s1s2)/tmp1;
  t1 = (t2*D1D2 - D1s1s2)/D1D1;

  // debug
  //printf("t1 %g, t2 %d\n");

  if (t1 <= 1 && t1 >= 0) {
    if (t2 <= 1 && t2 >= 0) {
      getr(&x1,&y1,&z1,s1,t1,D1);
      getr(&x2,&y2,&z2,s2,t2,D2);
      return dist_between(x1,y1,z1,x2,y2,z2);
    }
    else
// debug
  if(t2 == 1.0) {
    printf("Internal error in dist_betw_fils: t2 is %g!\n", t2);
  }
//printf("t2 %d\n");

     /* nearest point along fil2 is outside line segment defining filament */
      return dist_betw_pt_and_fil(fil1, D1, s1, D1D1, fil2,t2);
  }
  else {
    if (t2 <= 1 && t2 >= 0) {

 // debug
  if(t1 == 1.0) {
    printf("Internal error in dist_betw_fils: t1 is %g!\n", t1);
  }
//printf("t1 %g\n");
     /* nearest point along fil1 is outside line segment defining filament */
      return dist_betw_pt_and_fil(fil2, D2, s2, D2D2, fil1,t1);
    }
    else
      /* both point are out of range, just compare endpoints */
      return min_endpt_sep(fil1,fil2);
  }

}

getD(fil, D)
FILAMENT *fil;
double *D;
{
  D[XX] = fil->x[1] - fil->x[0];
  D[YY] = fil->y[1] - fil->y[0];
  D[ZZ] = fil->z[1] - fil->z[0];
}

getr(x,y,z,s,t,D)
double *x,*y,*z;
double *s,t,*D;
{
  *x = s[XX] + t*D[XX];
  *y = s[YY] + t*D[YY];
  *z = s[ZZ] + t*D[ZZ];
}

double vdotp(v1, v2)
double *v1,*v2;
{
  return v1[XX]*v2[XX] + v1[YY]*v2[YY] + v1[ZZ]*v2[ZZ];
}

double dist_between(x1,y1,z1,x2,y2,z2)
double x1,y1,z1,x2,y2,z2;
{
  return sqrt(SQUARE(x1 - x2) + SQUARE(y1 - y2) + SQUARE(z1 - z2));
}

/* returns the minimum distance between an endpt of fil1 and one of fil2 */
double min_endpt_sep(fil1,fil2)
FILAMENT *fil1, *fil2;
{
  double min, tmp;
  int idx1, idx2;

  idx1 = 0;
  idx2 = 0;
  min = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);

  idx1 = 1;
  idx2 = 0;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  idx1 = 0;
  idx2 = 1;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  idx1 = 1;
  idx2 = 1;
  tmp = dist_between(fil1->x[idx1],fil1->y[idx1],fil1->z[idx1],
		     fil2->x[idx2],fil2->y[idx2],fil2->z[idx2]);
  if (tmp < min) min = tmp;

  return min;
}

/* this finds the shortest distance between the line defined by
   r = s + tnew*D, t in [0,1] (which is along fil_line) and the endpoint
   on fil which is closer to fil_line (determined by the value of t) */
double dist_betw_pt_and_fil(fil_line, D, s, DD, fil,t)
FILAMENT *fil_line, *fil;
double *D, *s, t, DD;
{
  double e[3], sme[3], x,y,z;
  double tnew, Dsme;
  int idx;

  if (t < 0) {
    e[XX] = fil->x[0];
    e[YY] = fil->y[0];
    e[ZZ] = fil->z[0];
  }
  else if (t > 1) {
    e[XX] = fil->x[1];
    e[YY] = fil->y[1];
    e[ZZ] = fil->z[1];
  }
  else {
    fprintf(stderr, "Internal err: dist_bet_pt_and_fil: why is t = %lg?\n", t);
    exit(1);
  }

  sme[XX] = s[XX] - e[XX];
  sme[YY] = s[YY] - e[YY];
  sme[ZZ] = s[ZZ] - e[ZZ];

  Dsme = vdotp(D,sme);

  tnew = -Dsme/DD;

  if (tnew <= 1 && tnew >= 0) {
    /* This will be the case when a small fil is near a big one (fil_line).*/
    /* Calculate r = (s - e) + tnew*D */
    getr(&x,&y,&z,sme,tnew,D);
    return sqrt(x*x + y*y + z*z);
  }
  else {
    /* just find the distance between the nearest endpt and e[] */
    if (tnew < 0)
      idx = 0;
    else
      idx = 1;

    return dist_between(e[XX], e[YY], e[ZZ],
			fil_line->x[idx], fil_line->y[idx], fil_line->z[idx]);
  }
}

double aspectratio(fil)
FILAMENT *fil;
{
  double rat;

  rat = fil->height/fil->width;
  if (rat >= 1.0)
    return rat;
  else
    return 1.0/rat;
}

fill_Gquad()
{
  int i,j;

  Gweight = (double **)MattAlloc(MAXsubfils+1, sizeof(double *));
  Gpoint = (double **)MattAlloc(MAXsubfils+1, sizeof(double *));

  for(i = 1; i <= MAXsubfils; i++) {
    Gweight[i] = (double *)MattAlloc(i, sizeof(double));
    Gpoint[i] = (double *)MattAlloc(i, sizeof(double));
  }

  Gweight[1][0] = 2.0;
  Gpoint[1][0] = 0.0;

  for(i = 2; i <= MAXsubfils; i++)
    /* subtract 1 from pointers since this function starts with p[1] */
    gquad_weights(i,Gpoint[i] - 1, Gweight[i] - 1);
}

findnfils(fil, subfils, nfils)
FILAMENT *fil, *subfils;
int nfils;
{
  double hx,hy,hz,mag,wx,wy,wz,dx,dy,dz;
  int i,j;

  if (fil->segm->widthdir != NULL) {
    wx = fil->segm->widthdir[XX];
    wy = fil->segm->widthdir[YY];
    wz = fil->segm->widthdir[ZZ];
  }
  else {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(fil->y[1] - fil->y[0])*1.0;
    wy = (fil->x[1] - fil->x[0])*1.0;
    wz = 0;
    if ( fabs(wx/fil->segm->length) < EPS && fabs(wy/fil->segm->length) < EPS) {
      /* if all of x-y is perpendic to length, then choose x direction */
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }

  hx = -wy*(fil->z[1] - fil->z[0]) + (fil->y[1] - fil->y[0])*wz;
  hy = -wz*(fil->x[1] - fil->x[0]) + (fil->z[1] - fil->z[0])*wx;
  hz = -wx*(fil->y[1] - fil->y[0]) + (fil->x[1] - fil->x[0])*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;

  if (fil->width > fil->height) {
    dx = fil->width*wx/2;
    dy = fil->width*wy/2;
    dz = fil->width*wz/2;
  }
  else {
    dx = fil->height*hx/2;
    dy = fil->height*hy/2;
    dz = fil->height*hz/2;
  }

  /* all mutualfil needs are the filament coordinates and length */
  for (j = 0; j < nfils; j++) {
    for(i = 0; i < 2; i++) {
      subfils[j].x[i] = fil->x[i] + dx*Gpoint[nfils][j];
      subfils[j].y[i] = fil->y[i] + dy*Gpoint[nfils][j];
      subfils[j].z[i] = fil->z[i] + dz*Gpoint[nfils][j];
    }
    subfils[j].length = fil->length;
  }
}

#if 1==0
/* main for testing if these work */
main()
{
  FILAMENT fil1, fil2;

  while(1) {
    printf("fil1 1: ");
    scanf("%lf %lf %lf", &(fil1.x[0]), &(fil1.y[0]), &(fil1.z[0]));
    printf("fil1 2: ");
    scanf("%lf %lf %lf", &(fil1.x[1]), &(fil1.y[1]), &(fil1.z[1]));

    printf("fil2 1: ");
    scanf("%lf %lf %lf", &(fil2.x[0]), &(fil2.y[0]), &(fil2.z[0]));
    printf("fil2 2: ");
    scanf("%lf %lf %lf", &(fil2.x[1]), &(fil2.y[1]), &(fil2.z[1]));

    printf("dist = %lg\n",dist_betw_fils(&fil1, &fil2));
  }
}
#endif

gquad_weights(N,p,w)
int N;
double *p,*w;
{
  int i;

  if ( !((N>=2)&&(N<=6)) &&(N!=19)&&(N!=30)&&(N!=25)&&(N!=15)&&(N!=10)) {
    fprintf(stderr,"Hey, bad number of Guassian Quad points\n");
    exit(0);
  }
  if (N==10) {
    p[1]= -.9739065285; p[10]= -p[1];
    p[2]= -.8650633666; p[9] = -p[2];
    p[3]= -.6794095682; p[8] = -p[3];
    p[4]= -.4333953941; p[7] = -p[4];
    p[5]= -.1488743389; p[6] = -p[5];
    w[1]=w[10]=.0666713443;
    w[2]=w[9]= .1494513491;
    w[3]=w[8]= .2190863625;
    w[4]=w[7]= .2692667193;
    w[5]=w[6]= .2955242247;
  }
  if (N==15) {
    p[1]= -.9879925180;
    p[2]= -.9372733924;
    p[3]= -.8482065834;
    p[4]= -.7244177313;
    p[5]= -.5709721726;
    p[6]= -.3941513470;
    p[7]= -.2011940939;
    p[8]= 0.0;
    p[9]= -p[7];
    p[10]= -p[6];
    p[11]= -p[5];
    p[12]= -p[4];
    p[13]= -p[3];
    p[14]= -p[2];
    p[15]= -p[1];
    w[1]=w[15]= .0307532419;
    w[2]=w[14]= .0703660474;
    w[3]=w[13]= .1071592204;
    w[4]=w[12]= .1395706779;
    w[5]=w[11]= .1662692058;
    w[6]=w[10]= .1861610000;
    w[7]=w[9]= .1984314853;
    w[8]= .2025782419;
  }
  if (N==25){
    p[1]= -.9955568687; p[25]= -p[1];
    p[2]= -.9766639214; p[24]= -p[2];
    p[3]= -.9429745712; p[23]= -p[3];
    p[4]= -.8949919978; p[22]= -p[4];
    p[5]= -.8334426287; p[21]= -p[5];
    p[6]= -.7592592630; p[20]= -p[6];
    p[7]= -.6735663684; p[19]= -p[7];
    p[8]= -.5776629302; p[18]= -p[8];
    p[9]= -.4730027314; p[17]= -p[9];
    p[10]= -.3611723058; p[16]= -p[10];
    p[11]= -.2438668837; p[15]= -p[11];
    p[12]= -.1228646926; p[14]= -p[12];
    p[13]= 0.0;
    w[1]=w[25]= .0113937985;
    w[2]=w[24]= .0263549866;
    w[3]=w[23]= .0409391567;
    w[4]=w[22]= .0549046959;
    w[5]=w[21]= .0680383338;
    w[6]=w[20]= .0801407003;
    w[7]=w[19]= .0910282619;
    w[8]=w[18]= .1005359490;
    w[9]=w[17]= .1085196244;
    w[10]=w[16]= .1148582591;
    w[11]=w[15]= .1194557635;
    w[12]=w[14]= .1222424429;
    w[13]= .1231760537;
  }
  if (N==2) {
    p[1]= -.57735027;
    p[2]= .57735027;
    w[1]=w[2]=1.0;
  }
  if (N==3) {
    p[1]= -.77459667;
    p[2]= 0.;
    p[3]= .77459667;
    w[1]= .55555555;
    w[2]= .88888888;
    w[3]= .55555555;
  }
  if (N==4) {
    p[1]= -.86113631;
    p[2]= -.33998104;
    p[3]= .33998104;
    p[4]= .86113631;
    w[1]=w[4]=.34785485;
    w[2]=w[3]=.65214515;
  }
  if (N==5){
    p[1] = -.90617984;
    p[2] = -.53846931;
    p[3] = 0.;
    p[4] = .53846931;
    p[5] = .90617984;
    w[1]=w[5]=.23692688;
    w[2]=w[4]=.47862867;
    w[3]=.56888888;
  }
  if (N==6){
    p[1]= -.93246951;
    p[2]= -.66120938;
    p[3]= -.23861918;
    p[4]= -p[3];
    p[5]= -p[2];
    p[6]= -p[1];
    w[1]=w[6]= .17132449;
    w[2]=w[5]= .36076157;
    w[3]=w[4]= .46791393;
  }
  if (N==19) {
    p[1]= -.9931285991;
    p[2]= -.9602081521;
    p[3]= -.9031559036;
    p[4]= -.8227146565;
    p[5]= -.7209661773;
    p[6]= -.6005453046;
    p[7]= -.4645707413;
    p[8]= -.3165640999;
    p[9]= -.1603586456;
    p[10] = 0.;
    p[11] = -p[9];
    p[12] = -p[8];
    p[13] = -p[7];
    p[14] = -p[6];
    p[15] = -p[5];
    p[16] = -p[4];
    p[17] = -p[3];
    p[18] = -p[2];
    p[19] = -p[1];
    w[1]=w[19]= .0194617882;
    w[2]=w[18]= .0448142267;
    w[3]=w[17]= .0690445427;
    w[4]=w[16]=.0914900216;
    w[5]=w[15]=.1115666455;
    w[6]=w[14]=.1287539625;
    w[7]=w[13]=.1426067021;
    w[8]=w[12]=.1527660420;
    w[9]=w[11]=.1589688433;
    w[10]=.1610544498;
  }
  if (N==30) {
    p[1]= -.9968934840;
    p[2]= -.9836681232;
    p[3]= -.9600218649;
    p[4]= -.9262000474;
    p[5]= -.8825605357;
    p[6]= -.8295657623;
    p[7]= -.7677774321;
    p[8]= -.6978504947;
    p[9]= -.6206261829;
    p[10]= -.5366241481;
    p[11]= -.4470337695;
    p[12]= -.3527047255;
    p[13]= -.2546369261;
    p[14]= -.1538699136;
    p[15]= -.0514718425;
    p[16]= -p[15];
    p[17]= -p[14];
    p[18]= -p[13];
    p[19]= -p[12];
    p[20]= -p[11];
    p[21]= -p[10];
    p[22]= -p[9];
    p[23]= -p[8];
    p[24]= -p[7];
    p[25]= -p[6];
    p[26]= -p[5];
    p[27]= -p[4];
    p[28]= -p[3];
    p[29]= -p[2];
    p[30]= -p[1];
    w[1]=w[30]= .0079681924;
    w[2]=w[29]= .0184664683;
    w[3]=w[28]= .0287847078;
    w[4]=w[27]= .0387991925;
    w[5]=w[26]= .0484026728;
    w[6]=w[25]= .0574931562;
    w[7]=w[24]= .0659742298;
    w[8]=w[23]= .0737559747;
    w[9]=w[22]= .0807558952;
    w[10]=w[21]= .0868997872;
    w[11]=w[20]= .0921225222;
    w[12]=w[19]= .0963687371;
    w[13]=w[18]= .0995934205;
    w[14]=w[17]= .1017623897;
    w[15]=w[16]= .1028526528;
  }
}





