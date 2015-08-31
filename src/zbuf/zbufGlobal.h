/*
Copyright (c) 1990 Massachusetts Institute of Technology, Cambridge, MA.
All rights reserved.

This Agreement gives you, the LICENSEE, certain rights and obligations.
By using the software, you indicate that you have read, understood, and
will comply with the terms.

Permission to use, copy and modify for internal, noncommercial purposes
is hereby granted.  Any distribution of this program or any part thereof
is strictly prohibited without prior written consent of M.I.T.

Title to copyright to this software and to any associated documentation
shall at all times remain with M.I.T. and LICENSEE agrees to preserve
same.  LICENSEE agrees not to make any copies except for LICENSEE'S
internal noncommercial use, or to use separately any portion of this
software without prior written consent of M.I.T.  LICENSEE agrees to
place the appropriate copyright notice on any such copies.

Nothing in this Agreement shall be construed as conferring rights to use
in advertising, publicity or otherwise any trademark or the name of
"Massachusetts Institute of Technology" or "M.I.T."

M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By
way of example, but not limitation, M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS OR DOCUMENTATION WILL
NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
M.I.T. shall not be held liable for any liability nor for any direct,
indirect or consequential damages with respect to any claim by LICENSEE
or any third party on account of or arising from this Agreement or use
of this software.
*/


/* #include <stdio.h> */

/* zbuf data structures */
#include "zbufStruct.h"

#ifndef MIN
#define MIN(A,B)  ( (A) > (B) ? (B) : (A) )
#endif
#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef M_PI
#define M_PI 3.1415926535897931160E0
#endif

#ifndef ON
#define ON 1
#endif
#ifndef OFF
#define OFF 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define FONT 5.0                /* font size used to label lines in ps file */
#define CMDFONT 10.0		/* font used to write command line */

#define OFFSETX 34.0		/* offset on x from lower left (pnts) */
#define OFFSETY 34.0		/* offset on y from lower left (pnts) */
#define IMAGEX 540.0		/* x height of image (chosen to be points) */
#define IMAGEY 720.0		/* y height of image (chosen to be points)
				   comand line scale parameter muliplies
				   IMAGEX/Y to get final image size
				   - these values should allow 7.5x10 images */
#define KEYHGT IMAGEY/7.0	/* height of the shading key, -q option only */
#define KEYWID IMAGEX/7.5	/* width of the shading key, -q option only */
#define KEYBLKS 5		/* number of blocks in the key */
#define KEYPREC 5		/* precision of labels in key */
#define KEYFONT 10.0		/* font used to label key */

#define MARGIN 1e-10		/* 1 + MARGIN is considered 1; > -MARGIN=0;
				   MARGIN should be approx machine epsilon
				   except that PATRAN files can lead to probs*/
#define PARMGN 1.0		/* multiplies MARGIN in plane equivalence
				   checks (see is1stFaceDeeper()) to allow the
				   extra slop those checks apparently need */

#define LINE 1			/* used directly in setlinewidth as default */
#define LINCAP 1       		/* 0 = square, 1 = round cap, 2 = square cap */
#define LINJIN 1		/* 0 = miter, 1 = round, 2 = knock off cnrs */
#define GREYLEV 0.0		/* 1 = black fill */
#define DASHED -1		/* = width arg that means draw a dashed line */
#define DASWTH 2		/* width of dashed lines */
#define OVRWTH 0		/* overide line width - get with -w option
				   (not implemented) */
#define AXEWID 0.0		/* width of axis lines */
#define MAXSIDES 4		/* maximum #sides allowed for a face */

#define DEBUG 1			/* controls view point dump, etc. */
#define DEBUGX OFF		/* doLinesIntersect stats */
#define DMPINFO OFF		/* dumps face/line info when in a bind */
#define DMPMATLAB OFF		/* dmps info in MATLAB format if DMPINFO==ON */
#define RMWEDGE OFF		/* removes 1st quadrant dielectric panels */

#define ALEN 8			/* default arrow length in points */
#define AWID 4			/* default arrow width in points */
#define DOTSIZ 2		/* default dot radius, points */

#define POS 0			/* defined in whichSide() */
#define NEG 1
#define SPLIT 2
#define SAME 3

#define REVERSE 2		/* defined in is1stFaceDeeper() */

#define XOVTST OFF		/* cross overlap test enable--see zbufSort.c */

#ifndef XI
#define XI 0			/* upward-pointing axes types */
#endif
#ifndef YI
#define YI 1
#endif
#ifndef ZI
#define ZI 2
#endif

/* default command line options */
#define DEFAZM 50.0      	/* default azimuth, degrees (-a) */
#define DEFELE 50.0		/* default elevation, degrees (-e) */
#define DEFROT 0.0		/* default rotation rel z axis, degrees (-r) */
#define DEFDST 2.0		/* default view dist, 0 = on obj edge (-d) */
#define DEFSCL 1.0		/* default scale, fractions of IMAGEX,Y (-s) */
#define DEFWID 1.0		/* default line width, points (-w) */
#define DEFAXE 1.0		/* default axes length (-x) */
#define DEFUAX ZI		/* default upward-pointing axis (-u) */

/* extras.c */
void dump_face(FILE*, face*);

/* zbuf.c */
char *concat3(char*, char*, char*);
double *get_q(char*, charge*);

/* zbuf2fastcap.c */
void dump_ps_geometry(charge*, double*, int, int);

/* zbufInOut.c */
void setupLine(double***, int, double, double, double, double, double,
    double);
void figure_grey_levels(face**, double*, charge*, int);
void get_charge_densities(double*, char*, int);
void getAbsCoord(double*, charge*, int);
face **fastcap2faces(int*, charge*, double*, int);
void readLines(FILE*, line**, line**, int*);
line **getLines(char*, int*);
void getBndingBox(face**, int, line**, int, int*, int*, FILE*, double***);
void dumpAxes(double***, FILE*);
void copyBody(FILE*);
void numberFaces(face**, int, FILE*);
void numberFace(face*, FILE*);
void dumpAdjGraph(face**, int, FILE*);
void dumpFaceText(face**, int, FILE*);
void dump_line_as_ps(FILE*, char*, double, double, double);
void dump_shading_key(FILE*, int, int, double, int);
void numberLines(line**, int, FILE*);
void dumpLines(FILE*, line**, int);
void dumpPs(face**, int, line**, int, FILE*, char**, int, int);

/* zbufProj.c */
void image(face**, int, line**, int, double*, double, double*);
void initFaces(face**, int, double*);
void flatten(face**, int, line**, int, double, double, double*, double*);
void makePos(face**, int, line**, int);
void scale2d(face**, int, line**, int, double, double*);
double *getAvg(face**, int, line**, int, int);
double getSphere(double*, face**, int, line**, int);
double getNormal(double*, double, double*, double*, double);

/* zbufSort.c */
int diff_is_zero(double, double, double);
int diff_is_negative(double, double, double);
double dot(double*, double*);
void crossProd(double*, double*, double*);
double getPlane(double*, double*, double*, double*);
int whichSide(face*, face*);
int doLinesIntersect(double*, double*, double*, double*, double*);
int face_is_inside(double**, int, double**, int, double*);
int is1stFaceDeeper(face*, face*, double*, double, double*);
int isThereBoxOverlap(face*, face*, double*);
int chkCycle(face*, face*, FILE*);
void dumpCycles(face**, int, FILE*);
void setDepth(face*);
face **depthSortFaces(face**, int);
void getAdjGraph(face**, int, double*, double, double*);

