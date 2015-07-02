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


/* structs used by .ps file dumping part of code (zbuf)--somewhat redundant */

struct face {
  int numsides;			/* number of sides this face has */
  double **c;			/* corners of the face */
  double normal[3];		/* normal to the face's plane */
  double rhs;			/* rhs for the face's plane equation */
  int index;			/* input order index */
  int depth;			/* depth index - lower numbers are deeper */
  int mark;			/* flag for topological depth ordering */
  double greylev;		/* 0 (white) to 1 (black), default = GREYLEV */
  double width;			/* line width, default = LINE */
  int numbehind;		/* number of faces this face is behind */
  struct face **behind;		/* pntrs to faces this face is behind */
  struct face *prev;
  struct face *next;
};
typedef struct face face;

struct line {
  double from[3];
  double to[3];
  int index;
  int width;
  double arrow;			/* != 0.0 => put arrow hd on to end this sz */
  double dot;			/* != 0.0 => put dot on to end this sz */
  struct line *prev;
  struct line *next;
};
typedef struct line line;
