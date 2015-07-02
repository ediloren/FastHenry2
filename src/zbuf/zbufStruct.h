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
