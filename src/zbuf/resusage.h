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

/* header where rusage and time structs are defined */

#ifdef FOUR
#define NOTOTHER 1
#include <sys/time.h>
#include <sys/resource.h>
struct rusage timestuff;
#endif

#ifdef FIVE
#define NOTOTHER 1
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
struct tms timestuff;
#endif

/* define macros for time and resident memory usage checks */

static double dtime = 0.0;
static long sectime, utime;

#ifdef NOTOTHER

#ifdef FOUR			/* 4.2,3BSD (tested: Sun4, IBM6000, DEC5000) */
#define starttimer getrusage(RUSAGE_SELF, &timestuff); \
sectime = timestuff.ru_utime.tv_sec; \
utime = timestuff.ru_utime.tv_usec
#define stoptimer getrusage(RUSAGE_SELF, &timestuff); \
dtime = (double)(timestuff.ru_utime.tv_sec - sectime) \
        + 1.0e-6*(double)(timestuff.ru_utime.tv_usec - utime)
#define DUMPRSS			/*  */
#endif /* FOUR */

#ifdef FIVE			/* for System V (tested: HP300) */
#define starttimer times(&timestuff); \
utime = timestuff.tms_utime
#define stoptimer times(&timestuff); \
dtime = (timestuff.tms_utime)-utime; \
dtime /= HZ
#define DUMPRSS			/*  */
#endif /* FIVE */

#else				/* default - no timers */

#define starttimer		/*  */
#define stoptimer		/*  */
#define DUMPRSS			/*  */

#endif /* NOTOTHER */
