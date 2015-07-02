/* # ***** sort to /src/header
   # ***** */
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
