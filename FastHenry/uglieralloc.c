// 
// A new memory allocator added by Enrico Di Lorenzo, 2001/06/22
// because neither FastHenry nor the Sparse library free all the
// memory they dinamically reserve; however, the ualloc function
// only frees globally.
// The new allocator keeps trace of pointers and is 
// capable to free them one by one or all at one time
//
// The new functions are:
// dalloc, to be used instead of malloc or calloc (zeroes memory in any case)
// drealloc, instead of realloc
// dfree, instead of free
// dfreeall, to free all the (remaining) allocated memory
//
// in mulGlobal.h are defined the macros which use the above functions
// throughout the code; the macros are:
// DCALCORE, DMALCORE, DREALCORE, DFREECORE
//

/* # ***** sort to /src/misc
   # ***** */

/* A memory allocator with only the ability to free all memory allocated.
   It is a modification of the description below.    MK 11/95
   */

/* 
 memory allocator for fastcap
 - almost identical to Kerigan & Ritchie sec 8.7
 - differs in that there is no free list since fastcap never frees any memory
 - also the amount of memory sbrk()'d is exactly equal (within a unit) to 
   the requested amount if its more than NALLOC units
   this cuts down on unused blocks usually added to the free list
   while still keeping the number of sbrk()'s down 
 - the regular allocation functions are still availible for stdio to use
   although  stdio buffers are passed over when new memory is set up
 - machine dependancy should be confined to the type `ALIGN' and
   `MORECORE()' (i.e. non-UNIX machines won't have brk(), sbrk())
 - no attempt is made to make allocation efficient in terms of virtual pages
*/

// #defines avoid the compiler warnings about unsafe standard functions.
// 
// Remark: MUST be at the beginning of the file, before any stdc or crt include
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "FHWindow.h" // Enrico

#define NALLOC 8184		/* >= sizeof(HEADER)*NALLOC bytes sbrk()'d */
#define MAGICN 0xaaaaaaaaL	/* used to check fidelity of allocated blks */
#define UGDEBG 0		/* 0=> most efficient, no verification
				   1=> ualloc_verify() checks magicno, size
				   2=> ualloc_verify() prints lots to stdout
				   - can reduce efficiency due to big HEADER */

typedef double ALIGN;		/* all allocated mem starts on an ALIGN bdry */

union header {			/* allocated block header - all mem allocated
				   in units of sizeof(HEADER) bytes */
  struct {
#if UGDEBG == 1 || UGDEBG == 2
    long long magicno;		/* to check for overwrites */
    union header *ptr;		/* next allocated block */
#endif
#if UGDEBG == 2
    union header *enlarged;	/* limit of newly enlarged block */
    union header *enlgdfrom;	/* base pointer of newly added block */
    unsigned long long osize;		/* size before enlargement/truncation */
    unsigned long long request;      	/* nbytes requested from allocator */
#endif
    unsigned long long size;		/* block size (memory + 1 or 2 header units) */
  } s;
  union header *lastblock;      /* For freeing (MK  11/95) */
  ALIGN x;			/* for alignment only, never referenced */
};

typedef union header HEADER;

#define ALLOC_BUFFER_LEN 1024

	char cAllocBegin = 0;
	unsigned int uiBufferLevel;
	struct buffer {
		char *pBuffer[ALLOC_BUFFER_LEN];
		struct buffer *pNextBuffer;
	};
	struct buffer *pFirstBuf;
	struct buffer *pCurrBuf;


/* 
  uses calloc to get a new block of memory
  - any header space put in by calloc is wasted
  - zeros memory
  - an alternative to mocore() but should only be used if sbrk() doesnt zero
*/
#define MORECORE(SIZE) (HEADER *)calloc(1, SIZE*sizeof(HEADER))
//char *calloc();
//char *malloc();

static HEADER *base = NULL;    	/* base of allocated block list */
static HEADER *allocp = NULL;	/* last allocated block */
static unsigned int sizeofHDR = sizeof(HEADER);

static HEADER *lastblock = NULL;       /* pointer to last allocated block */
                                /* for freeing purposes */

/*
  asks operating system for more memory which is added to the top block
  - memory not zeroed out
*/
static HEADER *mocore(nu)
unsigned int nu;
{
  HEADER *cp;

  cp = (HEADER *)calloc(nu+1,sizeofHDR);
  if(cp == NULL) return(NULL);

  /* set pointer to last block for freeing */
  cp->lastblock = lastblock;
  lastblock = cp;

#ifdef MATTDEBUG2
  viewprintf(stderr,"%x\n", cp);
#endif

  return(cp + 1);
}

/* fills space with value */
fill_it(mem, k, nbytes)
void *mem;
char k;
int nbytes;
{
  /* printf("will fill: %x + 1 = %x\n",(HEADER *)mem - 1, (HEADER *)mem);*/
  memset( mem, k, nbytes);
}

/* frees lastblock and  the linked list of each block before it
   */
void ufree()
{
  HEADER *ptr, *next;

  ptr = lastblock;
  while(ptr != NULL) {
    next = ptr->lastblock;

#ifdef MATTDEBUG2
    viewprintf(stderr,"%x\n", ptr);
#endif

    free(ptr);
    ptr = next;
  }

  allocp = lastblock = NULL;

}
    
     
/* 
  ugly storage allocator 
  - no frees
  - allocates in blocks at least NALLOC*sizeof(HEADER) bytes long
  - ultimately uses mocore(), since no frees are done (sbrk() zeros added
    memory) this allocator performs like calloc() w/no explicit assigns to 0
*/

char *ualloc(nbytes)
unsigned long long nbytes;
{
  HEADER *mocore();
  HEADER *q;
  long long nunits;			// size in number of sizeof(HEADER)'s 
  long long brkunits;			// number to add to heap 
#if UGDEBG == 2
  HEADER *dummy = 0x0;
  HEADER *dummyret = 0x0;
#endif
  static unsigned long long topblksize;  	// size of current top block 

  if(nbytes == 0) return(NULL);	// should probably return something else 

#if UGDEBG == 1 || UGDEBG == 2
  nunits = 3+(nbytes-1)/sizeofHDR; // rm for 2 hdrs too 
#else
  nunits = 2+(nbytes-1)/sizeofHDR; // rm for 1 hdr (gets subtracted later) 
#endif

  if((q = allocp) == NULL) {	// no allocation yet 
    if((q = allocp = base = mocore(1)) == NULL) return(NULL);
#ifdef MATTDEBUG
    viewprintf(stderr,"started allocating:  size of header == %d\n",sizeof(HEADER));
#endif
    topblksize = 1;
#if UGDEBG == 1 || UGDEBG == 2
    base->s.size = 1;
    base->s.ptr = base;
    base->s.magicno = MAGICN;
#endif
  }

  
  //  check previously allocated block for room
  //  - if it's big enough, use head (not tail) part of it
  //  - if it's not, break out more memory and discard the previous block
  
  if (topblksize >= (unsigned int) nunits) {	// if it's big enough 
#if UGDEBG == 1 || UGDEBG == 2
	 HEADER *p;
    p = q;			// copy old top block pointer 
    q += (nunits - 1);	// make q new top block pointer 
    q->s.size = p->s.size - nunits + 1; // upper block size 
#else
    allocp = q + nunits - 1;
#endif
    topblksize -= (nunits - 1);
#if UGDEBG == 2
    p->s.osize = p->s.size;	// save size before new block carved out 
    p->s.request = nbytes;	// actual #bytes requested 
    p->s.enlarged = dummy;
    p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
    q->s.ptr = base;	// top block connected back to base 
    p->s.ptr = q;		// old top block pointer to new top block 
    if(allocp != p) allocp->s.ptr = p;	// link previous old block 
    q->s.magicno = MAGICN;
    p->s.size = nunits;	// set up returned block size 
    allocp = q;		// save the top block pointer 
    return((char *) (p+1));	// return pntr to first unit beyond header 
#else
    return((char *) q);	// use old header location 
#endif
  }
  else {			// get more memory, add to top block beyond 
				    // any stdio buffers made since last ualloc 
    brkunits = (nunits > NALLOC) ? nunits : NALLOC;
    if((q = mocore(brkunits)) == NULL) return(NULL);
#if UGDEBG == 1 || UGDEBG == 2
    q->s.size = brkunits;	// setup size, etc. of new top block 
    q->s.magicno = MAGICN;
#endif
    topblksize = brkunits;
#if UGDEBG == 2
    dummy = q + q->s.size - 1; // limit pointer of new block 
    dummyret = q;		// base pointer of new block 
#endif
    // repeat above code to aviod having a loop 
#if UGDEBG == 1 || UGDEBG == 2
    p = q;
    q += (nunits - 1);	// make q new top block pointer 
    q->s.size = p->s.size - nunits + 1; // upper block size 
#else
    allocp = q + nunits - 1;
#endif
    topblksize -= (nunits - 1);
#if UGDEBG == 2
    p->s.osize = p->s.size;	// save size before new block carved out 
    p->s.request = nbytes;	// actual #bytes requested 
    p->s.enlarged = dummy;
    p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
    q->s.ptr = base;	// top block connected back to base 
    p->s.ptr = q;		// old top block pointer to new top block 
    if(allocp != p) allocp->s.ptr = p;	// link previous old block 
    q->s.magicno = MAGICN;
    p->s.size = nunits;	// set up returned block size 
    allocp = q;		// save the top block pointer 
    return((char *) (p+1));	// return pntr to first unit beyond header 
#else
    return((char *) q);	// use old header location 
#endif
  }
} 

/*
  checks for overwrites at end of blocks
  - checks magic numbers and sizes (UGDEBG == 1 or 2)
  - checks if length corresponds to pointers (UGDEBG == 1 or 2)
  - prints information about list of allocated blocks (UGDEBG == 2)
*/
void ualloc_verify()
{
  int cnt = 1;

#if UGDEBG == 1 || UGDEBG == 2
  HEADER *p;
  for(p = base;; p = p->s.ptr, cnt++) {
    if(p == base && cnt > 1) break;
#endif
#if UGDEBG == 2
    viewprintf(stdout, 
	   "%d 0x%x 0x%x %u %u bytes (osize %u enlarged from 0x%x to 0x%x)\n", 
	    cnt, p, p->s.ptr, p->s.size, p->s.request,
	    p->s.osize, p->s.enlgdfrom, p->s.enlarged);
#endif
#if UGDEBG == 1 || UGDEBG == 2
    if(p->s.magicno != MAGICN) {
      fflush(stdout);
      viewprintf(stderr, "ualloc_verify: bad block %d magic number\n", cnt);
      fflush(stderr);
      FHExit(FH_NORMAL_END);
    }
    if(p->s.ptr - p + 1 < p->s.size) {
      fflush(stdout);
      viewprintf(stderr, "ualloc_verify: bad block size ", cnt);
      viewprintf(stderr, "(from 0x%x to 0x%x size %u)\n", p, p->s.ptr, p->s.size);
      fflush(stderr);
      FHExit(FH_NORMAL_END);
    }
  }
#endif
}

/*
  compares the total address range broken out to amount requested
  - efficiency is 100*memcount/(sbrk(0)-base) = 100*requested/(ttl broken out)
  - UGDEBG == 2 prints waste values: 
      ualloc waste = 100*(mem discarded)/(ttl broken out)
      stdio waste = 100*(mem discarded to apease stdio usage)/(ttl broken out)
       (stdio waste is often zero and ualloc waste does not include the
        memory lost in each header struct (a problem with many small things)
  - if base == NULL (not using ugly allocator), final break value is printed
*/
void uallocEfcy(memcount)
long long memcount;
{
#if UGDEBG == 2
  HEADER *p;
  unsigned int waste = 0;
  unsigned int stdiowaste = 0;
  int first = 1;
#endif
  long long total;
  char *sbrk();

  total = (int)(sbrk(0) - (char *)base);

  if(base == NULL) viewprintf(stdout, "(top of memory = 0x%x", sbrk(0));
  else viewprintf(stdout, "(%.3g%% efficiency",
	       100*((double)memcount)/((double)total));

#if UGDEBG == 2
  /* loop through blocks, get total btyes wasted (by allocator) */
  for(p = base;; p = p->s.ptr) {
    if(p == base && first == 0) break;
    if(p == NULL) break;	/* compatability when ualloc not used */
    if(p->s.request == 0) {
      if((p->s.ptr)->s.size >= p->s.size) { /* if block dropped bec. too sm. */
	if(p != allocp) waste += p->s.size - 1; 
	else waste += p->s.size;
      }
      else {
	if(p != allocp) stdiowaste += p->s.size - 1; 
      }	
    }
    first = 0;
  }
  if(p != NULL) 
      viewprintf(stdout, 
	      ", waste: %.2g%% (ualloc), %.2g%% (stdio)", 
	      100*(double)(waste*sizeof(HEADER))/((double)total),
	      100*(double)(stdiowaste*sizeof(HEADER))/((double)total));
#endif
  viewprintf(stdout, ")\n");

}

/* 
  debug storage allocator - Enrico
  is used to alloc memory instead of calloc and malloc, to keep trace
  of allocated memory

  uses map table to store pointers to allocated memory, so it will be
  able to globally free it (and to signal which blocks aren't freed) 
*/
char *dalloc(unsigned long long nbytes)
{
	struct buffer *tmpBuf;
	unsigned long long numUnits;
	char *pointer;

	if (cAllocBegin == 0) {
		pCurrBuf = pFirstBuf = (struct buffer *) malloc( sizeof(struct buffer) );
		uiBufferLevel = 0;
		cAllocBegin = 1;
	}

	if(uiBufferLevel > ALLOC_BUFFER_LEN - 1) {
		tmpBuf = pCurrBuf;
		pCurrBuf = (struct buffer *) malloc( sizeof(struct buffer) );
		tmpBuf->pNextBuffer = pCurrBuf;
		uiBufferLevel = 0;
	}

	numUnits = nbytes / sizeof(char) + nbytes % sizeof(char);
	pointer = (char*) pCurrBuf->pBuffer[uiBufferLevel++] = calloc(numUnits, sizeof(char));
	if(pointer == NULL)  {
		viewprintf(stderr, "Out of memory!\n");
		FHExit(FH_GENERIC_ERROR);
	}
	return(pointer);
}

void *drealloc(void *pointer, unsigned long long size) 
{
	struct buffer *tmpBuf;
	long long i, last;

	tmpBuf = pFirstBuf;
	do {
		if (tmpBuf == pCurrBuf) {
			last = uiBufferLevel;
		}
		else {
			last = ALLOC_BUFFER_LEN;
		}
		
		for(i=0; i<last; i++) {
			if (tmpBuf->pBuffer[i] == pointer) {
				tmpBuf->pBuffer[i] = realloc(pointer, size);
				return tmpBuf->pBuffer[i];
			}
		}

		if (tmpBuf != pCurrBuf) {
			tmpBuf = tmpBuf->pNextBuffer;
		}
		else 
			break;

	} while (1);

	viewprintf(stderr, "WARNING, drealloc: cannot find pointer to be reallocated");

	return NULL;
}

void dfree(void *pointer)
{
	struct buffer *tmpBuf;
	long long i, last;

	tmpBuf = pFirstBuf;
	do {
		if (tmpBuf == pCurrBuf) {
			last = uiBufferLevel;
		}
		else {
			last = ALLOC_BUFFER_LEN;
		}
		
		for(i=0; i<last; i++) {
			if (tmpBuf->pBuffer[i] == pointer) {
				free(pointer);
				tmpBuf->pBuffer[i] = NULL;
				return;
			}
		}

		if (tmpBuf != pCurrBuf) {
			tmpBuf = tmpBuf->pNextBuffer;
		}
		else 
			break;

	} while (1);

	viewprintf(stderr, "WARNING, dfree: cannot find pointer to be freed");
}

void dfreeall()
{
	struct buffer *tmpBuf, *tmpBuf2;
	long long i, last;

	tmpBuf = pFirstBuf;
	do {
		if (tmpBuf == pCurrBuf) {
			last = uiBufferLevel;
		}
		else {
			last = ALLOC_BUFFER_LEN;
		}
		
		for(i=0; i<last; i++) {
			if (tmpBuf->pBuffer[i] != NULL) {
				free(tmpBuf->pBuffer[i]);
			}
		}

		if (tmpBuf != pCurrBuf) {
			tmpBuf2 = tmpBuf;
			tmpBuf = tmpBuf->pNextBuffer;
			free(tmpBuf2);
		} 
		else {
			free(tmpBuf);
			break;
		}
	} while (1);
}

char stringd[10] = "0"; /* Enrico */
char *sbrk()        /* Enrico */
{return(stringd);} /* Enrico */


void InitUglierAllocVars(void)
{ 
	// global
	cAllocBegin = 0;
	strcpy(stringd, "0"); 

	// static, global or made global to be initializable
	base = NULL;
	allocp = NULL;
	sizeofHDR = sizeof(HEADER);
	lastblock = NULL;
}



