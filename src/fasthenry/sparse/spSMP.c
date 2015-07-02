/*
 *  Spice3 COMPATIBILITY MODULE
 *
 *  Author:                     Advising professor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This module contains routines that make Sparse1.3 a direct
 *  replacement for the SMP sparse matrix package in Spice3b1.
 *  Sparse1.3 is in general a faster and more robust package than SMP.
 *  These advantages become significant on large circuits.
 *
 *  >>> User accessible functions contained in this file:
 *  SMPaddElt
 *  SMPmakeElt
 *  SMPcClear
 *  SMPclear
 *  SMPcLUfac
 *  SMPluFac
 *  SMPcReorder
 *  SMPreorder
 *  SMPcSolve
 *  SMPsolve
 *  SMPmatSize
 *  SMPnewMatrix
 *  SMPdestroy
 *  SMPpreOrder
 *  SMPprint
 *  SMPgetError
 *  SMPcProdDiag
 */

/*
 *  To replace SMP with Sparse, rename the file spSpice3.h to
 *  spMatrix.h and place Sparse in a subdirectory of SPICE called
 *  `sparse'.  Then on UNIX compile Sparse by executing `make spice'.
 *  If not on UNIX, after compiling Sparse and creating the sparse.a
 *  archive, compile this file (spSMP.c) and spSMP.o to the archive,
 *  then copy sparse.a into the SPICE main directory and rename it
 *  SMP.a.  Finally link SPICE.
 *
 *  To be compatible with SPICE, the following Sparse compiler options
 *  (in spConfig.h) should be set as shown below:
 *
 *      REAL                            YES
 *      EXPANDABLE                      YES
 *      TRANSLATE                       NO
 *      INITIALIZE                      NO or YES, YES for use with test prog.
 *      DIAGONAL_PIVOTING               YES
 *      ARRAY_OFFSET                    YES
 *      MODIFIED_MARKOWITZ              NO
 *      DELETE                          NO
 *      STRIP                           NO
 *      MODIFIED_NODAL                  YES
 *      QUAD_ELEMENT                    NO
 *      TRANSPOSE                       NO
 *      SCALING                         NO
 *      DOCUMENTATION                   YES
 *      MULTIPLICATION                  NO
 *      DETERMINANT                     YES
 *      STABILITY                       NO
 *      CONDITION                       NO
 *      PSEUDOCONDITION                 NO
 *      FORTRAN                         NO
 *      DEBUG                           YES
 *      spCOMPLEX                       1
 *      spSEPARATED_COMPLEX_VECTORS     1
 *      spCOMPATIBILITY                 0
 *
 *      spREAL  double
 */

/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and its
 *  documentation for any purpose and without fee is hereby granted, provided
 *  that the above copyright notice appear in all copies and supporting
 *  documentation and that the authors and the University of California
 *  are properly credited.  The authors and the University of California
 *  make no representations as to the suitability of this software for
 *  any purpose.  It is provided `as is', without express or implied warranty.
 */

#ifndef lint
static char copyright[] =
    "Sparse1.3: Copyright (c) 1985,86,87,88 by Kenneth S. Kundert";
static char RCSid[] =
    "@(#)$Header: spSMP.c,v 1.2 88/06/24 05:02:42 kundert Exp $";
#endif




/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spMatrix.h
 *     Sparse macros and declarations.
 *  SMPdefs.h
 *     Spice3's matrix macro definitions.
 */

#include "spMatrix.h"
#include "../include/SMPdefs.h"

#define NO   0
#define YES  1


/*
 * SMPaddElt()
 */
int
SMPaddElt( Matrix, Row, Col, Value )
SMPmatrix *Matrix;
int Row, Col;
double Value;
{
    *spGetElement( (char *)Matrix, Row, Col ) = Value;
    return spError( (char *)Matrix );
}

/*
 * SMPmakeElt()
 */
double *
SMPmakeElt( Matrix, Row, Col )
SMPmatrix *Matrix;
int Row, Col;
{
    return spGetElement( (char *)Matrix, Row, Col );
}

/*
 * SMPcClear()
 */
void
SMPcClear( Matrix )
SMPmatrix *Matrix;
{
    spClear( (char *)Matrix );
}

/*
 * SMPclear()
 */
void
SMPclear( Matrix )
SMPmatrix *Matrix;
{
    spClear( (char *)Matrix );
}

/*
 * SMPcLUfac()
 */
/*ARGSUSED*/
int
SMPcLUfac( Matrix, PivTol )
SMPmatrix *Matrix;
double PivTol;
{
    spSetComplex( (char *)Matrix );
    return spFactor( (char *)Matrix );
}

/*
 * SMPluFac()
 */
/*ARGSUSED*/
int
SMPluFac( Matrix, PivTol, Gmin )
SMPmatrix *Matrix;
double PivTol, Gmin;
{
    spSetReal( (char *)Matrix );
    LoadGmin( (char *)Matrix, Gmin );
    return spFactor( (char *)Matrix );
}

/*
 * SMPcReorder()
 */
int
SMPcReorder( Matrix, PivTol, PivRel, NumSwaps )
SMPmatrix *Matrix;
double PivTol, PivRel;
int *NumSwaps;
{
    *NumSwaps = 0;
    spSetComplex( (char *)Matrix );
    return spOrderAndFactor( (char *)Matrix, (spREAL*)NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}

/*
 * SMPreorder()
 */
int
SMPreorder( Matrix, PivTol, PivRel, Gmin )
SMPmatrix *Matrix;
double PivTol, PivRel, Gmin;
{
    spSetComplex( (char *)Matrix );
    LoadGmin( (char *)Matrix, Gmin );
    return spOrderAndFactor( (char *)Matrix, (spREAL*)NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}

/*
 * SMPcSolve()
 */
void
SMPcSolve( Matrix, RHS, iRHS, Spare, iSpare)
SMPmatrix *Matrix;
double RHS[], iRHS[], Spare[], iSpare[];
{
    spSolve( (char *)Matrix, RHS, RHS, iRHS, iRHS );
}

/*
 * SMPsolve()
 */
void
SMPsolve( Matrix, RHS, Spare )
SMPmatrix *Matrix;
double RHS[], Spare[];
{
    spSolve( (char *)Matrix, RHS, RHS, (spREAL*)NULL, (spREAL*)NULL );
}

/*
 * SMPmatSize()
 */
int
SMPmatSize( Matrix )
SMPmatrix *Matrix;
{
    return spGetSize( (char *)Matrix, 1 );
}

/*
 * SMPnewMatrix()
 */
int
SMPnewMatrix( pMatrix )
SMPmatrix **pMatrix;
{
int Error;
    *pMatrix = (SMPmatrix *)spCreate( 0, 1, &Error );
    return Error;
}

/*
 * SMPdestroy()
 */
void
SMPdestroy( Matrix )
SMPmatrix *Matrix;
{
    spDestroy( (char *)Matrix );
}

/*
 * SMPpreOrder()
 */
int
SMPpreOrder( Matrix )
SMPmatrix *Matrix;
{
    spMNA_Preorder( (char *)Matrix );
    return spError( (char *)Matrix );
}

/*
 * SMPprint()
 */
/*ARGSUSED*/
void
SMPprint( Matrix, File )
SMPmatrix *Matrix;
FILE *File;
{
    spPrint( (char *)Matrix, 0, 1, 1 );
}

/*
 * SMPgetError()
 */
void
SMPgetError( Matrix, Col, Row)
SMPmatrix *Matrix;
int *Row, *Col;
{
    spWhereSingular( (char *)Matrix, Row, Col );
}

/*
 * SMPcProdDiag()
 */
int
SMPcProdDiag( Matrix, pMantissa, pExponent)
SMPmatrix *Matrix;
SPcomplex *pMantissa;
int *pExponent;
{
    spDeterminant( (char *)Matrix, pExponent, &(pMantissa->real),
                                              &(pMantissa->imag) );
    return spError( (char *)Matrix );
}

/*
 *  LOAD GMIN
 *
 *  This routine adds Gmin to each diagonal element.  Because Gmin is
 *  added to the current diagonal, which may bear little relation to
 *  what the outside world thinks is a diagonal, and because the
 *  elements that are diagonals may change after calling spOrderAndFactor,
 *  use of this routine is not recommended.  It is included here simply
 *  for compatibility with Spice3.
 */
#include "spDefs.h"

LoadGmin( eMatrix, Gmin )

char *eMatrix;
register double Gmin;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;
register int I;
register ArrayOfElementPtrs Diag;

/* Begin `spLoadGmin'. */
    ASSERT( IS_SPARSE( Matrix ) );

    Diag = Matrix->Diag;
    for (I = Matrix->Size; I > 0; I--)
        Diag[I]->Real += Gmin;
    return;
}
