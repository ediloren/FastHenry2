/*
 *  MATRIX COMPATIBILITY MODULE
 *
 *  Author:                     Advising professor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This module contains routines that make Sparse1.3 upward compatible
 *  from Sparse1.2.  These routines are not suggested for use in new
 *  software.
 *
 *  >>> User accessible functions contained in this file:
 *  AllocateMatrix
 *  DeallocateMatrix
 *  CleanMatrix
 *  ClearMatrix
 *  AddElementToMatrix
 *  AddRealElementToMatrix
 *  AddImagElementToMatrix
 *  AddComplexElementToMatrix
 *  AddAdmittanceToMatrix
 *  AddQuadToMatrix
 *  AddOnesToMatrix
 *  AddRealQuadElementToMatrix
 *  AddImagQuadElementToMatrix
 *  AddComplexQuadElementToMatrix
 *  OrderAndDecomposeMatrix
 *  DecomposeMatrix
 *  SolveMatrix
 *  SolveTransposedMatrix
 *  DeleteRowAndColFromMatrix
 *  PrintMatrix
 *  OutputMatrixToFile
 *  PreorderForModifiedNodal
 *  ScaleMatrix
 *  MultiplyMatrix
 *  Determinant
 *  MatrixRoundoffError
 *  MatrixError
 *  ClearMatrixError
 *  GetMatrixSize
 *  SetMatrixReal
 *  SetMatrixComplex
 *  MatrixFillinCount
 *  MatrixElementCount
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
    "@(#)$Header: spCompat.c,v 1.3 88/06/24 05:00:47 kundert Exp $";
#endif




/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spConfig.h
 *     Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *     Macros and declarations to be imported by the user.
 *  spDefs.h
 *     Matrix type and macro definitions for the sparse matrix routines.
 */

#define spINSIDE_SPARSE
#include "spConfig.h"
#if spCOMPATIBILITY
#include "spMatrix.h"
#include "spDefs.h"

#ifdef ultrix
extern void perror();
#endif



/*
 * FILE NAME DEFINITIONS
 */

#define  MATRIX_FILE            "./matrix"
#define  STATISTICS_FILE        "./stats"




/*
 * ALLOCATE MATRIX
 */

char *
AllocateMatrix( Size, Complex, pError )

int Size, *pError;
BOOLEAN Complex;
{
/* Begin `AllocateMatrix'. */
    return spCreate( Size, Complex, pError );
}





/*
 * DESTROY MATRIX
 */

void
DeallocateMatrix( Matrix )

char *Matrix;
{
/* Begin `DeallocateMatrix'. */
    spDestroy( Matrix );
    return;
}





#if STRIP
/*
 * CLEAN MATRIX
 */

void
CleanMatrix( Matrix )

char *Matrix;
{
/* Begin `CleanMatrix'. */
    spStripFills( Matrix );
    return;
}
#endif





/*
 * CLEAR MATRIX
 */

void
ClearMatrix( Matrix )

char *Matrix;
{
/* Begin `ClearMatrix'. */
    spClear( Matrix );
    return;
}






/*
 *  SINGLE ELEMENT ADDITION TO MATRIX BY INDEX
 */

/*VARARGS4*/

RealNumber *
AddElementToMatrix( Matrix, Row, Col, Real, Imag )

char *Matrix;
int  Row, Col;
RealNumber  Real, Imag;
{
RealNumber  *pElement;

/* Begin `AddElementToMatrix'. */
    pElement = spGetElement( Matrix, Row, Col );

/* Add Data to element. */
    spADD_REAL_ELEMENT( pElement, Real );
#if spCOMPLEX
    if (((MatrixPtr)Matrix)->Complex)
        spADD_IMAG_ELEMENT( pElement, Imag );
#endif

    return pElement;
}







/*
 *  ADD ELEMENT TO MATRIX BY POINTER ROUTINES
 *
 *  These routines add elements to the matrix using pointers.  There use
 *  is discouraged if the calling program is written in C.  If that is the
 *  case, then the equivalent macros should be used.  These routines are
 *  functionally equivalent to the macros, but will be slower.
 *
 *  >>> Arguments:
 *  Element  <input>  (RealNumber *)
 *     Pointer to the element that the data is to be added to.
 *  Real  <input>  (RealNumber)
 *     Real data to be added to elements.
 *  Imag  <input>  (RealNumber)
 *     Imag data to be added to elements.  If matrix is real, this argument
 *     may be deleted.
 */

void
AddRealElementToMatrix( Element, Real )

RealNumber *Element;
RealNumber  Real;
{
/* Begin `AddRealElementToMatrix'. */

   spADD_REAL_ELEMENT( Element, Real );
   return;
}



#if spCOMPLEX

void
AddImagElementToMatrix( Element, Imag )

RealNumber *Element;
RealNumber  Imag;
{
/* Begin `AddImagElementToMatrix'. */

   spADD_IMAG_ELEMENT( Element, Imag );
   return;
}




void
AddComplexElementToMatrix( Element, Real, Imag )

RealNumber *Element;
RealNumber  Real, Imag;
{
/* Begin `AddComplexElementToMatrix'. */

   spADD_COMPLEX_ELEMENT( Element, Real, Imag );
   return;
}

#endif








#if QUAD_ELEMENT
/*
 *  ADDITION OF ADMITTANCE TO MATRIX BY INDEX
 */

/*VARARGS5*/

void
AddAdmittanceToMatrix( Matrix, Node1, Node2, Template, Real, Imag )

char  *Matrix;
int  Node1, Node2;
struct  spTemplate  *Template;
RealNumber  Real, Imag;
{
/* Begin `AddAdmittanceToMatrix'. */
    Template->Element1 = spGetElement( Matrix, Node1, Node1 );
    Template->Element2 = spGetElement( Matrix, Node2, Node2 );
    Template->Element3Negated = spGetElement( Matrix, Node2, Node1 );
    Template->Element4Negated = spGetElement( Matrix, Node1, Node2 );

    if (Node1 == 0)
        SWAP( RealNumber *, Template->Element1, Template->Element2 );

/* Add Data to elements. */
    spADD_REAL_QUAD( *Template, Real );
#if spCOMPLEX
    if (((MatrixPtr)Matrix)->Complex)
        spADD_IMAG_QUAD( *Template, Imag );
#endif

    return;
}
#endif /* QUAD_ELEMENT */









#if QUAD_ELEMENT
/*
 *  ADDITION OF FOUR ELEMENTS TO MATRIX BY INDEX
 */

/*VARARGS7*/

void
AddQuadToMatrix( Matrix, Row1, Row2, Col1, Col2, Template, Real, Imag )

char  *Matrix;
int  Row1, Row2, Col1, Col2;
struct  spTemplate  *Template;
RealNumber  Real, Imag;
{
/* Begin `AddQuadToMatrix'. */
    Template->Element1 = spGetElement( Matrix, Row1, Col1 );
    Template->Element2 = spGetElement( Matrix, Row2, Col2 );
    Template->Element3Negated = spGetElement( Matrix, Row2, Col1 );
    Template->Element4Negated = spGetElement( Matrix, Row1, Col2 );

    if (Template->Element1 == &((MatrixPtr)Matrix)->TrashCan.Real)
        SWAP( RealNumber *, Template->Element1, Template->Element2 );

/* Add Data to elements. */
    spADD_REAL_QUAD( *Template, Real );
#if spCOMPLEX
    if (((MatrixPtr)Matrix)->Complex)
        spADD_IMAG_QUAD( *Template, Imag );
#endif

    return;
}
#endif /* QUAD_ELEMENT */









#if QUAD_ELEMENT
/*
 *  ADDITION OF FOUR STRUCTURAL ONES TO MATRIX BY INDEX
 *
 *  Performs similar function to AddQuadToMatrix except this routine is
 *  meant for components that do not have an admittance representation.
 *
 *  The following stamp is used:
 *         Pos  Neg  Eqn
 *  Pos  [  .    .    1  ]
 *  Neg  [  .    .   -1  ]
 *  Eqn  [  1   -1    .  ]
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to the matrix that component is to be entered in.
 *  Pos  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Pos be less than zero.
 *  Neg  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Neg be less than zero.
 *  Eqn  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Eqn be less than zero.
 *  Template  <output>  (struct spTemplate *)
 *     Collection of pointers to four elements that are later used to directly
 *     address elements.  User must supply the template, this routine will
 *     fill it.
 *
 *  Possible errors:
 *  NO_MEMORY
 *  RANGE
 *  Error is not cleared in this routine.
 */

void
AddOnesToMatrix(Matrix, Pos, Neg, Eqn, Template)

char  *Matrix;
int  Pos, Neg, Eqn;
struct  spTemplate  *Template;
{
/* Begin `AddOnesToMatrix'. */
    Template->Element4Negated = spGetElement( Matrix, Neg, Eqn );
    Template->Element3Negated = spGetElement( Matrix, Eqn, Neg );
    Template->Element2 = spGetElement( Matrix, Pos, Eqn );
    Template->Element1 = spGetElement( Matrix, Eqn, Pos );

    spADD_REAL_QUAD( *Template, 1.0 );
    return;
}
#endif /* QUAD_ELEMENT */












#if QUAD_ELEMENT
/*
 *  QUADRUPLE REAL ELEMENT ADDITION TO MATRIX BY POINTER
 *
 *  Adds Data to four real elements specified by the Template.
 *  It is better to use the equivalent macro.
 *
 *  >>> Arguments:
 *  Template  <input>  (struct spTemplate *)
 *     Pointer to the group of pointers to four elements.
 *  Real  <input>  (RealNumber)
 *     Real data to be added to elements.
 */

void
AddRealQuadElementToMatrix( Template, Real )

struct  spTemplate  *Template;
RealNumber  Real;
{
/* Begin `AddRealQuadElementToMatrix'. */
    spADD_REAL_QUAD( *Template, Real );
    return;
}
#endif /* QUAD_ELEMENT */











#if QUAD_ELEMENT AND spCOMPLEX
/*
 *  QUADRUPLE IMAGINARY ELEMENT ADDITION TO MATRIX BY POINTER
 *
 *  Adds data to four imaginary elements specified by the Template.
 *  It is better to use the equivalent macro.
 *
 *  >>> Arguments:
 *  Template  <input>  (struct spTemplate *)
 *     Pointer to the group of pointers to four elements.
 *  Imag  <input>  (RealNumber)
 *     Imaginary data to be added to elements.
 */

void
AddImagQuadElementToMatrix( Template, Imag )

struct  spTemplate  *Template;
RealNumber  Imag;
{
/* Begin `AddImagQuadElementToMatrix'. */
    spADD_IMAG_QUAD( *Template, Imag );
    return;
}
#endif /* QUAD_ELEMENT AND spCOMPLEX */











#if QUAD_ELEMENT AND spCOMPLEX
/*
 *  QUADRUPLE COMPLEX ELEMENT ADDITION TO MATRIX BY POINTER
 *
 *  Adds Data to four complex elements specified by the Template.
 *  It is better to use the equivalent macro.
 *
 *  >>> Arguments:
 *  Template  <input>  (struct spTemplate *)
 *     Pointer to the group of pointers to four elements.
 *  Real  <input>  (RealNumber)
 *     Real data to be added to elements.
 *  Imag  <input>  (RealNumber)
 *     Imaginary data to be added to elements.
 */

void
AddComplexQuadElementToMatrix( Template, Real, Imag )

struct  spTemplate  *Template;
RealNumber  Real, Imag;
{
/* Begin `AddComplexQuadElementToMatrix'. */
    spADD_COMPLEX_QUAD( *Template, Real, Imag );
    return;
}
#endif /* QUAD_ELEMENT AND spCOMPLEX */





/*
 *   ORDER AND DECOMPOSE MATRIX
 */

#if NOT STABILITY AND NOT PSEUDOCONDITION
/*VARARGS4*/
#endif

int
OrderAndDecomposeMatrix( eMatrix, RHS, RelThreshold, AbsThreshold,
                         Growth, PseudoCondition, LargestElement )

char *eMatrix;
RealNumber  RHS[], RelThreshold, AbsThreshold;
RealNumber  *Growth, *PseudoCondition, *LargestElement;
{
int Error;
RealNumber LargestBefore, LargestAfter;

/* Begin `OrderAndDecomposeMatrix'. */

#if STABILITY
    LargestBefore = spLargestElement( eMatrix );
#endif
    Error = spOrderAndFactor( eMatrix, RHS, RelThreshold, AbsThreshold, YES );
    if (Error >= FATAL) return Error;
#if PSEUDOCONDITION
    if (PseudoCondition != NULL) *PseudoCondition = spPseudoCondition(eMatrix);
#endif
#if STABILITY
    LargestAfter = spLargestElement( eMatrix );
    if (Growth != NULL)
    {   if (LargestBefore != 0.0)
            *Growth = LargestAfter/LargestBefore;
        else
            *Growth = 0.0;
    }
    if (LargestElement != NULL) *LargestElement = LargestAfter;
#endif
    return Error;
}




/*
 *   DECOMPOSE MATRIX
 */

#if NOT STABILITY AND NOT PSEUDOCONDITION
/*VARARGS1*/
#endif

int
DecomposeMatrix( eMatrix, Growth, PseudoCondition, LargestElement )

RealNumber  *Growth, *PseudoCondition, *LargestElement;
char *eMatrix;
{
int Error;
RealNumber LargestBefore, LargestAfter;

/* Begin `DecomposeMatrix'. */

#if STABILITY
    LargestBefore = spLargestElement( eMatrix );
#endif
    Error = spFactor( eMatrix );
    if (Error >= FATAL) return Error;
#if PSEUDOCONDITION
    if (PseudoCondition != NULL) *PseudoCondition = spPseudoCondition(eMatrix);
#endif
#if STABILITY
    LargestAfter = spLargestElement( eMatrix );
    if (Growth != NULL)
    {   if (LargestBefore != 0.0)
            *Growth = LargestAfter/LargestBefore;
        else
            *Growth = 0.0;
    }
    if (LargestElement != NULL) *LargestElement = LargestAfter;
#endif
    return Error;
}





/*
 *   SOLVE MATRIX
 */

/*VARARGS3*//*ARGSUSED*/

void
SolveMatrix( eMatrix, RHS, Solution, iRHS, iSolution )

char *eMatrix;
RealVector  RHS, Solution, iRHS, iSolution;
{
#if spCOMPLEX AND spSEPARATED_COMPLEX_VECTORS
    spSolve( eMatrix, RHS, Solution, iRHS, iSolution );
#else
    spSolve( eMatrix, RHS, Solution );
#endif
}


#if TRANSPOSE
/*VARARGS3*//*ARGSUSED*/

void
SolveTransposedMatrix( eMatrix, RHS, Solution, iRHS, iSolution )

char *eMatrix;
RealVector  RHS, Solution, iRHS, iSolution;
{
#if spCOMPLEX AND spSEPARATED_COMPLEX_VECTORS
    spSolveTransposed( eMatrix, RHS, Solution, iRHS, iSolution );
#else
    spSolveTransposed( eMatrix, RHS, Solution );
#endif
}
#endif




#if TRANSLATE AND DELETE
/*
 *  DELETE A ROW AND COLUMN FROM THE MATRIX
 */

void
DeleteRowAndColFromMatrix( eMatrix, Row, Col )

char *eMatrix;
int  Row, Col;
{
    spDeleteRowAndCol( eMatrix, Row, Col );
    return;
}
#endif




/*
 *   PRINT MATRIX
 */

void
PrintMatrix( eMatrix, Compressed, PrintReordered )

char *eMatrix;
BOOLEAN  Compressed, PrintReordered;
{
    spPrint( eMatrix, PrintReordered, NOT Compressed, YES );
}





/*
 *  OUTPUT MATRIX TO FILE
 */

void
OutputMatrixToFile( eMatrix, Label, Reordered, Data, Header )

char *eMatrix, *Label;
BOOLEAN Reordered, Data, Header;
{
char ErrMsg[BUFSIZ];

/* Begin `OutputMatrixToFile'. */

    if (NOT spFileMatrix( eMatrix, MATRIX_FILE, Label, Reordered, Data, Header))
    {   (void)sprintf(ErrMsg, "sparse: file `%s'", MATRIX_FILE);
        perror( ErrMsg );
    }
    return;
}




/*
 *  OUTPUT VECTOR TO FILE
 */
/*VARARGS*//*ARGSUSED*/

void
OutputVectorToFile( eMatrix, RHS, iRHS )

char *eMatrix;
RealVector RHS, iRHS;
{
char ErrMsg[BUFSIZ];

/* Begin `OutputVectorToFile'. */
    if (RHS == NULL) return;

    if (NOT spFileVector( eMatrix, MATRIX_FILE, RHS IMAG_RHS ))
    {   (void)sprintf(ErrMsg, "sparse: file `%s'", MATRIX_FILE);
        perror( ErrMsg );
    }
    return;
}




/*
 *  OUTPUT STATISTICS TO FILE
 */

void
OutputStatisticsToFile( eMatrix, Label )

char *eMatrix, *Label;
{
char ErrMsg[BUFSIZ];
extern int pFileStats;

/* Begin `OutputStatisticsToFile'. */

    if (NOT spFileStats( eMatrix, STATISTICS_FILE, Label ))
    {   (void)sprintf(ErrMsg, "sparse: file `%s'", STATISTICS_FILE);
        perror( ErrMsg );
    }
    return;
}





#if MODIFIED_NODAL
/*
 *  PREORDER MODIFIED NODE ADMITTANCE MATRIX TO REMOVE ZEROS FROM DIAGONAL
 */

void
PreorderForModifiedNodal( eMatrix )

char *eMatrix;
{
/* Begin `PreorderForModifiedNodal'. */

    spMNA_Preorder( eMatrix );
}
#endif






#if SCALING
/*
 *  SCALE MATRIX
 */

void
ScaleMatrix( eMatrix, RHS_ScaleFactors, SolutionScaleFactors )

char *eMatrix;
RealVector  RHS_ScaleFactors, SolutionScaleFactors;
{
/* Begin `ScaleMatrix'. */

    spScale( eMatrix, RHS_ScaleFactors, SolutionScaleFactors );
}
#endif





#if MULTIPLICATION
/*
 *  MATRIX MULTIPLICATION
 */
 /*VARARGS3*//*ARGSUSED*/

void
MatrixMultiply( eMatrix, RHS, Solution, iRHS, iSolution )

char *eMatrix;
RealVector RHS, Solution, iRHS, iSolution;
{
/* Begin `MatrixMultiply'. */

#if spCOMPLEX AND spSEPARATED_COMPLEX_VECTORS
    spMultiply( eMatrix, RHS, Solution, iRHS, iSolution );
#else
    spMultiply( eMatrix, RHS, Solution );
#endif
}
#endif





#if DETERMINANT
/*
 *   CALCULATE DETERMINANT
 */

void
Determinant (eMatrix, pExponent, pDeterminant, piDeterminant )

char *eMatrix;
register  RealNumber *pDeterminant, *piDeterminant;
int  *pExponent;
{
/* Begin `Determinant'. */

#if spCOMPLEX
    spDeterminant (eMatrix, pExponent, pDeterminant, piDeterminant );
#else
    spDeterminant (eMatrix, pExponent, pDeterminant );
#endif
}
#endif






#if STABILITY
/*
 *  CALCULATE ROUNDOFF ERROR ESTIMATION
 *
 *  This function is not implemented.
 */

RealNumber
MatrixRoundoffError( eMatrix )

char *eMatrix;
{
/* Begin `MatrixRoundoffError'. */
    return 0.0;
/*  return spRoundoff( eMatrix, Rho );  */
}
#endif






/*
 *  RETURN MATRIX ERROR STATUS
 */

int
MatrixError( Matrix )

char  *Matrix;
{
/* Begin `MatrixError'. */

    return spError( Matrix );
}






/*
 *  CLEAR MATRIX ERROR FLAG
 */

int
ClearMatrixError( Matrix )

char  *Matrix;
{
int  Error;

/* Begin `ClearMatrixError'. */
    if (Matrix == NULL) return RANGE;
    Error = ((MatrixPtr)Matrix)->Error;
    ((MatrixPtr)Matrix)->Error = spOKAY;
    return Error;
}






/*
 *   MATRIX SIZE
 */

int
GetMatrixSize( Matrix, External )

char  *Matrix;
BOOLEAN  External;
{
/* Begin `GetMatrixSize'. */

    return spGetSize( Matrix, External );
}







/*
 *   SET MATRIX COMPLEX OR REAL
 */

void
SetMatrixReal( Matrix )

char *Matrix;
{
/* Begin `SetMatrixReal'. */

    spSetReal( Matrix );
}


void
SetMatrixComplex( Matrix )

char  *Matrix;
{
/* Begin `SetMatrixComplex'. */
    spSetComplex( Matrix );
}








/*
 *   ELEMENT OR FILL-IN COUNT
 */

int
MatrixFillinCount( Matrix )

char *Matrix;
{
/* Begin `MatrixFillinCount'. */
    return spFillinCount( Matrix );
}


int
MatrixElementCount( Matrix )

char  *Matrix;
{
/* Begin `MatrixElementCount'. */
    return spElementCount( Matrix );
}

#endif /* spCOMPATIBILITY */
