.LP
.pn 1
.if t .nr PS 11
.if t .nr VS 13
.if t .po 1.5i
.if t .nr PO 1.5i
.if n .ll 7.5i
.if n .nr LL 7.5i
.rs
.ls 1
.\
.\
.\  >>> Revision information:
.\  $Author: kundert $
.\  $Date: 88/06/26 05:01:56 $
.\  $Revision: 1.4 $
.\
.\
.\  Define registers for section and subsection
.nr sc 0 1
.nr ss 0 1
.\
.\
.\  Define macro for separating procedure descriptions
.de Lb 
.KS
.nf
.ta 6i
.tc \(ru
	
.tc
.LP
.sp 0.3v
.ta 6iR
\s+2\fB\\n(sc.\\n+(ss:\ \ \\$1(\|)\fR\s-2       \"	\s+4\fI\\$1()\fR\s-4
.ta 0.5i
.LP
..
.de Le 
.LP
.ta 6i
.tc \(ru
	
.tc
.KE
.sp 0.3v
..
.EQ
delim $$
gsize 11
tdefine norm % \^ \(br fwd 10 \(br %
ndefine norm % || %
ndefine delta % d %
ndefine epsilon % e %
ndefine kappa % K %
.EN
.hw SEP-AR-ATED_-COM-PLEX_-VEC-TORS
.\
.\
.\  Begin the document

.ps 24
.sp 1.75i
.ce 
.B "Sparse User's Guide"
.ps
.sp .2i
.LG
.LG
.LG
.I
.ce 1
A Sparse Linear Equation Solver
.NL
.NL
.R
.sp 0.25i
.ce
Version 1.3a
.sp .05i
.ce 1
1 April 1988
.NL
.sp 0.75i
.R
.LG
.R
.ce 2
.I "Kenneth S. Kundert"
.I "Alberto Sangiovanni-Vincentelli"
.sp 1.0i
.ce 4
Department of 
.sp 0.25v
Electrical Engineering and Computer Sciences
.sp 0.25v
University of California, Berkeley
.sp 0.25v
Berkeley, Calif. 94720
.NL
.bp
.RT
.ne 1i
.LG
.B "\n+(sc:  INTRODUCTION"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Introduction
.XE
.PP
\fISparse1.3\fP is a flexible package of subroutines written in C used
to quickly and accurately solve large sparse systems of linear
equations.  The package is able to handle arbitrary real and complex
square matrix equations.  Besides being able to solve linear systems,
it is also able to quickly solve transposed systems, find determinants,
and estimate errors due to ill-conditioning in the system of equations
and instability in the computations.  \fISparse\fR also provides a test
program that is able read matrix equations from a file, solve them, and
print useful information about the equation and its solution.
.PP
\fISparse1.3\fR is generally as fast or faster than other popular
sparse matrix packages when solving many matrices of similar
structure.  \fISparse\fR does not require or assume symmetry and is
able to perform numerical pivoting to avoid unnecessary error in the
solution.  It handles its own memory allocation, which allows the user
to forgo the hassle of providing adequate memory.  It also has a
natural, flexible, and efficient interface to the calling program.
.PP
\fISparse\fR was originally written for use in circuit simulators and
is particularly apt at handling node- and modified-node admittance
matrices.  The systems of linear generated in a circuit simulator stem
from solving large systems of nonlinear equations using Newton's method
and integrating large stiff systems of ordinary differential
equations.  However, \fISparse\fR is also suitable for other uses, one
in particular is solving the very large systems of linear equations
resulting from the numerical solution of partial differential
equations.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Features of Sparse1.3"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Features of Sparse1.3
.XE
.PP
Beyond the basic capability of being able to create, factor and 
solve systems of equations, this package features several other capabilities
that enhance its utility.  These features are:
.IP \(bu
Ability to handle both real and complex systems of equations.  Both types
may resident and active at the same time.  In fact, the same matrix
may alternate between being real and complex.
.IP \(bu
Ability to quickly solve the transposed system.  This feature is useful
when computing the sensitivity of a circuit using the adjoint method.
.IP \(bu
Memory for elements in the matrix is allocated dynamically, so the
size of the matrix is only limited by the amount of memory available
to \fISparse\fP and the range of the integer data type, which is used to
hold matrix indices.
.IP \(bu
Ability to efficiently compute the condition number of the matrix and
an a posteriori estimate of the error caused by growth in the size of
the elements during the factorization.
.IP \(bu
Much of the matrix initialization can be performed by \fISparse\fP,
providing advantages in speed and simplified coding of the calling
program.
.IP \(bu
Ability to preorder modified node admittance matrices to enhance 
accuracy and speed.
.IP \(bu
Ability to exploit sparsity in the right-hand side vector to 
reduce unnecessary computation.
.IP \(bu
Ability to scale matrices prior to factoring to reduce uncertainty 
in the solution.
.IP \(bu
The ability to create and build a matrix without knowing its final size.
.IP \(bu
The ability to add elements, and rows and columns,
to a matrix after the matrix has been reordered.  
.IP \(bu
The ability to delete rows and columns from a matrix.  
.IP \(bu
The ability to strip the fill-ins from a matrix.  This can improve the
efficiency of a subsequent reordering.
.IP \(bu
The ability to handle matrices that have rows and columns missing from
their input description.  
.IP \(bu
Ability to output the matrix in forms readable by either by people or
by the \fISparse\fP package.  Basic statistics on the matrix can also be
output.
.IP \(bu
By default all arithmetic operations and number storage use double precision.
Thus, \fISparse\fP usually gives accurate results, even on highly
ill-conditioned systems.  If so desired, \fISparse\fP can be easily
configured to use single precision arithmetic.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Enhancements of Sparse1.3 over Sparse1.2"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Enhancements of Sparse1.3 over Sparse1.2
.XE
.PP
Most notable of the enhancements provided by \fISparse1.3\fR is that
it is considerably faster on dense matrices.  Also, external names
have been made unique to 7 characters and the \fISparse\fR prefix
\fBsp\fR has been prepended to all externally accessible names to
avoid conflicts.  In addition, a routine that efficiently estimates
the condition number of a matrix has been added and the code that
estimates the growth in the factorization has been split off from the
actual factorization so that it is computed only when needed.
.PP
It is now possible for the user program to store information in the
matrix elements.  It is also possible to provide a subroutine to
\fISparse\fP that uses that information to initialize the matrix.  This
can greatly simplify the user's code.
.PP
Though the interface between \fISparse1.3\fR and the calling program
has changed considerable from previous versions of \fISparse\fR, it is
possible to compile additional code that provides backward compatibility
to \fISparse1.2\fR at the expense of a slight loss of efficiency by
setting a compiler option.  \fISparse1.3\fR now also has an FORTRAN
interface.  Routines written in FORTRAN can access almost all of the
features \fISparse1.3\fR.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Copyright Information"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Copyright Information
.XE
.PP
\fISparse1.3\fP has been copyrighted.  Permission to use, copy, modify,
and distribute this software and its documentation without fee is hereby
granted, provided that the copyright notice appear in all copies and
that the authors, \fISparse\fP, and the University of California at
Berkeley are referenced in all documentation for the program or product
in which \fISparse\fP is to be installed.  The authors and the
University of California make no representations as to the suitability
of the software for any purpose.  It is provided `as is', without
express or implied warranty.
.bp
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  PRIMER"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Primer
.XE
.sp
.RT
.LG
.B "\n(sc.\n+(ss:  Solving Matrix Equations"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Solving Matrix Equations
.XE
.PP
\fISparse\fP contains a collection of C subprograms that can be used to solve
linear algebraic systems of equations.  These systems are of the form:
.EQ I
bold { Ax ~=~ b }
.EN
where $bold A$ is an $n times  n$ matrix, $bold x$ is the vector of $n$
unknowns and $bold b$ is the vector of $n$ right-hand side terms.
Through out this package $bold A$ is denoted \fIMatrix\fP, $bold x$ is
denoted \fISolution\fP and $bold b$ is denoted \fIRHS\fP (for
right-hand side).  The system is solved using LU factorization, so the
actual solution process is broken into two steps, the factorization or
decomposition of the matrix, performed by \fIspFactor()\fP, and the
forward and backward substitution, performed by \fIspSolve()\fP.
\fIspFactor()\fP factors the given matrix into upper and lower
triangular matrices independent of the right-hand side.  Once this is
done, the solution vector can be determined efficiently for any number
of right-hand sides without refactoring the matrix.
.PP
This package exploits the fact that large matrices usually are sparse
by not storing or operating on elements in the matrix that are zero.
Storing zero elements is avoided by organizing the matrix into an
orthogonal linked-list.  Thus, to access an element if only its indices
are known requires stepping through the list, which is slow.  This
function is performed by the routine \fIspGetElement()\fP.  It is used to
initially enter data into a matrix and to build the linked-list.
Because it is common to repeatedly solve matrices with identical
zero/nonzero structure, it is possible to reuse the linked-list.  Thus,
the linked list is left in memory and the element values are simply
cleared by \fIspClear()\fP before the linked-list is reused.  To speed the
entering of the element values into successive matrices,
\fIspGetElement()\fP returns a pointer to the element in the matrix.
This pointer can then be used to place data directly into the matrix
without having to traverse through the linked-list.
.PP
The order in which the rows and columns of the matrix are factored is
very important.  It directly affects the amount of time required for
the factorization and the forward and backward substitution.  It also
affects the accuracy of the result.  The process of choosing this order
is time consuming, but fortunately it usually only has to be done once
for each particular matrix structure encountered.  When a matrix with a
new zero/nonzero structure is to be factored, it is done by using
\fIspOrderAndFactor()\fP.  Subsequent matrices of the same structure are
factored with \fIspFactor()\fP.  The latter routine does not have the
ability to reorder matrix, but it is considerably faster.  It may be
that a order chosen may be unsuitable for subsequent factorizations.
If this is known to be true a priori, it is possible to use
\fIspOrderAndFactor()\fP for the subsequent factorizations, with a
noticeable speed penalty.  \fIspOrderAndFactor()\fP monitors the numerical
stability of the factorization and will modify an existing ordering to
maintain stability.  Otherwise, an a posteriori measure of the numerical
stability of the factorization can be computed, and the matrix
reordered if necessary.
.PP
The \fISparse\fP routines allow several matrices of different
structures to be resident at once.  When a matrix of a new structure is
encountered, the user calls \fIspCreate()\fP.  This routine creates the
basic frame for the linked-list and returns a pointer to this frame.
This pointer is then passed as an argument to the other \fISparse\fP
routines to indicate which matrix is to be operated on.  The number of
matrices that can be kept in memory at once is only limited by the
amount of memory available to the user and the size of the matrices.
When a matrix frame is no longer needed, the memory can be reclaimed by
calling \fIspDestroy()\fP.
.PP
A more complete discussion of sparse systems of equations, methods
for solving them, their error mechanisms, and the algorithms used in
\fISparse\fP can be found in Kundert [kundert86].  A particular emphasis
is placed on matrices resulting from circuit simulators.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Error Control"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Error Control
.XE
.PP
There are two separate mechanisms that can cause errors during the
factoring and solution of a system of equations.  The first is
ill-conditioning in the system.  A system of equations is
ill-conditioned if the solution is excessively sensitive to
disturbances in the input data, which occurs when the system is nearly
singular.  If a system is ill-conditioned then uncertainty in the
result is unavoidable, even if $bold A$ is accurately factored into
$bold L$ and $bold U$.  When ill-conditioning is a problem, the problem
as stated is probably ill-posed and the system should be reformulated
such that it is not so ill-conditioned.  It is possible to measure the
ill-conditioning of matrix using \fIspCondition()\fP.  This function
returns an estimate of the reciprocal of the condition number of the
matrix ($kappa ( bold A )$) [strang80].  The condition number can be
used when computing a bound on the error in the solution using the
following inequality [golub83].
.EQ
{norm delta bold x norm} over {norm bold x norm} ~<=~ kappa ( bold A)
left ( {norm delta bold A norm} over {norm bold A norm} ~+~
{norm delta bold b norm} over {norm bold b norm} right ) ~+~ "higher order terms"
.EN
where $delta bold A$ and $delta bold b$ are the uncertainties in the
matrix and right-hand side vector and are assumed small.
.PP
The second mechanism that causes uncertainty is the build up of
roundoff error.  Roundoff error can become excessive if there is
sufficient growth in the size of the elements during the
factorization.  Growth is controlled by careful pivoting.  In
\fISparse\fP, the pivoting is controlled by the 
relative threshold parameter.  In conventional full matrix techniques the
pivot is chosen to be the largest element in a column.  When working 
with sparse matrices it is important to choose pivots to minimize the
reduction in sparsity.  The best pivot to retain sparsity is often not the
best pivot to retain accuracy.  Thus, some compromise must be made.  In
threshold pivoting, as used in this package, the best pivot to retain 
sparsity is used unless it is smaller than the relative threshold
times the largest element in the column.  Thus, a relative threshold close 
to one emphasizes accuracy so it will produce a minimum
amount of growth, unfortunately it also slows the factorization.  A
very small relative threshold emphasizes maintenance of sparsity and so
speeds the factorization, but can result in a
large amount of growth.  In our experience, we have found that a relative
threshold of 0.001 seems to result in a satisfactory compromise between
speed and accuracy, though other authors suggest a more conservative
value of 0.1 [duff86].
.PP
The growth that occurred during a factorization can be computed by
taking the ratio of the largest matrix element in any stage of the
factorization to the largest matrix element before factorization.  The
two numbers are estimated using \fIspLargestElement()\fP.  If the growth
is found to be excessive after \fIspOrderAndFactor()\fP, then the
relative threshold should be increased and the matrix reconstructed and
refactored.  Once the matrix has been ordered and factored without
suffering too much growth, the amount of growth that occurred should be
recorded.  If, on subsequent factorizations, as performed by
\fIspFactor()\fP, the amount of growth becomes significantly larger, then
the matrix should be reconstructed and reordered using the same
relative threshold with \fIspOrderAndFactor()\fP.  If the growth is still
excessive, then the relative threshold should be raised again.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Building the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Building the Matrix
.XE
.PP
It is not necessary to specify the size of the matrix before beginning
to add elements to it.  When the compiler option EXPANDABLE is turned
on it is possible to initially specify the size of the matrix to any
size equal to or smaller than the final size of the matrix.
Specifically, the matrix size may be initially specified as zero.  If
this is done then, as the elements are entered into the matrix, the
matrix is enlarged as needed.  This feature is particularly useful in
circuit simulators because it allows the building of the matrix as the
circuit description is parsed.  Note that once the matrix has been
reordered by the routines \fIspMNA_Preorder()\fP, \fIspFactor()\fP or
\fIspOrderAndFactor()\fP the size of the matrix becomes fixed and may no
longer be enlarged unless the compiler option TRANSLATE is enabled.
.PP
The TRANSLATE option allows \fISparse\fP to translate a non-packed set
of row and column numbers to an internal packed set.  In other words,
there may be rows and columns missing from the external description of
the matrix.  This feature provides two benefits.  First, if two
matrices are identical in structure, except for a few missing rows and
columns in one, then the TRANSLATE option allows them to be treated
identically.  Similarly, rows and columns may be deleted from a matrix
after it has been built and operated upon.  Deletion of rows and
columns is performed by the function \fIspDeleteRowAndCol()\fP.  Second,
it allows the use of the functions \fIspGetElement()\fP,
\fIspGetAdmittance()\fP, \fIspGetQuad()\fP, and \fIspGetOnes()\fP after the
matrix has been reordered.  These functions access the matrix by using
row and column indices, which have to be translated to internal indices
once the matrix is reordered.  Thus, when TRANSLATE is used in
conjunction with the EXPANDABLE option, rows and columns may be added
to a matrix after it has been reordered.
.PP
Another provided feature that is useful with circuit simulators is the
ability to add elements to the matrix in row zero or column zero.
These elements will have no affect on the matrix or the results.  The
benefit of this is that when working with a nodal formulation, grounded 
components do not have to be treated special when building the matrix.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Initializing the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Initializing the Matrix
.XE
.PP
Once a matrix has been factored, it is necessary to clear the matrix
before it can be reloaded with new values.  The straight forward way of
doing that is to call \fIspClear()\fP, which sets the value of every
element in the matrix to zero.  \fISparse\fP also provides a more
flexible way to clear the matrix.  Using \fIspInitialize()\fP, it is
possible to clear and reload at least part of the matrix in one step.
.PP
\fISparse\fR allows the user to keep initialization information with
each structurally nonzero matrix element.  Each element has a pointer
that is set and used by the user.  The user can set this pointer using
\fIspInstallInitInfo()\fR and may read it using \fIspGetInitInfo()\fR.
The function \fIspInitialize()\fR is a user customizable way to
initialize the matrix.  Passed to this routine is a function pointer.
\fIspInitialize()\fR sweeps through every element in the matrix and
checks the \fIpInitInfo\fR pointer (the user supplied pointer).  If the
\fIpInitInfo\fR is NULL, which is true unless the user changes it
(always true for fill-ins), then the element is zeroed.  Otherwise, the
function pointer is called and passed the \fIpInitInfo\fR pointer as
well as the element pointer and the external row and column numbers,
allowing the user to initialize the matrix element and the right-hand
side.
.PP
Why \fIspInitialize()\fR would be used over \fIspClear()\fR can be
illustrated by way of an example.  Consider a circuit simulator that
handles linear and nonlinear resistors and capacitors performing a
transient analysis.  For the linear resistors, a constant value is
loaded into the matrix at each time step and for each Newton
iteration.  For the linear capacitor, a value is loaded into the matrix
that is constant over Newton iterations, but is a function of the time
step and the integration method.  The nonlinear components contribute
values to the matrix that change on every time step and Newton
iteration.
.PP
\fISparse\fP allows the user to attach a data structure to each element
in the matrix.  For this example, the user might attach a structure
that held several pieces of information, such as the conductance of the
linear resistor, the capacitance of the linear capacitor, the
capacitance of the nonlinear capacitor, and perhaps past values of
capacitances.  The user also provides a subroutine to
\fIspInitialize()\fP that is called for each user-created element in the
matrix.  This routine would, using the information in the attached data
structure, initialize the matrix element and perhaps the right-hand
side vector.
.PP
In this example, the user supplied routine might load the linear
conductance into the matrix and multiply it by some voltage to find a
current that could be loaded into the right-hand side vector.  For the
capacitors, the routine would first apply an integration method and
then load the matrix and the right-hand side.
.PP
This approach is useful for two reasons.  First, much of the work of
the device code in the simulator can be off-loaded onto the matrix
package.  Since there are usually many devices, this usually results
overall in a simpler system.  Second, the integration method can be
hidden from the simulator device code.  Thus the integration method can
be changed simply by changing the routine handed to \fIspInitialize()\fP,
resulting in a much cleaner and more easily maintained simulator.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Indices"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Indices
.XE
.PP
By far the most common errors made when using \fISparse\fR are related
to array indices.  \fISparse\fR itself contributes to the problem by
having several different indexing schemes.  There are three different
options that affect index bounds or the way indices are interpreted.
The first is ARRAY_OFFSET, which only affects array indices.
ARRAY_OFFSET is a compiler flag that selects whether arrays start at
index zero or index one.  Note that if ARRAY_OFFSET is zero then
\fIRHS\fP[0] corresponds to row one in the matrix and \fISolution\fP[0]
corresponds to column one.  Further note that when ARRAY_OFFSET is set
to one, then the allocated length of the arrays handed to the
\fISparse\fR routines should be at least the external size of the matrix
plus one.  The main utility of ARRAY_OFFSET is that it allows natural
array indexing when \fISparse\fR is coupled to programs in other
languages.  For example; in FORTRAN arrays always start at one whereas
in C array always start at zero.  Thus the first entry in a FORTRAN
array corresponds to the zero'th entry in a C array.  Setting
ARRAY_OFFSET to zero allows the arrays in FORTRAN to start at one rather
than two.  For the rest of this discussion, assume that ARRAY_OFFSET is
set so that arrays start at one in the program that calls \fISparse\fR.
.sp 0.3v
.PP
The second option that affects indices is EXPANDABLE.  When EXPANDABLE
is set false the upper bound on array and matrix indices is \fISize\fP,
where \fISize\fP is a parameter handed to \fIspCreate()\fP.  When
EXPANDABLE set true, then there is essentially no upper bound on array
indices.  Indeed, the size of the matrix is determined by the largest
row or column number handed to \fISparse\fR.  The upper bound on the
array indices then equals the final size determined by \fISparse\fR.
This size can be determined by calling \fIspGetSize()\fP.
.sp 0.3v
.PP
The final option that affects indices is TRANSLATE.  This option was
provided to allow row and columns to be deleted, but it also allows row
and column numbers to be missing from the input description for a
matrix.  This means that the size of the matrix is not determined by
the largest row or column number entered into the matrix.  Rather, the
size is determined by the total number of rows or column entered.  For
example, if the elements [2,3], [5,3], and [7,2] are entered into the
matrix, the internal size of the matrix becomes four while the external
size is seven.  The internal size equals the number of rows and columns
in the matrix while the external size equals the largest row or column
number entered into the matrix.  Note that if a row is entered into the
matrix, then its corresponding column is also entered, and vice versa.
The indices used in the \fIRHS\fP and \fISolution\fP vectors
correspond to the row and column indices in the matrix.  Thus, for this
example, valid data is expected in \fIRHS\fP at locations 2, 3, 5
and 7.  Data at other locations is ignored.  Similarly, valid data is
returned in \fISolution\fP at locations 2, 3, 5, and 7.  The other
locations are left unmolested.  This shows that the length of the
arrays correspond to the external size of the matrix.  Again, this
value can be determined by \fIspGetSize()\fP.
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Configuring Sparse"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Configuring Sparse
.XE
.PP
It is possible at compile-time to customize \fISparse\fR for your
particular application.  This is done by changing the compiler options,
which are kept in the personality file, \fBspConfig.h\fP.  There are
three classes of choices available.  First are the \fISparse\fP options,
which specify the dominant personality characteristics, such as if real
and/or complex systems of equations are to be handled.  The second
class is the \fISparse\fP constants, such as the default pivot threshold and
the amount of memory initially allocated per matrix.  The last class is
the machine constants.  These numbers must be updated when \fISparse\fP
is ported to another machine.
.PP
As an aid in the setup and testing of \fISparse\fP a test routine and
several test matrices and their solutions have been provided.  The test
routine is capable of reading files generated by \fIspFileMatrix()\fP
and \fIspFileVector()\fP.
.PP
By default \fISparse\fP stores all real numbers and performs all
computations using double precision arithmetic.  This can be changed by
changing the definition of \fIspREAL\fP from \fBdouble\fP to
\fBfloat\fP.  \fIspREAL\fP is defined in \fBspExports.h\fP.
.bp
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  INTRODUCTION TO THE SPARSE ROUTINES"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Introduction to the Sparse Routines
.XE
.sp
In this section the routines are grouped by function and briefly described.

.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Creating the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Creating the Matrix
.XE
.LP
.I spCreate()
.RS
.br
Allocates and initializes the data structure for a matrix.  Necessarily the
first routine run for any particular matrix.
.RE
.KE
.KS
.LP
.I spDestroy()
.RS
.br
Destroys the data structure for a matrix and frees the memory.
.RE
.KE
.KS
.LP
.I spSetReal()
.br
.I spSetComplex()
.RS
.br
These routines toggle a flag internal to \fISparse\fP that indicates
that the matrix is either real or complex.  This is useful if both real
and complex matrices of identical structure are expected.
.RE
.KE
.sp 2
.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Building the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Building the Matrix
.XE
.LP
.I spGetElement()
.RS
.br
Assures that the specified element exists in the matrix data structure
and returns a pointer to it.
.RE
.KE
.KS
.LP
.I spGetAdmittance()
.br
.I spGetQuad()
.br
.I spGetOnes()
.RS
.br
These routines add a group of four related elements to the matrix.
\fIspGetAdmittance()\fP adds the four elements associated with a two
terminal admittance.  \fIspGetQuad()\fP is a more general routine that is
useful for entering controlled sources to the matrix.  \fIspGetOnes()\fP
adds the four structural ones to the matrix that are often encountered
with elements that do not have admittance representations.
.RE
.KE
.KS
.LP
.I spDeleteRowAndCol()
.RS
.br
This function is used to delete a row and column from the matrix.
.RE
.KE
.sp 2
.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Clearing the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Clearing the Matrix
.XE
.LP
.I spClear()
.RS
.br
Sets every element in the matrix to zero.
.RE
.KE
.KS
.LP
.I spInitialize()
.RS
.br
Runs a user provided initialization routine on each element in the matrix.
This routine would be used in lieu of \fIspClear()\fR.
.RE
.KE
.KS
.LP
.I spGetInitInfo()
.br
.I spInstallInitInfo()
.RS
.br
These routines allow the user to install and read a user-provided pointer
to initialization data for a particular matrix element.
.RE
.KE
.KS
.LP
.I  spStripFills()
.RS
.br
This routine returns a matrix to a semi-virgin state by removing all
fill-ins.  This can be useful if a matrix is to be reordered and it has
changed significantly since it was previously ordered.  This may be the
case if a few rows and columns have been added or deleted or if the
previous ordering was done on a matrix that was numerically quite
different than the matrix currently being factored.  Stripping and
reordering a matrix may speed subsequent factorization if the current
ordering is inferior, whereas simply reordering will generally only
enhance accuracy and not speed.
.RE
.KE
.sp 2
.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Placing Data in the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Placing Data in the Matrix
.XE
.LP
.I spADD_REAL_ELEMENT()
.br
.I spADD_IMAG_ELEMENT()
.br
.I spADD_COMPLEX_ELEMENT()
.RS
.br
Adds a value to a particular matrix element.
.RE
.KE
.KS
.LP
.I spADD_REAL_QUAD()
.br
.I spADD_IMAG_QUAD()
.br
.I spADD_COMPLEX_QUAD()
.RS
.br
Adds a value to a group of four matrix elements.
.RE
.KE
.sp 2
.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Influencing the Factorization"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Influencing the Factorization
.XE
.LP
.I spMNA_Preorder()
.RS
.br
This routine preorders modified node admittance matrices so that
\fISparse\fP can take full advantage of their structure.  In
particular, this routine tries to remove zeros from the diagonal so
that diagonal pivoting can be used more successfully.
.RE
.KE
.KS
.LP
.I spPartition()
.RS
.br
\fISparse\fR partitions the matrix in an attempt to make \fIspFactor()\fP
run as fast as possible.  The partitioning is a relatively expensive
operation that is not needed in all cases.  \fIspPartition()\fR allows
the user specify a simpler and faster partitioning.
.RE
.KE
.KS
.LP
.I spScale()
.br
.RS
It is sometimes desirable to scale the rows and columns of a matrix in
to achieve a better pivoting order.  This is particularly true in
modified node admittance matrices, where the size of the elements in a
matrix can easily vary through ten to twelve orders of magnitude.  This
routine performs scaling on a matrix.
.RE
.KE
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Factoring the Matrix"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Factoring the Matrix
.XE
.KS
.LP
.I spOrderAndFactor()
.RS
.br
This routine chooses a pivot order for the matrix and factors it into
LU form.  It handles both the initial factorization and subsequent
factorizations when a reordering is desired.
.RE
.KE
.KS
.LP
.I spFactor()
.RS
.br
Factors a matrix that has already been ordered by \fIspOrderAndFactor()\fR.
If \fIspFactor()\fR is passed a matrix that needs ordering, it will
automatically pass the matrix to \fIspOrderAndFactor()\fR.
.RE
.KE
.sp 2
.KS
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Solving the Matrix Equation"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Solving the Matrix Equation
.XE
.LP
.I spSolve()
.br
.RS
Solves the matrix equation
.EQ I
bold "A" bold x ~=~ bold b
.EN
given the matrix $bold A$ factored into LU form and $bold b$.
.RE
.KE
.KS
.LP
.I spSolveTransposed()
.br
.RS
When working with adjoint systems, such as in sensitivity analysis, it
is desirable to quickly solve
.EQ I
bold "A" sup roman T bold x ~=~ bold b
.EN
Once $bold A$ has been factored into LU form, this
routine can be used to solve the transposed system without having to
suffer the cost of factoring the matrix again.
.RE
.KE
.sp 2
.KS
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Numerical Error Estimation"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Numerical Error Estimation
.XE
.LP
.I spCondition()
.RS
.br
Estimates the L-infinity condition number of the matrix.  This number is
a measure of the ill-conditioning in the matrix equation.  It is also useful
for making estimates of the error in the solution.
.RE
.KE
.KS
.LP
.I spNorm()
.RS
.br
Returns the L-infinity norm (the maximum absolute row sum) of an unfactored
matrix.
.RE
.KE
.KS
.LP
.I spPseudoCondition()
.RS
.br
Returns the ratio of the largest pivot to the smallest pivot of a factored
matrix.  This is a rough indicator of ill-conditioning in the matrix.
.RE
.KE
.KS
.LP
.I spLargestElement()
.RS
.br
If passed an unfactored matrix, this routine returns the absolute value of
the largest element in the matrix.  If passed a factored matrix, it returns
an estimate of the largest element that occurred in any of the reduced
submatrices during the factorization.  The ratio of these two numbers
(factored/unfactored) is the growth, which is used to determine if the
pivoting order is numerically satisfactory.
.RE
.KE
.KS
.LP
.I spRoundoff()
.RS
.br
Returns a bound on the magnitude of the largest element in
$bold E ~=~ bold A - bold LU$, where $bold E$ represents error in
the matrix resulting from roundoff error during the factorization.
.RE
.KE
.sp 2
.KS
.RT
.LG
.B "\n(sc.\n+(ss:  Matrix Operations"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Matrix Operations
.XE
.LP
.I spDeterminant()
.RS
.br
This routine simply calculates and returns the determinant of the
factored matrix.
.RE
.KE
.KS
.LP
.I spMultiply()
.br
.RS
This routine multiplys the matrix by a vector on the
right.  This is useful for forming the product $ bold "Ax" ~=~ bold b $
in order to determine if a calculated solution is correct.
.RE
.KE
.KS
.LP
.I spMultTransposed()
.br
.RS
Multiplys the transposed matrix by a vector on the right.  This is useful
for forming the product $ bold "A sup {roman T} x" ~=~ bold b $
in order to determine if a calculated solution is correct.
.RE
.KE
.KS
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Matrix Statistics and Documentation"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Matrix Statistics and Documentation
.XE
.LP
.I spError()
.RS
.br
Determines the error status of a particular matrix.  While most of the
\fISparse\fR routines do return an indication that an error has
occurred, some do not and so \fIspError()\fP provides the only way of
uncovering these errors.
.RE
.KE
.KS
.LP
.I spWhereSingular()
.RS
.br
Returns the row and column number where the matrix was detected as 
singular or where a zero pivot was found.
.RE
.KE
.KS
.LP
.I spGetSize()
.RS
.br
A function that returns the size of the matrix.  Either the internal or
external size of the matrix is returned.  The internal size of the
matrix is the actual size of the matrix whereas the external size is
the value of the largest row or column number.  These two numbers may
differ if the TRANSLATE option is used.
.RE
.KE
.KS
.LP
.I spElementCount()
.br
.I spFillinCount()
.RS
.br
Functions that return the total number of elements in the matrix, and the
number of fill-ins in the matrix.  These functions are useful for gathering
statistics on matrices.
.RE
.KE
.KS
.LP
.I spPrint()
.RS
.br
This routine outputs the matrix as well as some statistics to standard
output in a format that is readable by people.  The matrix can be
printed in either a compressed or standard format.  In the standard
format, a numeric value is given for each structurally nonzero element,
whereas in the compressed format, only the existence or nonexistence of
an element is indicated.  This routine is not suitable for use on large
matrices.
.RE
.KE
.KS
.LP
.I spFileMatrix()
.br
.I spFileVector()
.RS
.br
These two routines send a copy of the matrix and its right-hand side
vector to a file.  This file can then be read by the test program that
is included with \fISparse\fP.  Only those elements of the matrix that
are structurally nonzero are output, so very large matrices can be sent
to a file.
.RE
.KE
.KS
.LP
.I spFileStats()
.RS
.br
This routine calculates and sends some useful statistics concerning a
matrix to a file.
.RE
.KE
.bp
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  SPARSE ROUTINES"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Routines
.XE
.sp
This section contains a complete list of the \fISparse\fR routines that are
available to the user.  Each routine is described as to its
function and how to use it.  The routines are listed in alphabetic
order.
.sp 2
.Lb spClear
Sets every element in the matrix to zero.  The \fISparse\fR error state
is cleared to spOKAY in this routine.
.sp 0.3v
\fBvoid\fI spClear( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix that is to be cleared.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspClear()\fP
.XE
.Lb spCondition
\fIspCondition()\fP computes an estimate of the condition number using a
variation on the LINPACK condition number estimation algorithm.  This
quantity is an measure of ill-conditioning in the matrix.  To avoid
problems with overflow, the reciprocal of the condition number is
returned.  If this number is small, and if the matrix is scaled such
that uncertainties in the RHS and the matrix entries are equilibrated,
then the matrix is ill-conditioned.  If the this number is near one,
the matrix is well conditioned.  This routine must only be used after a
matrix has been factored by \fIspOrderAndFactor()\fR or \fIspFactor()\fR
and before it is cleared by \fIspClear()\fR or \fIspInitialize()\fR.
.LP
Unlike the LINPACK condition number estimator, this routines returns
the L infinity condition number.  This is an artifact of Sparse placing
ones on the diagonal of the upper triangular matrix rather than the
lower.  This difference should be of no importance.
.sp 0.3v
.nf
\fIspREAL spCondition( Matrix, NormOfMatrix, Error )
.XP
\(bu Returns:
.br
An estimate of the L infinity condition number of the matrix.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix for which the condition number is desired.
.XP 
\fINormOfMatrix\fP  input  (spREAL)
.br
The L-infinity norm of the unfactored matrix as computed by \fIspNorm()\fP.
.XP 
\fIError\fP  output  (\fBint *\fP)
.br
Used to return the error code.
.RE
.XP
\(bu Possible errors:
.nf
spSINGULAR
spNO_MEMORY
Error is not cleared in this routine.
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
CONDITION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspCondition()\fP
.XE
.Lb spCreate
Allocates and initializes the data structures associated with a matrix.  This
routine is necessarily the first routine run for any particular matrix.
.sp 0.3v
.nf
\fBchar *\fIspCreate( Size, Complex, Error )\fR
.XP
\(bu Returned:
.br
A pointer to the matrix is returned cast into the form of a pointer to
a character.  This pointer is then passed and used by the other matrix
routines to refer to a particular matrix.  If an error occurs, the NULL
pointer is returned.
.XP
\(bu Arguments:
.RS
.XP
\fISize\fP  input  (\fBint\fP)
.br
Size of matrix.  When the compiler option EXPANDABLE is turned on,
\fISize\fP is used as a lower bound on the size of the matrix.
\fISize\fP must not be negative.
.XP
\fIComplex\fP  input  (\fBint\fP)
.br
Type of matrix.  If \fIComplex\fP is 0 then the matrix is real,
otherwise the matrix will be complex.  Note that if the routines are
not set up to handle the type of matrix requested, then a spPANIC error
will occur.
.XP
\fIError\fP  output  (\fBint *\fP)
.br
Returns error flag, needed because function \fIspError()\fP will not
work correctly if \fIspCreate()\fP returns NULL.
.RE
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
spPANIC
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspCreate()\fP
.XE
.Lb spDeleteRowAndCol
This function is used to delete a row and column from the matrix.  The
elements removed from the matrix are never used again and are not freed
until the matrix is destroyed and so the pointers to these elements remain
valid.
.sp 0.3v
.nf
\fBvoid\fI spDeleteRowAndCol( Matrix, Row, Col )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix from which the row and column are to be deleted.
.XP
\fIRow\fP  input  (\fBint\fP)
.br
The row to be deleted.
.XP
\fICol\fP  input  (\fBint\fP)
.br
The column to be deleted.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DELETE
TRANSLATE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspDeleteRowAndCol()\fP
.XE
.Lb spDestroy
Destroys a matrix frame and reclaims the memory.
.sp 0.3v
.nf
\fBvoid\fI spDestroy( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix frame which is to be removed from memory.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspDestroy()\fP
.XE
.Lb spDeterminant
This routine in capable of calculating the determinant of the matrix
once the LU factorization has been performed.  Hence, only use this
routine after \fIspFactor()\fP or \fIspOrderAndFactor()\fP and before
\fIspClear()\fP or \fIspInitialize()\fR.  Note that the determinants
of matrices can be very large or very small.  On large matrices, the
determinant can be far larger or smaller than can be represented by a
floating point number.  For this reason the mantissa and exponent of
the determinant are returned separately.
.sp 0.3v
.nf
\fBvoid\fI spDeterminant( Matrix, Exponent, Determinant )\fR
.br
\fBvoid\fI spDeterminant( Matrix, Exponent, Determinant, iDeterminant )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix for which the determinant is desired.
.XP
\fIExponent\fP  output  (\fBint *\fP)
.br
The logarithm base 10 of the scale factor for the determinant.  To find
the actual determinant, Exponent should be added to the exponent of 
\fIDeterminant\fP and \fIiDeterminant\fP.
.XP
\fIDeterminant\fP  output  (spREAL \fB*\fR)
.br
The real portion of the determinant.  If the matrix is real, then the 
magnitude of this number is scaled to be greater than or equal to 1.0 
and less than 10.0. Otherwise the magnitude of the complex determinant
will be scaled such.
.XP
\fIiDeterminant\fP  output  (spREAL \fB*\fR)
.br
The imaginary portion of the determinant.  When the matrix is real
this pointer need not be supplied; nothing will be returned.   
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DETERMINANT
.XP
\(bu Bugs:
.br
The sign of determinant may be in error if rows and columns have been added
or deleted from matrix.
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspDeterminant()\fP
.XE
.Lb spElementCount
Returns the total number of structurally nonzero elements in the matrix.
.sp 0.3v
.nf
\fBint\fI spElementCount( Matrix )\fR
.XP
\(bu Returns:
.br
The total number of structurally nonzero elements.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspElementCount()\fP
.XE
.Lb spError
This function returns the error status of a matrix.
.sp 0.3v
.nf
\fBint\fI MatrixError( Matrix )\fR
.XP
\(bu Returned:
.br
The error status of the given matrix.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix for which the error status is desired.
.RE
.XP
\(bu Possible errors:
.nf
spOKAY
spILL_CONDITIONED
spZERO_PIVOT
spSINGULAR
spNO_MEMORY
spPANIC
Error is not cleared in this routine.
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspError()\fP
.XE
.Lb spFactor
This routine factors the matrix into LU form and is the companion
routine to \fIspOrderAndFactor()\fP.  Unlike \fIspOrderAndFactor()\fP,
\fIspFactor()\fP cannot change the ordering.  Its utility is that it is
considerably faster.  The standard way to use these two routines is to
first use \fIspOrderAndFactor()\fP for the initial factorization.  For
subsequent factorizations, \fIspFactor()\fP is used.  If
\fIspFactor()\fP is called for the initial factorization of the matrix,
then it will automatically call \fIspOrderAndFactor()\fP with the
default thresholds.  If \fIspFactor()\fP finds a zero on the diagonal,
it will terminate early and complain.  This does not necessarily mean
that matrix is singular.  Before a matrix is condemned as being
singular, it should be run through \fIspOrderAndFactor()\fP, which can
reorder the matrix and remove the offensive zero from the diagonal.
.sp 0.3v
.nf
\fBint\fI spFactor( Matrix )\fR
.XP
\(bu Returned:
.br
The error code is returned.  Possible errors are listed below.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix to be factored.
.RE
.XP
\(bu Possible errors:
.nf
spZERO_PIVOT
spNO_MEMORY
spSINGULAR
spILL_CONDITIONED
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspFactor()\fP
.XE
.Lb spFileMatrix
Writes matrix to file in format suitable to be read back in by the
matrix test program.  Normally, \fIspFileMatrix()\fP should be executed
before the matrix is factored, otherwise matrix is output in factored
form.  If the matrix is sent to a file without the header or data, it
will be in a form that is easily plotted by typical plotting programs.
.sp 0.3v
.nf
\fBint\fI spFileMatrix( Matrix, File, Label, Reordered, Data, Header )\fR
.XP
\(bu Returns:
.br
One is returned if routine was successful, otherwise zero is returned.
The calling function can query \fIerrno\fR (the system global error
variable) as to the reason why this routine failed.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix that is to be sent to file.
.XP
\fIFile\fP  input  (\fBchar *\fP)
.br
Name of output file.
.XP
\fILabel\fP  input  (\fBchar *\fP)
.br
String that is transferred to file and used as a label.  String should
fit on one line and have no embedded line feeds.
.XP
\fIReordered\fP  input  (\fBint\fP)
.br
Specifies whether the matrix should be output using the original order
or in reordered form.  Zero specifies original order.
.XP
\fIData\fP  input  (\fBint\fP)
.br
Indicates that the element values should be output along with the
indices for each element.  Element values are not output if \fIData\fP
is zero.  This parameter must be nonzero if matrix is to be read by the
\fISparse\fR test program.
.XP
\fIHeader\fP  input  (\fBint\fP)
.br
If nonzero a header is output that includes that size of the matrix and
the label.  This parameter must be nonzero if matrix is to be read by
the \fISparse\fP test program.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DOCUMENTATION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspFileMatrix()\fP
.XE
.Lb spFileStats
Appends useful information concerning the matrix to the end of a file.
If file does not exist, it is created.  This file should not be the same
as one used to hold the matrix or vector if the matrix is to be read by
the \fISparse\fP test program.  Should be executed after the matrix is
factored.
.sp 0.3v
.nf
\fBint\fI spFileStats( Matrix, File, Label )\fR
.XP
\(bu Returns:
.br
One is returned if routine was successful, otherwise zero is returned.
The calling function can query \fIerrno\fR (the system global error
variable) as to the reason why this routine failed.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix for which statistics are desired.
.XP
\fIFile\fP  input  (\fBchar *\fP)
.br
Name of output file.
.XP
\fILabel\fP  input  (\fBchar *\fP)
.br
String that is transferred to file and is used as a label. String
should fit on one line and have no embedded line feeds.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DOCUMENTATION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspFileStats()\fP
.XE
.Lb spFileVector
Appends the \fIRHS\fP vector to the end of a file in a format suitable
to be read back in by the matrix test program.  If file does not exist,
it is created.  To be compatible with the test program, if
\fIspFileVector()\fP is run, it must be run after \fIspFileMatrix()\fP
and use the same file.
.sp 0.3v
.nf
\fBint\fI spFileVector( Matrix, File, RHS )\fR
.br
\fBint\fI spFileVector( Matrix, File, RHS, iRHS )\fR
.XP
\(bu Returns:
.br
One is returned if routine was successful, otherwise zero is returned.
The calling function can query \fIerrno\fR (the system global error
variable) as to the reason why this routine failed.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix that corresponds to the vector to be output.
.XP
\fIFile\fP  input  (\fBchar *\fP)
.br
Name of file where output is to be written.
.XP
\fIRHS\fP  input  (spREAL[])
.br
The right-hand side vector.  \fIRHS\fP contains only the real portion
of the right-hand side vector if the matrix is complex and
spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiRHS\fP  input  (spREAL[])
.br
Right-hand side vector, imaginary portion.  Not necessary if matrix is
real or if spSEPARATED_COMPLEX_VECTORS is set false.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DOCUMENTATION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspFileVector()\fP
.XE
.Lb spFillinCount
Returns the total number of fill-ins in the matrix.  A fill-in is an
element that is originally structurally zero, but becomes nonzero during the
factorization.
.sp 0.3v
.nf
\fBint\fI spFillinCount( Matrix )\fR
.XP
\(bu Returns:
.br
The total number of fill-ins.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspFillinCount()\fP
.XE
.Lb spGetAdmittance
Performs same function as \fIspGetElement()\fP except rather than one
element, all four matrix elements for a floating admittance are
reserved.  This routine also works if the admittance is grounded (zero
is the ground node).  This function returns a group of pointers to the
four elements through \fITemplate\fP, which is an output.  They are
used by the \fIspADD_QUAD()\fP macros to directly access matrix
elements during subsequent loads of the matrix.
\fIspGetAdmittance()\fR arranges the pointers in \fITemplate\fR so that
the \fIspADD_QUAD()\fR routines add the admittance to the elements at
[\fINode1,Node1\fP] and [\fINode2,Node2\fP] and subtract the admittance
from the elements at [\fINode1,Node2\fP] and [\fINode2,Node1\fP]. This
routine is only to be used before \fIspMNA_Preorder()\fP,
\fIspFactor()\fP or \fIspOrderAndFactor()\fP unless the compiler flag
TRANSLATE is enabled.
.sp 0.3v
.nf
\fBint\fI spGetAdmittance( Matrix, Node1, Node2, Template )\fR
.XP
\(bu Returned:
.br
The error code is returned.  Possible errors are listed below.
\fIspGetAdmittance()\fR does not clear the error state, so it is
possible to ignore the return code of each \fIspGetAdmittance()\fR
call, and check for errors after constructing the whole matrix by
calling \fIspError()\fR.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix that admittance is to be installed.
.XP
\fINode1\fP  input  (\fBint\fP)
.br
One node number for the admittance.  \fINode1\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler
flags are set true.  In either case \fINode1\fP must not be negative.
.XP
\fINode2\fP  input  (\fBint\fP)
.br
Other node number for the admittance.  \fINode2\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler
flags are set true.  In either case \fINode2\fP must not be negative.
.XP
\fITemplate\fP  output  (\fBstruct \fIspTemplate\fB *\fR)
.br
Collection of pointers to four elements that are later used to directly
address elements.  User must supply the template, this routine will
fill it.
.RE
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
Error is not cleared in this routine.
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
QUAD_ELEMENT
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetAdmittance()\fP
.XE
.Lb spGetElement
Reserves an element at [\fIRow,Col\fP] and returns a pointer to it.  If
element is not found then it is created and spliced into matrix.  A
pointer to the real portion of the element is returned.  This pointer
is later used by the \fIspADD_ELEMENT()\fR macros to directly access
the element.  This routine is only to be used before
\fIspMNA_Preorder()\fP, \fIspFactor()\fP or \fIspOrderAndFactor()\fP
unless the compiler option TRANSLATE is set true.
.sp 0.3v
.nf
spREAL \fB*\fIspGetElement( Matrix, Row, Col )\fR
.XP
\(bu Returned:
.br
Returns a pointer to the element.  This pointer is then used to
directly access the element during successive builds.  Returns NULL if
insufficient memory is available.  \fIspGetElement()\fR does not clear
the error state, so it is possible to ignore the return code of each
\fIspGetElement()\fR call, and check for errors after constructing the
whole matrix by calling \fIspError()\fR.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix that the element is to be added to. 
.XP
\fIRow\fP  input  (\fBint\fP)
.br
Row index for element. \fIRow\fP must be in the range [0..\fISize\fR]
unless either the TRANSLATE or EXPANDABLE compiler flags are set true.
In either case \fIRow\fP must not be negative though it may be zero.
If zero then the element is not entered into the matrix, but is otherwise
treated normally.  
.XP
\fICol\fP  input  (\fBint\fP)
.br
Column index for element. \fICol\fP must be in the range [0..\fISize\fR]
unless either the TRANSLATE or EXPANDABLE compiler flags are set true.
In either case \fICol\fP must not be negative though it may be zero.
If zero then the element is not entered into the matrix, but is
otherwise treated normally.
.RE
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
Error is not cleared in this routine.
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetElement()\fP
.XE
.Lb spGetInitInfo
With the INITIALIZE compiler option enabled \fISparse\fP allows the
user to keep initialization information with each structurally nonzero
matrix element.  Each element has a pointer (referred to as
\fIpInitInfo\fR) that is set and used by the user.  This routine returns
\fIpInitInfo\fR from a particular matrix element.
.sp 0.3v
.nf
\fBchar *\fIspGetInitInfo( pElement )\fR
.XP
\(bu Returned:
.br
The user installed pointer \fIpInitInfo\fR.
.XP
\(bu Argument:
.RS
.XP
\fIpElement\fR  input  (spREAL \fB*\fR)
.br
Pointer to the element to which \fIpInitInfo\fR is attached.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
INITIALIZE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetInitInfo()\fP
.XE
.Lb spGetOnes
Performs a similar function to \fIspGetAdmittance()\fP except that the
four reserved matrix elements are assumed to be structural ones
generated by components without admittance representations during a
modified nodal analysis.  Positive ones are placed at [\fIPos,Eqn\fP]
and [\fIEqn,Pos\fP] and negative ones are placed at [\fINeg,Eqn\fP] and
[\fIEqn,Neg\fP].  This function returns a group of pointers to the four
elements through \fITemplate\fP, which is an output.  They are used by
the \fIspADD_QUAD()\fP macros to add the ones directly to the matrix
elements during subsequent loads of the matrix.  This routine is only
to be used before \fIspMNA_Preorder()\fP, \fIspFactor()\fP or
\fIspOrderAndFactor()\fP unless the compiler flag TRANSLATE is set
true.
.sp 0.3v
.nf
\fBint\fI spGetOnes( Matrix, Pos, Neg, Eqn, Template )\fR
.XP
\(bu Returned:
.br
The error code is returned.  Possible errors are listed below.
\fIspGetOnes()\fR does not clear the error state, so it is possible to
ignore the return code of each \fIspGetOnes()\fR call, and check for
errors after constructing the whole matrix by calling \fIspError()\fR.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix that ones are to be entered in.
.XP
\fIPos\fP  input  (\fBint\fP)
.br
Number of positive node.  Must be in the range of [0..\fISize\fR]
unless either the options EXPANDABLE or TRANSLATE are used.  Zero is the
ground row.  In no case may \fIPos\fP be less than zero.  
.XP
\fINeg\fP input  (\fBint\fP)
.br
Number of negative node.  Must be in the range of [0..\fISize\fP]
unless either the options EXPANDABLE or TRANSLATE are used. Zero is the
ground row.  In no case may \fINeg\fP be less than zero.  
.XP
\fIEqn\fP input  (\fBint\fP)
.br
Row that contains the branch equation.  Must be in the range of [1..\fISize\fP]
unless either the options EXPANDABLE or TRANSLATE are used. In no case may 
\fIEqn\fP be less than one.  
.XP
\fITemplate\fP  output  (\fBstruct \fIspTemplate\fB *\fR)
.br
Collection of pointers to four elements that are later used to directly 
address elements.  User must supply the template, this routine will 
fill it.
.RE
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
Error is not cleared in this routine.
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
QUAD_ELEMENT
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetOnes()\fP
.XE
.Lb spGetQuad
Similar to \fIspGetAdmittance()\fP, except that \fIspGetAdmittance()\fP
only handles 2-terminal components, whereas \fIspGetQuad()\fP handles
simple 4-terminals as well.  These 4-terminals are simply generalized
2-terminals with the option of having the sense terminals different
from the source and sink terminals.  \fIspGetQuad()\fP installs four
elements into the matrix and returns their pointers in the
\fITemplate\fR structure, which is an output.  The pointers are
arranged in \fITemplate\fR such that when passed to one of the
\fIspADD_QUAD()\fR macros along with an admittance, the admittance will be
added to the elements at [\fIRow1,Col1\fP] and [\fIRow2,Col2\fP] and
subtracted from the elements at [\fIRow1,Col2\fP] and
[\fIRow2,Col1\fP].  The routine works fine if any of the rows and
columns are zero.  This routine is only to be used before
\fIspMNA_Preorder()\fP, \fIspFactor()\fP or \fIspOrderAndFactor()\fP
unless TRANSLATE is set true.
.sp 0.3v
.nf
\fBint\fI spGetQuad( Matrix, Row1, Row2, Col1, Col2, Template )\fR
.XP
\(bu Returned:
.br
The error code is returned.  Possible errors are listed below.
\fIspGetQuad()\fR does not clear the error state, so it is possible to
ignore the return code of each \fIspGetQuad()\fR call, and check for
errors after constructing the whole matrix by calling \fIspError()\fR.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix that quad is to be entered in.
.XP
\fIRow1\fP  input  (\fBint\fP)
.br
First row index for the elements.  \fIRow1\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler flags are
set true.  In either case \fIRow1\fP must not be negative.
.XP
\fIRow2\fP  input  (\fBint\fP)
.br
Second row index for the elements.  \fIRow2\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler flags are
set true.  In either case \fIRow2\fP must not be negative.
.XP
\fICol1\fP  input  (\fBint\fP)
.br
First column index for the elements.  \fICol1\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler flags are
set true.  In either case \fICol1\fP must not be negative.
.XP
\fICol2\fP  input  (\fBint\fP)
.br
Second column index for the elements.  \fICol2\fP must be in the range
[0..\fISize\fP] unless either the TRANSLATE or EXPANDABLE compiler flags are
set true.  In either case \fICol2\fP must not be negative.
.XP
\fITemplate\fP  output  (\fBstruct \fIspTemplate\fB *\fR)
.br
Collection of pointers to four elements that are later used to directly 
address elements.  User must supply the template, this routine will 
fill it.
.RE
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
Error is not cleared in this routine.
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
QUAD_ELEMENT
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetQuad()\fP
.XE
.Lb spGetSize
Returns the size of the matrix, either the internal or external size of
the matrix is returned.  The internal size is the actual number of rows
and columns in the matrix.  The external size is equal to the largest
row or column number.  These numbers will be the same unless the
TRANSLATE option is enabled.
.sp 0.3v
.nf
\fBint\fI spGetSize( Matrix, External )\fR
.XP
\(bu Returned:
.br
The size of the matrix.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix for which the size is desired.
.XP
\fIExternal\fP  input  (\fBint\fP)
.br
If \fIExternal\fR is nonzero, the external size of the matrix is returned,
otherwise the internal size of the matrix is returned.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspGetSize()\fP
.XE
.Lb spInitialize
\fIspInitialize()\fR is a user customizable way to initialize the
matrix.  Passed to this routine is a function pointer.
\fIspInitialize()\fR sweeps through every element in the matrix and
checks the \fIpInitInfo\fR pointer (the user supplied pointer).  If the
\fIpInitInfo\fR is NULL, which is true unless the user changes it
(always true for fill-ins), then the element is zeroed.  Otherwise, the
function pointer is called and passed the \fIpInitInfo\fR pointer as
well as the element pointer and the external row and column numbers
allowing the user to set the value of each element and perhaps the
right-hand side vector.
.LP
The user function (\fIpInit()\fR) is expected to return a nonzero integer
if there is a fatal error and zero otherwise.  Upon encountering a nonzero
return code, \fIspInitialize()\fR terminates and returns the error code.
.LP
The \fISparse\fR error state is cleared to spOKAY in this routine.
.sp 0.3v
.nf
\fBint\fI spInitialize( Matrix, pInit )\fR
.XP
\(bu Returns:
.br
The error code returned by \fIpInit\fR.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix that is to be initialized.
.XP
\fIpInit\fP  input  (\fB(*int)()\fP)
.br
Pointer to a function that, given a pointer to an element, a pointer to
the users data structure containing initialization information for that
element, and the row and column number of the element, initializes it.
.RE

.LP
\fBint\fI pInit( pElement, pInitInfo, Row, Col )\fR
.XP
\(bu Returns:
.br
Nonzero if fatal error, zero otherwise.
.XP
\(bu Arguments:
.RS
.XP
\fIpElement\fP  input  (spREAL \fB*\fP)
.br
The pointer to the real portion of the element.  The real portion can
be accessed using either \fI*pElement\fR or \fIpElement[0]\fR.  The
imaginary portion can be accessed using either \fI*(pElement+1)\fR or
\fIpElement[1]\fR.
.XP
\fIpInitInfo\fP  input  (\fBchar *\fP)
.br
The user-installed pointer to the initialization data structure.
.XP
\fIRow\fP  input  (\fBint\fP)
.br
The external row number of the element.
.XP
\fICol\fP  input  (\fBint\fP)
.br
The external column number of the element.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
INITIALIZE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspInitialize()\fP
.XE
.Lb spInstallInitInfo
With the INITIALIZE compiler option enabled \fISparse\fP allows the
user to keep initialization information with each structurally nonzero
matrix element.  Each element has a pointer (referred to as
\fIpInitInfo\fR) that is set and used by the user.  This routine installs
the pointer \fIpInitInfo\fR into a particular matrix element.
.sp 0.3v
.nf
\fBvoid\fI spInstallInitInfo( pElement, pInitInfo )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIpElement\fR  input  (spREAL \fB*\fR)
.br
Pointer to the element to which \fIpInitInfo\fR is to be attached.
.XP
\fIpInitInfo\fR  input  (\fBchar *\fR)
.br
The pointer \fIpInitInfo\fR.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
INITIALIZE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspInstallInitInfo()\fP
.XE
.Lb spLargestElement
If this routine is called before the matrix is factored, it returns the
absolute value of the largest element in the matrix.  If called after
the matrix has been factored, it returns a lower bound on the absolute
value of the largest element that occurred in any of the reduced
submatrices during the factorization.  The ratio of these two numbers
(factored/unfactored) is the growth, which can be used to determine if
the pivoting order is adequate.  A large growth implies that
considerable error has been made in the factorization and that it is
probably a good idea to reorder the matrix.  If a large growth in
encountered after using \fIspFactor()\fR, reconstruct the matrix and
refactor using \fIspOrderAndFactor()\fR.  If a large growth is
encountered after using \fIspOrderAndFactor()\fR, refactor using
\fIspOrderAndFactor()\fR with the pivot threshold increased, say to
0.1.
.sp 0.3v
.nf
spREAL\fI spLargestElement( Matrix )\fR
.XP
\(bu Returns:
.br
If matrix is unfactored, returns the magnitude of the largest element in
the matrix.  If the matrix is factored, a bound on the magnitude of the
largest element in any of the reduced submatrices is returned.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fR  input  (\fBchar *\fR)
.br
Pointer to the matrix.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
STABILITY
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspLargestElement()\fP
.XE
.Lb spMNA_Preorder
This routine massages modified node admittance matrices to improve the
performance of \fIspOrderAndFactor()\fP.  It tries to remove structural
zeros from the diagonal by exploiting the fact that the row and column
associated with a zero diagonal usually have structural ones placed
symmetrically.  For this routine to work, the structural ones must be
exactly equal to either one or negative one.  This routine should be
used only on modified node admittance matrices and must be executed
after the matrix has been built but before \fIspScale()\fP, \fIspNorm()\fP,
\fIspMultiply()\fP, \fIspFactor()\fP, \fIspOrderAndFactor()\fP or
\fIspDeleteRowAndCol()\fP are executed.  It should be executed for the
initial factorization only.
.sp 0.3v
.nf
\fBvoid\fI spMNA_Preorder( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.IP
Pointer to the matrix to be preordered.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
MODIFIED_NODAL
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspMNA_Preorder()\fP
.XE
.Lb spMultiply
Multiplies \fIMatrix\fP by \fISolution\fP on the right to find \fIRHS\fP.
Assumes matrix has not been factored.  This routine can be  used as a
test to see if solutions are correct.
.sp 0.3v
\fBvoid\fI spMultiply( Matrix, RHS, Solution )\fR
.br
\fBvoid\fI spMultiply( Matrix, RHS, Solution, iRHS, iSolution )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix.
.XP
\fIRHS\fP  output  (spREAL[])
.br
\fIRHS\fP is the right hand side vector.  This is what is being
solved for.  \fIRHS\fP contains only the real portion of the
right-hand side if spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fISolution\fP  input  (spREAL[])
.br
\fISolution\fP is the vector being multiplied by the matrix.
\fISolution\fP contains only the real portion of that vector if
spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiRHS\fP  output  (spREAL[])
.br
\fIiRHS\fP is the imaginary portion of the right hand side. This
is what is being solved for.  It is only necessary to supply
\fIiRHS\fP if the matrix is complex and
spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiSolution\fP  input  (spREAL[])
.br
\fIiSolution\fP is the imaginary portion of the vector being
multiplied by the matrix.  It is only necessary to supply
\fIiRHS\fP if the matrix is complex and
spSEPARATED_COMPLEX_VECTORS is set true.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
MULTIPLICATION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspMultiply()\fP
.XE
.Lb spMultTransposed
Multiplies transposed \fIMatrix\fP by \fISolution\fP on the right to
find \fIRHS\fP.  Assumes matrix has not been factored.  This routine can
be used as a test to see if solutions are correct.
.sp 0.3v
\fBvoid\fI spMultTransposed( Matrix, RHS, Solution )\fR
.br
\fBvoid\fI spMultTransposed( Matrix, RHS, Solution, iRHS, iSolution )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix.
.XP
\fIRHS\fP  output  (spREAL[])
.br
\fIRHS\fP is the right hand side vector.  This is what is being
solved for.  \fIRHS\fP contains only the real portion of the
right-hand side if spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fISolution\fP  input  (spREAL[])
.br
\fISolution\fP is the vector being multiplied by the matrix.
\fISolution\fP contains only the real portion of that vector if
spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiRHS\fP  output  (spREAL[])
.br
\fIiRHS\fP is the imaginary portion of the right hand side. This
is what is being solved for.  It is only necessary to supply
\fIiRHS\fP if the matrix is complex and
spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiSolution\fP  input  (spREAL[])
.br
\fIiSolution\fP is the imaginary portion of the vector being
multiplied by the matrix.  It is only necessary to supply
\fIiRHS\fP if the matrix is complex and
spSEPARATED_COMPLEX_VECTORS is set true.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
MULTIPLICATION
TRANSPOSE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspMultTransposed()\fP
.XE
.Lb spNorm
Computes and returns the L-infinity norm of an unfactored matrix.
This number is used in computing the condition number of the matrix.
It is a fatal error to pass this routine a factored matrix.
.sp 0.3v
spREAL\fI spNorm( Matrix )\fR
.XP
\(bu Returns:
.br
The largest absolute row sum (the L-infinity norm) of the matrix.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
CONDITION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspNorm()\fP
.XE
.Lb spOrderAndFactor
This routine chooses a pivot order for the matrix and factors it into LU
form.  It handles both the initial factorization and subsequent
factorizations when a reordering or threshold pivoting is desired.  This
is handled in a manner that is transparent to the user.
.sp 0.3v
\fBint\fI spOrderAndFactor( Matrix, RHS, Threshold, AbsoluteThreshold, DiagPivoting )\fR
.KE
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspOrderAndFactor()\fP
.XE
.XP
\(bu Returned:
.br
The error code is returned.  Possible errors are listed below.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix to be factored. 
.XP
\fIRHS\fP  input  (spREAL[])
.br
Representative \fIRHS\fP vector that is used to determine pivoting order
when the right-hand side vector is sparse.  If a term in \fIRHS\fP
is zero, it is assumed that it will usually be zero.  Conversely, a
nonzero term in \fIRHS\fP indicates that the term will often be
nonzero.  If \fIRHS\fP is a NULL pointer then the right-hand side vector is
assumed to be full and it is not used when determining the pivoting
order.
.XP
\fIThreshold\fR  input  (spREAL)
.br
This is the pivot threshold, which should be between zero and one.  If
it is one then the pivoting method becomes complete pivoting, which is
very slow and tends to fill up the matrix.  If it is set close to zero
the pivoting method becomes strict Markowitz with no threshold.  The
pivot threshold is used to eliminate pivot candidates that would cause
excessive element growth if they were used.  Element growth is the cause
of roundoff error, which can occur even in well-conditioned matrices.
Setting the threshold large will reduce element growth and roundoff
error, but setting it too large will cause execution time to be
excessive and will result in a large number of fill-ins.  If this
occurs, accuracy can actually be degraded because of the large number of
operations required on the matrix due to the large number of fill-ins.
A good value for diagonal pivoting seems to be 0.001 while a good value
for complete pivoting appears to be 0.1.  The default is chosen by
giving a value larger than one or less than or equal to zero.  Once the
pivot threshold is set, the value becomes the new default for later
calls to \fIspOrderAndFactor\fR.  The threshold value should be
increased and the matrix re-solved if growth is found to be excessive.
Changing the pivot threshold does not improve performance on matrices
where growth is low, as is often the case with ill-conditioned matrices.
The default value of \fIThreshold\fR was choosen for use with nearly
diagonally dominant matrices such as node- and modified-node admittance
matrices.  For these matrices it is usually best to use diagonal
pivoting.  For matrices without a strong diagonal, it is usually best to
use a larger threshold, such as 0.01 or 0.1.
.XP
\fIAbsoluteThreshold\fP  input  (spREAL)
.br
The absolute magnitude an element must have to be considered as a pivot
candidate, except as a last resort.  This number should be set
significantly smaller than the smallest diagonal element that is is
expected to be placed in the matrix.  If there is no reasonable
prediction for the lower bound on these elements, then
\fIAbsoluteThreshold\fP should be set to zero.  \fIAbsoluteThreshold\fP
is used to reduce the possibility of choosing as a pivot an element
that has suffered heavy cancellation and as a result mainly consists of
roundoff error.  Note that if \fIAbsoluteThreshold\fP is set too large,
it could drastically increase the time required to factor and solve
the matrix.  \fIAbsoluteThreshold\fP should be nonnegative.
If no element in the matrix is larger than \fIAbsoluteThreshold\fR, the
warning spILL_CONDITIONED is returned.
.XP
\fIDiagPivoting\fP  input  (int)
.br
A flag indicating that pivot selection should be confined to the
diagonal if possible.  If \fIDiagPivoting\fR is nonzero and if
DIAGONAL_PIVOTING is enabled pivots will be chosen only from the
diagonal unless there are no diagonal elements that satisfy the
threshold criteria.  Otherwise, the entire reduced submatrix is
searched when looking for a pivot.  The diagonal pivoting in
\fISparse\fR is efficient and well refined, while the complete
pivoting is not.  For symmetric and near symmetric matrices, it is best
to use diagonal pivoting because it results in the best performance
when reordering the matrix and when factoring the matrix without
ordering.  If there is a considerable amount of nonsymmetry in the
matrix, then complete pivoting may result in a better equation
ordering simply because there are more pivot candidates to choose
from.  A better ordering results in faster subsequent factorizations.
However, the initial pivot selection process takes considerably longer
for complete pivoting.
.RE
.KS
.XP
\(bu Possible errors:
.nf
spNO_MEMORY
spSINGULAR
spILL_CONDITIONED
.Le
.RE
.Lb spPartition
This routine determines the cost to factor each row using both direct
and indirect addressing and decides, on a row-by-row basis, which
addressing mode is fastest.  This information is used in
\fIspFactor()\fR to speed the factorization.
.LP
When factoring a previously ordered matrix using \fIspFactor()\fR,
f\ISparse\fR operates on a row-at-a-time basis.  For speed, on each
step, the row being updated is copied into a full vector and the
operations are performed on that vector.  This can be done one of two
ways, either using direct addressing or indirect addressing.  Direct
addressing is fastest when the matrix is relatively dense and indirect
addressing is best when the matrix is quite sparse.  The user selects
the type of partition used with \fIMode\fR.  If \fIMode\fR is set to
spDIRECT_PARTITION, then the all rows are placed in the direct
addressing partition.  Similarly, if \fIMode\fR is set to
spINDIRECT_PARTITION, then the all rows are placed in the indirect
addressing partition.  By setting \fIMode\fR to spAUTO_PARTITION, the
user allows \fISparse\fR to select the partition for each row
individually.  \fIspFactor()\fR generally runs faster if \fISparse\fR
is allowed to choose its own partitioning, however choosing a partition
is expensive.  The time required to choose a partition is of the same
order of the cost to factor the matrix.  If you plan to factor a large
number of matrices with the same structure, it is best to let
\fISparse\fR choose the partition.  Otherwise, you should choose the
partition based on the predicted density of the matrix.  By default
(i.e., if \fIspPartition()\fR is never called), \fISparse\fR chooses
the partition for each row individually.
.sp 0.3v
\fBvoid\fI spPartition( Matrix, Mode )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix to be partitioned.
.XP
\fIMode\fP  input  (\fBint\fR)
.br
Mode must be one of three special codes: spDIRECT_PARTITION,
spINDIRECT_PARTITION, or spAUTO_PARTITION.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspPartition()\fP
.XE
.Lb spPrint
Formats and send the matrix to standard output.  Some elementary
statistics are also output.  The matrix is output in a format that is
readable by people.  This routine should not be used on large matrices.
.sp 0.3v
\fBvoid\fI spPrint( Matrix, PrintReordered, Data, Header )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix to be printed.
.XP
\fIPrintReordered\fP  input  (\fBint\fP)
.br
Indicates whether the matrix should be printed out in its original
form, as input by the user, or whether it should be printed in its
reordered form, as used internally by the matrix routines.  A zero
indicates that the matrix should be printed as inputed, a one indicates
that it should be printed reordered.
.XP
\fIData\fP  input  (\fBint\fP)
.br
Boolean flag that when false indicates that output should be compressed
such that only the existence of an element should be indicated rather
than giving the actual value.  Thus 10 times as many elements can be
printed on a row.  A zero indicates that the matrix should be printed
compressed.  A one signifies that the matrix should be printed in all
its glory.
.XP
\fIHeader\fP  input  (\fBint\fP)
.br
A flag indicating that extra information should be printed, such as
row and column numbers.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
DOCUMENTATION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspPrint()\fP
.XE
.Lb spPseudoCondition
Computes the magnitude of the ratio of the largest to the smallest
pivots.  This quantity is an indicator of ill-conditioning in the
matrix.  If this ratio is large, and if the matrix is scaled such that
uncertainties in the right-hand side vector and the matrix entries are
equilibrated, then the matrix is ill-conditioned.  However, a small
ratio does not necessarily imply that the matrix is well-conditioned.
This routine must only be used after a matrix has been factored by
\fIspOrderAndFactor()\fR or \fIspFactor()\fR and before it is cleared
by \fIspClear()\fR or \fIspInitialize()\fR.  The pseudocondition is
faster to compute than the condition number calculated by
\fIspCondition()\fR, but is not as informative.
.sp 0.3v
spREAL \fI spPseudoCondition( Matrix )\fR
.XP
\(bu Returns:
.br
The magnitude of the ratio of the largest to smallest
pivot used during previous factorization.  If the matrix was singular,
zero is returned.
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
PSEUDOCONDITION
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspPseudoCondition()\fP
.XE
.Lb spRoundoff
Returns a bound on the magnitude of the largest element in
$bold E ~=~ bold A - bold LU$, where $bold E$ represents error in the
matrix resulting from roundoff during the factorization.
.sp 0.3v
spREAL \fI spRoundoff( Matrix, Rho )\fR
.XP
\(bu Returns:
.br
Returns a bound on the magnitude of the largest element in
$bold E ~=~ bold A - bold LU$.
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix.  Matrix must be factored.
.XP
\fIRho\fP  input  (spREAL)
.br
The bound on the magnitude of the largest element in any of the reduced
submatrices.  This is the number computed by the function
\fIspLargestElement()\fR when given a factored matrix.  If this number
is negative, the bound will be computed automatically.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
STABILITY
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspRoundoff()\fP
.XE
.Lb spScale
This function scales the matrix to enhance the possibility of finding a
good pivoting order.  Note that scaling enhances accuracy of the
solution only if it affects the pivoting order, so it only makes sense
to scale the matrix before \fIspOrderAndFactor()\fP.  There are several
things to take into account when choosing the scale factors.  First,
the scale factors are directly multiplied times the elements in the
matrix.  To prevent roundoff, each scale factor should be equal to an
integer power of the number base of the machine.  Since most machines
operate in base two, scale factors should be a power of two.  Second,
the matrix should be scaled such that the matrix of element
uncertainties is equilibrated.  Third, this function multiplies the
scale factors times the elements, so if one row tends to have
uncertainties 1000 times smaller than the other rows, then its scale
factor should be 1024, not 1/1024.  Fourth, to save time, this function
does not scale rows or columns if their scale factors are equal to
one.  Thus, the scale factors should be normalized to the most common
scale factor.  Rows and columns should be normalized separately.  For
example, if the size of the matrix is 100 and 10 rows tend to have
uncertainties near 1e-6 and the remaining 90 have uncertainties near
1e-12, then the scale factor for the 10 should be 1/1,048,576 and the
scale factors for the remaining 90 should be 1. Fifth, since this
routine directly operates on the matrix, it is necessary to apply the
scale factors to the \fIRHS\fP and \fISolution\fP vectors.  It may be
easier to simply use \fIspOrderAndFactor()\fP on a scaled matrix to
choose the pivoting order, and then throw away the matrix.  Subsequent
factorizations, performed with \fIspFactor()\fP, will not need to have
the \fIRHS\fP and \fISolution\fP vectors descaled.
.sp 0.3v
\fBvoid\fI spScale( Matrix, RHS_ScaleFactors, SolutionScaleFactors )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to the matrix to be scaled.
.XP
\fIRHS_ScaleFactors\fP  input  (spREAL[])
.br
The array of RHS scale factors.  These factors scale the rows.  All
scale factors are real-valued.
.XP
\fISolutionScaleFactors\fP  input  (spREAL[])
.br
The array of Solution scale factors.  These factors scale the columns.
All scale factors are real-valued.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
SCALING
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspScale()\fP
.XE
.Lb spSetComplex
The type of the matrix may then be toggled back and forth between
complex and real.  This function changes the type of matrix to
complex.  For the matrix to be set complex, the compiler option spCOMPLEX
must be set true.
.sp 0.3v
\fBvoid\fI spSetComplex( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.IP
The matrix that is to be to be complex.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspSetComplex()\fP
.XE
.Lb spSetReal
The type of the matrix may then be toggled back and forth between
complex and real.  This function changes the type of matrix to
real.  For the matrix to be set real, the compiler option REAL
must be set true.
.sp 0.3v
\fBvoid\fI spSetReal( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix that is to be real.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspSetReal()\fP
.XE
.Lb spSolve
Performs the forward and backward elimination to find the unknown 
\fISolution\fP vector from \fIRHS\fP and the factored matrix.  
.sp 0.3v
\fBvoid\fI spSolve( Matrix, RHS, Solution )\fR
.br
\fBvoid\fI spSolve( Matrix, RHS, Solution, iRHS, iSolution )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix. 
.XP
\fIRHS\fP  input  (spREAL[])
.br
\fIRHS\fP is the input data array, the right-hand side vector. 
\fIRHS\fP contains only the real portion of the right-hand side vector if
spSEPARATED_COMPLEX_VECTORS is set true.  \fIRHS\fP is undisturbed and may
be reused for other solves. 
.XP
\fISolution\fP  output  (spREAL[])
.br
\fISolution\fP is the output data array, the unknown vector. This
routine is constructed such that \fIRHS\fP and \fISolution\fP can be
the same array.  \fISolution\fP contains only the real portion of the
unknown vector if spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiRHS\fP  input  (spREAL[])
.br
\fIiRHS\fP is the imaginary portion of the input data array, the
right-hand side vector. This data is undisturbed and may be reused for other
solves.  This argument is unnecessary if the matrix is real or
spSEPARATED_COMPLEX_VECTORS is set false.
.XP
\fIiSolution\fP  output  (spREAL[])
.br
\fIiSolution\fP is the imaginary portion of the output data array.
This routine is constructed such that \fIiRHS\fP and
\fIiSolution\fP can be the same array.  This argument is
unnecessary if the matrix is real or spSEPARATED_COMPLEX_VECTORS is
set false.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspSolve()\fP
.XE
.Lb spSolveTransposed
Performs the forward and backward elimination to find the unknown
\fISolution\fP vector from \fIRHS\fP and the transposed factored
matrix. This routine is useful when performing sensitivity analysis on
a circuit using the adjoint method.
.sp 0.3v
\fBvoid\fI spSolveTransposed( Matrix, RHS, Solution )\fR
.br
\fBvoid\fI spSolveTransposed( Matrix, RHS, Solution, iRHS, iSolution )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix. 
.XP
\fIRHS\fP  input  (spREAL[])
.br
\fIRHS\fP is the input data array, the right-hand side vector.
\fIRHS\fP contains only the real portion of the right-hand side
vector if spSEPARATED_COMPLEX_VECTORS is set true.  \fIRHS\fP is
undisturbed and may be reused for other solves.
.XP
\fISolution\fP  output  (spREAL[])
.br
\fISolution\fP is the output data array, the unknown vector. This
routine is constructed such that \fIRHS\fP and \fISolution\fP can be
the same array.  \fISolution\fP contains only the real portion of the
unknown vector if spSEPARATED_COMPLEX_VECTORS is set true.
.XP
\fIiRHS\fP  input  (spREAL[])
.br
\fIiRHS\fP is the imaginary portion of the input data array, the
right-hand side vector. This data is undisturbed and may be reused for other
solves.  This parameter is unnecessary if the matrix is real or
spSEPARATED_COMPLEX_VECTORS is set false.
.XP
\fIiSolution\fP  output  (spREAL[])
.br
\fIiSolution\fP is the imaginary portion of the output data array.
This routine is constructed such that \fIiRHS\fP and
\fIiSolution\fP can be the same array.  This parameter is
unnecessary if the matrix is real or spSEPARATED_COMPLEX_VECTORS is
set false.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
TRANSPOSE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspSolveTransposed()\fP
.XE
.Lb spStripFills
\fIspStripFills()\fP strips all accumulated fill-ins from a matrix.  This
is often a useful thing to do before reordering a matrix to help insure
that subsequent factorizations will be as efficient as possible.
.sp 0.3v
.nf
\fBvoid\fI spStripFills( Matrix )\fR
.XP
\(bu Argument:
.RS
.XP 
\fIMatrix\fP  input  (\fBchar *\fP)
.br
The matrix to be stripped.
.RE
.XP
\(bu Compiler options that must be set for this routine to exist:
.nf
STRIP
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspStripFills()\fP
.XE
.Lb spWhereSingular
This function returns the row and column number where the matrix was
detected as singular or where a zero pivot was found.
.sp 0.3v
\fBvoid\fI spWhereSingular( Matrix, Row, Col )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIMatrix\fP  input  (\fBchar *\fP)
.br
Pointer to matrix.
.XP
\fIRow\fP  output  (\fBint *\fP)
.br
The row number.
.XP
\fIRow\fP  output  (\fBint *\fP)
.br
The column number.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspWhereSingular()\fP
.XE
.bp
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  MACRO FUNCTIONS"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Macro Functions
.XE
These macro functions are used to quickly enter data into the matrix
using pointers.  These pointers are originally acquired by the user
from \fIspGetElement()\fP, \fIspGetAdmittance()\fP, \fIspGetQuad()\fP,
and \fIspGetOnes()\fP during the initial loading of the
matrix.  These macros work correctly even if the elements they are to add
data to are in row or column zero.
.PP
The macros reside in the file \fBspExports.h\fP.  To use them,
this file must be included in the file of the calling routine and that
routine must be written in C.
.Lb spADD_REAL_ELEMENT
Macro function that adds a real value to an element in the matrix by a pointer. 
.sp 0.3v
\fIspADD_REAL_ELEMENT( pElement , Real )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIpElement\fP  input  (spREAL \fB*\fP)
.br
A pointer to the element to which \fIReal\fR is to be added.
.XP
\fIReal\fP  input  (spREAL)
.br
The real value that is to be added to the element.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_REAL_ELEMENT()\fP
.XE
.Lb spADD_IMAG_ELEMENT
Macro function that adds a imaginary value to an element in the matrix
by a pointer.
.sp 0.3v
\fIspADD_IMAG_ELEMENT( pElement , Imag )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIpElement\fP  input  (spREAL \fB*\fP)
.br
A pointer to the element to which \fIImag\fR is to be added.
.XP
\fIImag\fP  input  (spREAL)
.br
The imaginary value that is to be added to the element.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_IMAG_ELEMENT()\fP
.XE
.Lb spADD_COMPLEX_ELEMENT
Macro function that adds a complex value to an element in the matrix
by a pointer. 
.sp 0.3v
\fIspADD_COMPLEX_ELEMENT( pElement, Real, Imag )\fR
.XP
\(bu Arguments:
.RS
.XP
\fIpElement\fP  input  (spREAL \fB *\fP)
.br
A pointer to the element to which \fIReal\fR and \fIImag\fR are to be added.
.XP
\fIReal\fP  input  (spREAL)
.br
The real value that is to be added to the element.
.XP
\fIImag\fP  input  (spREAL)
.br
The imaginary value that is to be added to the element.
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_COMPLEX_ELEMENT()\fP
.XE
.Lb spADD_REAL_QUAD
Macro that adds a real value to the four elements specified by
\fITemplate\fP.  The value is added to the first two elements in
\fITemplate\fR, and subtracted from the last two.
.sp 0.3v
\fIspADD_REAL_QUAD( Template, Real )\fR
.XP
\(bu Arguments:
.RS
.XP
\fITemplate\fP  input  (\fBstruct \fIspTemplate\fR)
.br
Data structure containing the pointers to four matrix elements.
.XP
\fIReal\fP  input  (spREAL)
.br
Real value to be added to the elements. 
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_REAL_QUAD()\fP
.XE
.Lb spADD_IMAG_QUAD
Macro that adds an imaginary value to the four elements specified by
\fITemplate\fP.  The value is added to the first two elements in
\fITemplate\fR, and subtracted from the last two.
.sp 0.3v
\fIspADD_IMAG_QUAD( Template, Imag )\fR
.XP
\(bu Arguments:
.RS
.XP
\fITemplate\fP  input  (\fBstruct \fIspTemplate\fR)
.br
Data structure containing the pointers to four matrix elements.  
.XP
\fIImag\fP  input  (spREAL)
.br
Imaginary value to be added to the elements. 
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_IMAG_QUAD()\fP
.XE
.Lb spADD_COMPLEX_QUAD
Macro that adds a complex value to the four elements specified by
\fITemplate\fP.  The value is added to the first two elements in
\fITemplate\fR, and subtracted from the last two.
.sp 0.3v
\fIspADD_COMPLEX_QUAD( Template, Real, Imag )\fR
.XP
\(bu Arguments:
.RS
.XP
\fITemplate\fP  input  (\fBstruct \fIspTemplate\fR)
.br
Data structure containing the pointers to four matrix elements.  
.XP
\fIReal\fP  input  (spREAL)
.br
Real value to be added to the elements. 
.XP
\fIImag\fP  input  (spREAL)
.br
Imaginary value to be added to the elements. 
.RE
.Le
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  \fIspADD_COMPLEX_QUAD()\fP
.XE
.bp
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  CONFIGURING SPARSE"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Configuring \fISparse\fR
.XE
.PP
\fISparse\fR has a extensive set of options and parameters that can be
set at compile time to alter the personality of the program.  They also are
used to eliminate routines that are not needed so as to reduce the
amount of memory required to hold the object code.  These options and
parameters consist of macros definitions and are contained in the
file \fBspConfig.h\fR.  To configure \fISparse\fR, \fBspConfig.h\fP must be
edited and then \fISparse\fP must be recompiled.
.PP
Some terminology should be defined.  The Markowitz row count is the number
of non-zero elements in a row excluding the one being considered as pivot.
There is one Markowitz row count for every row.  The Markowitz column count
is defined similarly for columns.  The Markowitz product for an element
is the product of its row and column counts. It is a measure of how much
work would be required on the next step of the factorization if that 
element were chosen to be pivot.  A small Markowitz product is desirable.
For a more detailed explanation, see Kundert [kundert86].
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Sparse Options"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Sparse Options
.XE
.KS
.LP
REAL
.LP
This specifies that the routines are expected to handle real systems
of equations.  The routines can be compiled to handle both real and
complex systems at the same time, but there is a slight speed and 
memory advantage if the routines are complied to handle only real
systems of equations.
.KE
.sp 0.6v
.KS
.LP
spCOMPLEX
.LP
This specifies that the routines will be complied to handle complex
systems of equations.
.KE
.sp 0.6v
.KS
.LP
EXPANDABLE
.LP
Setting this compiler flag true makes the matrix expandable 
before it has been reordered.  If the matrix is expandable,
then if an element is added that would be considered out of bounds
in the current matrix, the size of the matrix is increased to hold
that element.  As a result, the size of the matrix need not be known
before the matrix is built.  The matrix can be allocated with size
zero and expanded.  It is possible to expand the size of a matrix
after it is been reordered if TRANSLATE and EXPANDABLE are both set true.
.KE
.sp 0.6v
.KS
.LP
TRANSLATE
.LP
This option allows the set of external row and column numbers to be
non-packed.  In other words, the row and column numbers need not be
contiguous.  The priced paid for this flexibility is that when
TRANSLATE is set true, the time required to initially build the matrix
will be greater because the external row and column number must be
translated into internal equivalents.  This translation brings about
other benefits though.  First, the \fIspGetElement()\fR,
\fIspGetAdmittance()\fR, \fIspGetQuad()\fR, and \fIspGetOnes()\fR
routines may be used after the matrix has been factored.  Further,
elements, and even rows and columns, may be added to the matrix, and
rows and columns may be deleted from the matrix, after it has been
reordered.  Note that when the set of row and column number is not a
packed set, neither are the \fIRHS\fP and \fISolution\fP vectors.  Thus
the size of these vectors must be at least as large as the external
size, which is the value of the largest given row or column numbers.
.KE
.sp 0.6v
.KS
.LP
INITIALIZE
.LP
Causes the \fIspInitialize()\fR, \fIspGetInitInfo()\fR, and
\fIspInstallInitInfo()\fR routines to be compiled.  These routines
allow the user to store and read one pointer in each nonzero element in
the matrix.  \fIspInitialize()\fR then calls a user specified function
for each structural nonzero in the matrix, and includes this pointer as
well as the external row and column numbers as arguments.  This allows
the user to write custom matrix and right-hand side vector
initialization routines.
.KE
.sp 0.6v
.KS
.LP
DIAGONAL_PIVOTING
.LP
Many matrices, and in particular node- and modified-node admittance
matrices, tend to be nearly symmetric and nearly diagonally dominant.
For these matrices, it is a good idea to select pivots from the
diagonal.  With this option enabled, this is exactly what happens,
though if no satisfactory pivot can be found on the diagonal, an
off-diagonal pivot will be used.  If this option is disabled,
\fISparse\fR does not preferentially search the diagonal.  Because of
this, \fISparse\fR has a wider variety of pivot candidates available,
and so presumably fewer fill-ins will be created.  However, the initial
pivot selection process will take considerably longer.  If working with
node admittance matrices, or other matrices with a strong diagonal, it
is probably best to use DIAGONAL_PIVOTING for two reasons.  First,
accuracy will be better because pivots will be chosen from the large
diagonal elements, thus reducing the chance of growth and hence,
roundoff.  Second, a near optimal ordering will be chosen quickly.  If
the class of matrices you are working with does not have a strong
diagonal, do not use DIAGONAL_PIVOTING, but consider using a larger
threshold.  When DIAGONAL_PIVOTING is turned off, the following options
and constants are not used: MODIFIED_MARKOWITZ, MAX_MARKOWITZ_TIES, and
TIES_MULTIPLIER.
.KE
.sp 0.6v
.KS
.LP
ARRAY_OFFSET
.LP
This determines whether arrays start at an index of zero or one.  This
option is necessitated by the fact that standard C convention dictates
that arrays begin with an index of zero but the standard mathematic
convention states that arrays begin with an index of one.  So if you
prefer to start your arrays with zero, or you're calling \fISparse\fP
from some other programming language, use an ARRAY_OFFSET of 0.
Otherwise, use an ARRAY_OFFSET of 1.  Note that if you use an offset of
one, the arrays that you pass to \fISparse\fR must have an allocated
length of one plus the external size of the matrix.  ARRAY_OFFSET must
be either 0 or 1, no other offsets are valid.
.KE
.sp 0.6v
.KS
.LP
spSEPARATED_COMPLEX_VECTORS
.LP
This specifies the format for complex vectors.  If this is set false
then a complex vector is made up of one double sized array of spREALs
in which the real and imaginary numbers are placed alternately in the
array.  In other words, the first entry would be \fIComplex[1].Real\fP,
then comes \fIComplex[1].Imag\fP, then \fIComplex[2].Real\fP, etc.  If
spSEPARATED_COMPLEX_VECTORS is set true, then each complex vector is
represented by two arrays of spREALs, one with the real terms, the
other with the imaginary.
.KE
.sp 0.6v
.KS
.LP
MODIFIED_MARKOWITZ
.LP
This specifies that the modified Markowitz method of pivot
selection is to be used.  The modified Markowitz method differs
from standard Markowitz in two ways.  First, under modified
Markowitz, the search for a pivot can be terminated early if a
adequate (in terms of sparsity) pivot candidate is found.
Thus, when using modified Markowitz, the initial factorization
can be faster, but at the expense of a suboptimal pivoting
order that may slow subsequent factorizations.  The second
difference is in the way modified Markowitz breaks Markowitz
ties.  When two or more elements are pivot candidates and they
all have the same Markowitz product, then the tie is broken by
choosing the element that is best numerically.  The numerically
best element is the one with the largest ratio of its magnitude
to the magnitude of the largest element in the same column,
excluding itself.  The modified Markowitz method results in
marginally better accuracy.
.KE
.sp 0.6v
.KS
.LP
DELETE
.LP
This specifies that the \fIspDeleteRowAndCol()\fP routine should
be compiled.  Note that for this routine to be compiled, both DELETE
and TRANSLATE should be set true.
.KE
.sp 0.6v
.KS
.LP
STRIP
.LP
This specifies that the \fIspStripFills()\fP routine should be compiled.
.KE
.sp 0.6v
.KS
.LP
MODIFIED_NODAL
.LP
This specifies that the \fIspMNA_Preorder()\fR, the routine that
preorders modified node admittance matrices, should be compiled.
This routine results in greater speed and accuracy if used with this
type of matrix.
.KE
.sp 0.6v
.KS
.LP
QUAD_ELEMENT
.LP
This specifies that the routines that allow four related elements to be
entered into the matrix at once should be compiled.  The routines affected by
QUAD_ELEMENT are \fIspGetAdmittance()\fP, \fIspGetQuad()\fP, and
\fIspGetOnes()\fP.
.KE
.sp 0.6v
.KS
.LP
TRANSPOSE
.LP
This specifies that \fIspSolveTranspose()\fR and perhaps
\fIspMultTransposed()\fR, which operate on
the matrix as if it was transposed, should be compiled.
.KE
.sp 0.6v
.KS
.LP
SCALING
.LP
This specifies that the routine that performs scaling on the matrix
should be complied.  Scaling is not strongly supported.  The routine to
scale the matrix is provided, but no routines are provided to scale and
descale the \fIRHS\fP and \fISolution\fP vectors.  It is suggested that
if scaling is desired, it only be performed when the pivot order is
being chosen, which is done in \fIspOrderAndFactor()\fP.  This, and
when the condition number of the matrix is calculated with
\fIspCondition()\fR, are the only times scaling has an effect.  The
scaling may then either be removed from the solution by the user or the
scaled factored matrix may simply be thrown away.
.KE
.sp 0.6v
.KS
.LP
DOCUMENTATION
.LP
This specifies that routines that are used to document the matrix,
\fIspPrint()\fP, \fIspFileMatrix()\fP, \fIspFileVector()\fP, and
\fIspFileStats()\fP, should be compiled.
.KE
.sp 0.6v
.KS
.LP
DETERMINANT
.LP
This specifies that the \fIspDeterminant()\fR routine should be complied.
.KE
.sp 0.6v
.KS
.LP
STABILITY
.LP
This specifies that \fIspLargestElement()\fR and \fIspRoundoff()\fR
should be compiled.  These routines are used to check the stability
(and hence the quality of the pivoting) of the factorization by
computing a bound on the size of the element is the matrix $bold E ~=~
bold A - bold LU$.  If this bound is very high after applying
\fIspOrderAndFactor()\fR, then the pivot threshold should be raised.
If the bound increases greatly after using \fIspFactor()\fR, then the
matrix should probably be reordered.
.KE
.sp 0.6v
.KS
.LP
CONDITION
.LP
This specifies that \fIspCondition()\fR and \fIspNorm()\fR, the code
that computes a good estimate of the condition number of the matrix,
should be compiled.
.KE
.sp 0.6v
.KS
.LP
PSEUDOCONDITION
.LP
This specifies that \fIspPseudoCondition()\fR, the code that computes
a crude and easily fooled indicator of the ill-conditioning in the matrix,
should be compiled.
.KE
.sp 0.6v
.KS
.LP
MULTIPLICATION
.LP
This specifies that \fIspMultiply()\fR and perhaps \fIspMultTransposed()\fR,
the routines that multiply an unfactored matrix by a vector, should be compiled.
.KE
.sp 0.6v
.KS
.LP
FORTRAN
.LP
This specifies that the FORTRAN interface to Sparse1.3 should be
compiled.  The ARRAY_OFFSET option should be set to NO when interfacing
to FORTRAN programs.
.KE
.sp 0.6v
.KS
.LP
DEBUG
.LP
This specifies that additional error checking should be compiled.
The type of errors checked are those that are common when the
matrix routines are first integrated into a user's program.  Once
the routines have been integrated in and are running smoothly, this
option should be turned off.  With DEBUG enabled, \fISparse\fP is
very defensive.  If a \fISparse\fP routine is called improperly,
a message will be printed describing the file and line number where
the error was found and execution is aborted.  One thing that \fISparse\fR
is particularly picky about is calling certain functions after an error
has occurred.  If an error has occurred, do not call \fIspMNA_Preorder()\fR,
\fIspScale()\fR, \fIspOrderAndFactor()\fR, \fIspFactor()\fR,
\fIspSolve()\fR, or \fIspSolveTransposed()\fR until the error has been
cleared by \fIspClear()\fR or \fIspInitialize()\fR.
.KE
.sp 0.6v
.KS
.LP
spCOMPATIBILITY
.LP
This specifies that Sparse1.3 should be configured to be upward
compatible from Sparse1.2.  This option is not suggested for use
in new software.  Sparse1.3, when configured to be compatible with
Sparse1.2, is not completely compatible.  In particular, if
recompiling the calling program, it is necessary to change the names
of the \fISparse1.2\fR include files.  This option will go away on any future
version of \fISparse\fR.
.KE
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Sparse Constants"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Sparse Constants
.XE
.PP
These constants are used throughout the sparse matrix routines.  They
should be set to suit the type of matrices being solved.
.sp 0.6v
.KS
.LP
DEFAULT_THRESHOLD
.LP
The threshold used if the user enters an invalid threshold.  Also the
threshold used by \fIspFactor()\fP when calling
\fIspOrderAndFactor()\fP.  The default threshold should not be less
than or equal to zero nor larger than one.
.KE
.sp 0.6v
.KS
.LP
DIAG_PIVOTING_AS_DEFAULT
.LP
This indicates whether \fIspOrderAndFactor()\fP should use diagonal
pivoting as default.  This issue only arises when \fIspOrderAndFactor()\fP
is called from \fIspFactor()\fP.
.KE
.sp 0.6v
.KS
.LP
SPACE_FOR_ELEMENTS
.LP
This number multiplied by the size of the matrix equals the number
of elements for which memory is initially allocated in \fIspCreate()\fP.
.KE
.sp 0.6v
.KS
.LP
SPACE_FOR_FILL_INS
.LP
This number multiplied by the size of the matrix equals the number
of elements for which memory is initially allocated and specifically
reserved for fill-ins in \fIspCreate()\fP.
.KE
.sp 0.6v
.KS
.LP
ELEMENTS_PER_ALLOCATION
.LP
The number of matrix elements requested from the \fImalloc\fP utility
on each call to it.  Setting this value greater than one reduces the
amount of overhead spent in this system call.
.KE
.sp 0.6v
.KS
.LP
MINIMUM_ALLOCATED_SIZE
.LP
The minimum allocated size of a matrix.  Note that this does not
limit the minimum size of a matrix.  This just prevents having to
resize a matrix many times if the matrix is expandable, large and
allocated with an estimated size of zero.  This number must not
be less than one.
.KE
.sp 0.6v
.KS
.LP
EXPANSION_FACTOR
.LP
The minimum increase in the allocated size of the matrix when it is
expanded.  This number must be greater than one but shouldn't be much
larger than two.
.KE
.sp 0.6v
.KS
.LP
MAX_MARKOWITZ_TIES
.LP
This number is used for two slightly different things, both of which
relate to the search for the best pivot.  First, it is the maximum 
number of elements that are Markowitz tied that will be sifted 
through when trying to find the one that is numerically the best.
Second, it creates an upper bound on how large a Markowitz product
can be before it eliminates the possibility of early termination
of the pivot search.  In other words, if the product of the smallest
Markowitz product yet found and TIES_MULTIPLIER is greater than
MAX_MARKOWITZ_TIES, then no early termination takes place.  Set 
MAX_MARKOWITZ_TIES to some small value if no early termination of 
the pivot search is desired. An array of spREALs is allocated of size 
MAX_MARKOWITZ_TIES so it must be positive and shouldn't be too large.
.KE
.sp 0.6v
.KS
.LP
TIES_MULTIPLIER
.LP
Specifies the number of Markowitz ties that are allowed to occur before
the search for the pivot is terminated early.  Set to some large value
if no early termination of the pivot search is desired.  This number is
multiplied by the Markowitz product to determine how many ties are
required for early termination.  This means that more elements will be
searched before early termination if a large number of fill-ins could
be created by accepting what is currently considered the best choice
for the pivot.  Setting this number to zero effectively eliminates all
pivoting, which should be avoided.  This number must be positive.
.KE
.sp 0.6v
.KS
.LP
DEFAULT_PARTITION
.LP
Which partition mode is used by \fIspPartition()\fR as default.
Possibilities include:
.RS
.XP
spDIRECT_PARTITION  \(em each row used direct addressing, best for
a few relatively dense matrices.
.XP
spINDIRECT_PARTITION  \(em each row used indirect addressing, best
for a few very sparse matrices.
.XP
spAUTO_PARTITION  \(em direct or indirect addressing is chosen on
a row-by-row basis, carries a large overhead, but speeds up
both dense and sparse matrices, best if there is a large
number of matrices that can use the same ordering.
.RE
.KE
.sp 0.6v
.KS
.LP
PRINTER_WIDTH
.LP
Gives the number of characters printable in one page width.
Set to 80 for terminals and 132 for line printers.
.KE
.sp 2
.ne 1.5i
.RT
.LG
.B "\n(sc.\n+(ss:  Machine Constants"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Machine Constants
.XE
.LP
These numbers 
.UL must 
be updated when the program is ported to a new machine.
.sp 0.6v
.KS
.LP
MACHINE_RESOLUTION
.LP
This is the smallest positive real double precision
number $epsilon$ such that $ 1 ~ + ~ epsilon ~ != ~ 1 $.
.KE
.sp 0.6v
.KS
.LP
LARGEST_REAL
.LP
The largest positive real number representable by a \fBdouble\fR.
.KE
.sp 0.6v
.KS
.LP
SMALLEST_REAL
.LP
The smallest positive real number representable by a \fBdouble\fR.
.KE
.sp 0.6v
.KS
.LP
LARGEST_SHORT_INTEGER
.LP
The largest positive integer representable by a \fBshort\fR.
.KE
.sp 0.6v
.KS
.LP
LARGEST_LONG_INTEGER
.LP
The largest positive integer representable by a \fBlong\fR.
.KE
.bp
.LP
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  EXPORTS"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  Exports
.XE

.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Error Codes"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Error Codes
.XE
.PP
Errors are indicated with a integer error code.  Macros definitions for
these error codes are set up and placed in the file \fBspMatrix.h\fP.
They may be imported into the users program to give readable
names to the possible matrix errors.  The possible error codes and
there corresponding macros are:
.sp 0.25i
.KS
.LP
spOKAY  \(em  0
.LP
No error has occurred.
.KE
.KS
.LP
spSMALL_PIVOT  \(em  1
.LP
When reordering the matrix, no element was found which satisfies the
absolute threshold criteria.  The largest element in the matrix was chosen
as pivot.  Nonfatal.
.KE
.KS
.LP
spZERO_DIAG  \(em  2
.LP
Fatal error.  A zero was encountered on the diagonal of the matrix.  This
does not necessarily imply that the matrix is singular.  When this error
occurs, the matrix should be reconstructed and factored using
\fIspOrderAndFactor()\fR.
.KE
.KS
.LP
spSINGULAR  \(em  3
.LP
Fatal error.  Matrix is singular, so no unique solution exists.
.KE
.KS
.LP
spNO_MEMORY  \(em  4
.LP
Fatal error.  Indicates that not enough memory is available from the system
to handle the matrix.
.KE
.KS
.LP
spPANIC  \(em  5
.LP
Fatal error indicating that the routines are being asked to do
something nonsensical or something they are not prepared for.  This
error may occur when the matrix is specified to be real and the
routines are not compiled for real matrices, or when the matrix is
specified to be complex and the routines are not compiled to handle
complex matrices.
.KE
.KS
.LP
spFATAL  \(em  2
.LP
Not an error flag, but rather the dividing line between fatal errors
and warnings.
.KE
.sp 2
.RT
.ne 1i
.LG
.B "\n(sc.\n+(ss:  Data Structures"
.NL
.XS \n%         \" Save section title for table of contents
	\n(sc.\n(ss:  Data Structures
.XE
.PP
There is only one data structure that may need to be imported from
\fISparse\fR by the user.  This data structure is used to hold pointers
to four related elements in matrix.  It is used in conjunction with the
routines
.nf
	\fIspGetAdmittance()\fP
	\fIspGetOnes()\fP
	\fIspGetQuad()\fP
.LP
\fIspGetAdmittance()\fP, \fIspGetOnes()\fP, and \fIspGetQuad()\fP stuff
the structure which is later used by the \fIspADD_QUAD()\fP macros.  It
is also possible for the user to collect four pointers returned by
\fIspGetElement()\fP and stuff them into the template.  The
\fIspADD_QUAD()\fP macros add a value into \fIElement1\fP and
\fIElement2\fP and subtract the value from \fIElement3\fP and
\fIElement4\fP.  The structure is:
.KS
.sp 0.25i
.nf 
\fBstruct\fP \fIspTemplate\fP
{	spREAL    *\fIElement1\fP;
	spREAL    *\fIElement2\fP;
	spREAL    *\fIElement3Negated\fP;
	spREAL    *\fIElement4Negated\fP;
};
.KE
.bp
.LP
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  FORTRAN COMPATIBILITY"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  FORTRAN Compatibility
.XE
.PP
The \fISparse1.3\fR package contains routines that interface to a
calling program written in FORTRAN.  Almost every externally available
\fISparse1.3\fR routine has a counterpart defined with the same name 
except that the `sp' prefix is changed to `sf'.  The
\fIspADD_ELEMENT()\fR and \fIspADD_QUAD()\fR macros are also replaced
with the \fIsfAdd1()\fR and \fIsfAdd4()\fR functions.
.PP
Any interface between two languages is going to have portibility
problems, this one is no exception.  To ease porting the FORTRAN
interface file to different operating systems, the names of the
interface functions can be easily redefined (search for `Routine
Renaming' in \fBspFortran.c\fR).  When interfacing to a FORTRAN program,
the FORTRAN option should be set to YES and the ARRAY_OFFSET option
should be set to NO (see \fBspConfig.h\fR).  For details on the return
value and argument list of a particular interface routine, see the
file \fBspFortran.c\fR.
.PP
A simple example of a FORTRAN program that calls \fISparse\fR follows.

.KS
.B Example:
.RS
.nf
      integer matrix, error, sfCreate, sfGetElement, spFactor
      integer element(10)
      double precision rhs(4), solution(4)
c
c create matrix
      matrix = sfCreate(4,0,error)
c
c reserve elements
      element(1) = sfGetElement(matrix,1,1)
      element(2) = sfGetElement(matrix,1,2)
      element(3) = sfGetElement(matrix,2,1)
      element(4) = sfGetElement(matrix,2,2)
      element(5) = sfGetElement(matrix,2,3)
      element(6) = sfGetElement(matrix,3,2)
      element(7) = sfGetElement(matrix,3,3)
      element(8) = sfGetElement(matrix,3,4)
      element(9) = sfGetElement(matrix,4,3)
      element(10) = sfGetElement(matrix,4,4)
c
c clear matrix
      call sfClear(matrix)
c
c load matrix
      call sfAdd1Real(element(1), 2d0)
      call sfAdd1Real(element(2), \-1d0)
      call sfAdd1Real(element(3), \-1d0)
      call sfAdd1Real(element(4), 3d0)
      call sfAdd1Real(element(5), \-1d0)
      call sfAdd1Real(element(6), \-1d0)
      call sfAdd1Real(element(7), 3d0)
      call sfAdd1Real(element(8), \-1d0)
      call sfAdd1Real(element(9), \-1d0)
      call sfAdd1Real(element(10), 3d0)
      call sfPrint(matrix, .false., .false., .true.)
      rhs(1) = 34d0
      rhs(2) = 0d0
      rhs(3) = 0d0
      rhs(4) = 0d0
c
c factor matrix
      error = sfFactor(matrix)
c
c solve matrix
      call sfSolve(matrix, rhs, solution)
      write (6, 10) solution(1), solution(2), solution(3), solution(4)
   10 format (f 10.2)
      end
.fi
.KE
.RE
.bp
.LP
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  SPARSE TEST PROGRAM"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  \fISparse\fR Test Program
.XE
.PP
The \fISparse\fR package includes a test program that is able to read
matrix equations from text files and print their solution along with
matrix statistics and timing information.  The program can also generate
files containing stripped versions of the unfactored and factored matrix
suitable for plotting using standard plotting programs, such as the
\s-1UNIX\s+1 graph and plot commands.
.LP
The \fISparse\fR test program is invoked using the following syntax.
.ne 2.25i
.IP
\fBsparse\fR [\fIoptions\fR] [\fIfile1\fR] [\fIfile2\fR] ...
.TS
l
l lw(4.5i).
Options:
\fB\-s\fR	T{
Print solution only.
T}
\fB\-r\fI x\fR	T{
Use \fIx\fR as relative threshold.
T}
\fB\-a\fI x\fR	T{
Use \fIx\fR as absolute threshold.
T}
\fB\-n\fI n\fR	T{
Print first \fIn\fR terms of solution vector.
T}
\fB\-i\fI n\fR	T{
Repeat build/factor/solve \fIn\fR times for better timing results.
T}
\fB\-b\fI n\fR	T{
Use column \fIn\fR of matrix as right-hand side vector.
T}
\fB\-p\fR	T{
Create plot files ``\fIfilename\fB.bef\fR'' and ``\fIfilename\fB.aft\fR''.
T}
\fB\-c\fR	T{
Use complete (as opposed to diagonal) pivoting.
T}
\fB\-x\fR	T{
Treat real matrix as complex with imaginary part zero.
T}
\fB\-t\fR	T{
Solve transposed system.
T}
\fB\-u\fR	T{
Print usage message.
T}
.TE
.LP
The presence of certain options is dependent on whether the appropriate
\fISparse\fR option has been enabled.
.LP
If no input files are specified, \fBsparse\fP reads from the standard
input.  The syntax of the input file is as follows.  The matrix begins
with one line of arbitrary text that acts as the label, followed by a
line with the integer size of the matrix and either the \fBreal\fR or
\fBcomplex\fR keywords.  After the header is an arbitrary number of
lines that describe the structural nonzeros in the matrix.  These lines
have the form \fIrow column data\fP, where \fIrow\fR and \fIcolumn\fR
are integers and \fIdata\fP is either one real number for real matrices
or a real/imaginary pair of numbers for complex matrices.  Only one
structural nonzero is described per line and the section ends when
either \fIrow\fR or \fIcolumn\fR are zero.  Following the matrix, an
optional right-hand side vector can be described.  The vector is given
one element per line, the number of element must equal the size of the
matrix.  Only one matrix and one vector are allowed per file, and the
vector, if given, must follow the matrix.

.KS
.B Example:
.RS
.nf
mat0  \(em  Simple matrix.
4	real
1	1	2.0
1	2	-1.0
2	1	-1.0
2	2	3.0
2	3	-1.0
3	2	-1.0
3	3	3.0
3	4	-1.0
4	3	-1.0
4	4	3.0
0	0	0.0
34.0
0.0
0.0
0.0
.fi
.KE
.RE
.bp
.LP
.nr ss 0 1	\" Reset subsection
.RT
.ne 1.5i
.LG
.B "\n+(sc:  SPARSE FILES"
.NL
.XS \n%         \" Save section title for table of contents
\n(sc:  \fISparse\fR Files
.XE
.PP
The following is a list of the files contained in the \fISparse\fR
package and a brief description of their contents.  Of the files, only
\fBspConfig.h\fP is expected to be modified by the user and only
\fBspMatrix.h\fP need be imported into the program that calls \fISparse\fR.

.KS
.LP
\fBspAlloc.c\fR
.LP
This file contains the routines for allocating and deallocating objects
associated with the matrices, including the matrices themselves.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spCreate()
spDestroy()
spError()
spWhereSingular()
spGetSize()
spSetReal()
spSetComplex()
spFillinCount()
spElementCount()
.R
.KE

.KS
.LP
\fBspBuild.c\fR
.LP
This file contains the routines for clearing and loading the matrix.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spClear()
spGetAdmittance()
spGetElement()
spGetInitInfo()
spGetOnes()
spGetQuad()
spInitialize()
spInstallInitInfo()
.R
.KE

.KS
.LP
\fBspCompat.c\fR
.LP
This file contains the routines for making \fISparse1.3\fR upward compatible
from \fISparse1.2\fP.  These routines are not suggested for use in new
software.  These routines will not be available in future versions of
\fISparse\fR.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
AddAdmittanceToMatrix()
AddComplexElementToMatrix()
AddComplexQuadElementToMatrix()
AddElementToMatrix()
AddImagElementToMatrix()
AddImagQuadElementToMatrix()
AddOnesToMatrix()
AddQuadToMatrix()
AddRealElementToMatrix()
AddRealQuadElementToMatrix()
CleanMatrix()
ClearMatrix()
ClearMatrixError()
CreateMatrix()
DecomposeMatrix()
DeleteRowAndColFromMatrix()
DestroyMatrix()
Determinant()
GetMatrixSize()
MatrixElementCount()
MatrixError()
MatrixFillinCount()
MatrixRoundoffError()
MultiplyMatrix()
OrderAndDecomposeMatrix()
OutputMatrixToFile()
PreorderForModifiedNodal()
PrintMatrix()
ScaleMatrix()
SetMatrixComplex()
SetMatrixReal()
SolveMatrix()
SolveTransposedMatrix()
.R
.KE

.KS
.LP
\fBspConfig.h\fR
.LP
This file contains the options that are used to customize the package.
For example, it is possible to specify whether only real or complex
systems of equations are to be solved.  Also included in this file are
the various constants used by the \fISparse\fR package, such as the
amount of memory initially allocated for each matrix and the largest
real number represented by the machine.  The user is expected to modify
this file to maximize the performance of the routines with his/her
matrices.
.KE

.KS
.LP
\fBspDefs.h\fR
.LP
This module contains common data structure definitions and macros for
the sparse matrix routines.  These definitions are meant to remain
hidden from the program that calls the sparse matrix routines.
.KE

.KS
.LP
\fBspDoc\fR
.LP
This reference manual.  \fBspDoc\fP contains the manual in a form
that is readable on-line and \fBspDoc.ms\fP contains the manual in a
form that is suitable for input into the text formatting program
\fItroff\fP using the -ms macros.
.KE

.KS
.LP
\fBspFactor.c\fR
.LP
This file contains the routines for factoring matrices into LU form.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spFactor()
spOrderAndFactor()
spPartition()
.R
.KE

.KS
.LP
\fBspFortran.c\fR
.LP
This file contains the routines for interfacing \fISparse1.3\fR to a
program written in FORTRAN.  The function and argument lists of the
routines in this file are almost identical to their C equivalents
except that they are suitable for calling from a FORTRAN program.
The names of these routines use the `sf' prefix to distinguish them
from their C counterparts.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
sfAdd1Complex()
sfAdd1Imag()
sfAdd1Real()
sfAdd4Complex()
sfAdd4Imag()
sfAdd4Real()
sfClear()
sfCondition()
sfCreate()
sfDeleteRowAndCol()
sfDestroy()
sfDeterminant()
sfElementCount()
sfError()
sfFactor()
sfFileMatrix()
sfFileStats()
sfFileVector()
sfFillinCount()
sfGetAdmittance()
sfGetElement()
sfGetOnes()
sfGetQuad()
sfGetSize()
sfLargestElement()
sfMNA_Preorder()
sfMultTransposed()
sfMultiply()
sfNorm()
sfOrderAndFactor()
sfPartition()
sfPrint()
sfPseudoCondition()
sfRoundoff()
sfScale()
sfSetComplex()
sfSetReal()
sfSolve()
sfSolveTransposed()
sfStripFills()
sfWhereSingular()
.R
.KE

.KS
.LP
\fBspMatrix.h\fR
.NL
.LP
This file contains definitions that are useful to the calling program.
In particular, this file contains error keyword definitions, some macro
functions that are used to quickly enter data into the matrix, the
definition of a data structure that acts as a template for entering
admittances into the matrix, and the type declarations of the various
\fISparse\fR functions.
.KE

.KS
.LP
\fBspOutput.c\fR
.LP
This file contains the output-to-file and output-to-screen routines for
the matrix package.  They are capable of outputting the matrix in either a
form readable by people or a form readable by the \fISparse\fR test program.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spFileMatrix()
spFileStats()
spFileVector()
spPrint()
.R
.KE

.KS
.LP
\fBspRevision\fR
.LP
The history of updates for the program.  This file also includes ordering
information for the \fISparse\fR package.
.KE

.KS
.LP
\fBspSolve.c\fR
.LP
This module contains the forward and backward substitution routines.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spSolve()
spSolveTransposed()
.R
.KE

.KS
.LP
\fBspTest.c\fR
.LP
This module contains a test program for the sparse matrix routines.
It is able to read matrices from files and solve them.  Because of
the large number of options and capabilities built into \fISparse\fR,
it is impossible to have one test routine thoroughly exercise
\fISparse\fR.  Thus, emphasis is on exercising as many capabilities as
is reasonable while also providing a useful tool.
.KE

.KS
.LP
\fBspUtil.c\fR
.LP
This module contains various optional utility routines.
.XP
\(bu User accessible functions contained in this module:
.nf
.I
spCondition()
spDeleteRowAndCol()
spDeterminant()
spLargestElement()
spMNA_Preorder()
spMultiply()
spMultTransposed()
spNorm()
spPseudoCondition()
spRoundoff()
spScale()
spStripFills()
.R
.KE

.KS
.LP
\fBMakefile\fR  
.LP
This file is used in conjunction with the UNIX program \fBmake\fR to
compile the matrix routines and their test program. 
.KE

.KS
.LP
.B make.com  
.LP
This file is used to automatically compile \fISparse\fR under the VMS
operating system.  It needs to modified slightly before being used, see
the installation notes.
.KE
.bp
.LP
.LG
.B "REFERENCES"
.NL
.XS \n%
References
.XE
.KS
.IP [duff86] 15
I. S. Duff, A. M. Erisman, J. K. Reid.
\fIDirect Methods for Sparse Matrices\fP.
Oxford University Press, 1986.
.KE
.KS
.IP [golub86] 15
G. H. Golub, C. F. V. Van Loan.
\fIMatrix Computations\fP.
The Johns Hopkins University Press, 1983.
.KE
.KS
.IP [kundert86] 15
Kenneth S. Kundert.
Sparse matrix techniques.
In \fICircuit Analysis, Simulation and Design\fP, Albert Ruehli
(editor).  North-Holland, 1986.
.KE
.KS
.IP [strang80] 15
Gilbert Strang.
\fILinear Algebra and Its Applications\fP.
Academic Press, 1980.
.KE
.sp 2
.LG
.B Acknowledgements
.NL
.ad
.PP
We would like to acknowledge and thank the those people that
contributed ideas that were incorporated into \fISparse\fP.  In
particular, Jacob White, Kartikeya Mayaram, Don Webber,
Tom Quarles, Howard Ko and Beresford Parlett.

.ds CH   \" stop printing page number
.bp
.PX

