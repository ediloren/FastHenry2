.LP
.po 1.5i
.nr PO 1.5i
.rs
.ls 1
.EQ
delim @@
.EN

.LG
.LG
.I
.ce 2
Sparse1.3
.sp .1i
A Sparse Linear Equation Solver 
.NL
.NL
.sp 0.25i
.LG
.R
.ce 2
.I Kenneth\ S.\ Kundert
.I Alberto\ Sangiovanni-Vincentelli
.NL
.sp 0.25i
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
.PP
The \fISparse1.3\fR package is currently available from the Department
of Electrical Engineering and Computer Sciences of the University of
California, Berkeley.  It was written in the C programming language by
Kenneth Kundert and versions exist for the \s-1UNIX\s+1 and
\s-1VAX/VMS\s+1 operating systems.  Be sure to specify the version when
ordering.  \fISparse1.3\fR has replaced \fISparse1.2\fR; providing
greater capability and speed, and a more refined interface to the
calling program.  Besides the native \fISparse1.3\fR interface, two
others are also provided.  The first is a \fISparse1.2\fR interface,
to make \fISparse1.3\fR backward compatible.  The second allows
\fISparse.13\fR to be called from FORTRAN programs.
.PP
\fISparse\fR is available for a $150.00 charge.  The package includes
the source code on tape, the user's guide, and a large selection of
test matrices.  To obtain a copy of \fISparse\fR, send a check or money
order payable to the \fIRegents of the University of California\fR to:
.sp 0.25
.nf
	EECS Industrial Liaison Program
	461 Cory Hall
	University of California
	Berkeley, CA 94720
.LP
Please allow four weeks for delivery.
