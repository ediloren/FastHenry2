
.LP
.EQ
delim $$
.EN
.if t .po 1.5i
.if t .nr PO 1.5i
.if n .ll 7.5i
.if n .nr LL 7.5i
.rs
.ls 1
.LP
.LG
.LG
.ce
.I "Sparse1.3 Installation Notes"
.NL
.NL
.sp .25i
.LP
.LG
.B "General Description"
.LP
When \fISparse1.3\fR is loaded from tape, a directory structure will be
created.  The head directory is \fBsparse\fR, it contains source and
various utility files.  The \fBdoc\fR subdirectory contains the users
manual and the \fBmatrices\fR subdirectory contains a collection of
test matrices and their solution.
.LP 
Before compiling, the file \fBspConfig.h\fR should be modified to suit
your needs.  Initially, I recommend turning on the DEBUG option and
testing the error state of \fISparse\fR after each call to a
\fISparse\fR routine.  After things are working smoothly you can remove
most of the error checking.  Also in \fBspConfig.h\fR is the machine
constants.  These should be checked and modified before compiling if
you are not using a \s-1VAX\s+1.
.LP
\fISparse\fR is compiled by using the makefile provided.  This is a
script file that automatically compiles each file and creates an object
file for all of \fISparse\fR.  Also created is an executable test
program called \fBsparse\fR.
.LP
The \fISparse1.3\fR manual is contained in \fBspDoc\fR and
\fBspDoc.ms\fR.  \fBspDoc\fR is preformatted and readable on line.
\fBspDoc.ms\fR is a raw \fItroff\fR file.
.sp .25i
.LP
.B UNIX
.LP
The
\fItar\fR
program was used to write
\fISparse\fR 
on the tape at 1600 bpi.  To read it, simply type:
.IP
.B "tar x"
.LP
This will load \fISparse\fR into the current directory and
automatically create the proper subdirectories.  To compile, move into
the \fBsparse\fR directory and type:
.IP
.B make
.LP
To create a printed manual, move into the \fBdoc\fR subdirectory and type:
.IP
.B "make"
.LP
You may have to modify the make file to indicate the name of the printer
the document should be printed on.
.sp .25i
.LP
.B
VMS
.LP
\fISparse\fR was written on the tape in the ansi format at 1600 bpi and
so is readable by \fIcopy\fR.  Once the tape has been read, a series of
files with nonsensical names will exist plus the file
\fBUNPACK.COM\fR.  This file, when executed, creates the directory
hierarchy and moves the \fISparse\fR files into the proper directory
and gives them the proper names.  Thus, \fISparse\fR is installed with
the following command sequence:
.IP
.B "MOUNT MFA0: SPARSE"
.br
.B "COPY MFA0:*.* *.*.*"
.br
.B "@UNPACK"
.LP
The file \fBUNPACK.COM\fR is then no longer needed.
.LP
Before compiling, it is necessary to edit the file \fBMAKE.COM\fR.  The
last line of the file must be changed so that the correct directory
name for \fISparse\fR is given.  This last line must be executed each
time you relogin and plan to use the \fISparse\fR test program; you may
want to copy it into \fBLOGIN.COM\fR.  To compile, type:
.IP
.B @MAKE
.sp .75i
.LP
.B "Good Luck ..."
.LP
Every effort has been made to make \fISparse\fR
a pleasant program to work with.  I hope you find it so.  If you have
any comments, criticisms, or praise, I would like to hear from you.  I 
would be especially be interested in reports on any bugs you may find,
however I cannot promise to reply. 

.nf
Ken Kundert
Dept. of Electrical Engineering and Computer Sciences
University of California
Berkeley, California 94720
June 1988

sparse@ic.berkeley.edu
ucbvax!ic!sparse
