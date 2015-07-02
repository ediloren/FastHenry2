$ def lnk$library sys$library:vaxcrtl.olb
$ cc spAllocate.c
$ cc spBuild.c
$ cc spCompat.c
$ cc spFactor.c
$ cc spFortran.c
$ cc spOutput.c
$ cc spSolve.c
$ cc spTest.c
$ cc spUtils.c
$ library/create/object sparse.olb -
spAllocate, -
spBuild, -
spCompat, -
spFactor, -
spFortran, -
spOutput, -
spSolve, -
spUtils
$ link/exe=sparse spTest.obj, sparse.olb/library
$ sparse :== $$disk2:[kundert.sparse]sparse.exe
