# Microsoft Developer Studio Project File - Name="FastHenry2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=FastHenry2 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "FastHenry2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "FastHenry2.mak" CFG="FastHenry2 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "FastHenry2 - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "FastHenry2 - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "FastHenry2 - Win32 Release"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 6
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "." /I "./sparse" /I "./FastHenry" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "WINDOWS" /D "NO_GETHOSTNAME" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /o /win32 "NUL"
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /o "FastHenry2.midl.log" /win32
# ADD BASE RSC /l 0x410 /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x410 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /machine:I386
# ADD LINK32 Htmlhelp.lib /nologo /subsystem:windows /machine:I386

!ELSEIF  "$(CFG)" == "FastHenry2 - Win32 Debug"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 6
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MDd /W3 /GX /ZI /Od /I "." /I "./sparse" /I "./FastHenry" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "NO_GETHOSTNAME" /D "WINDOWS" /YX /FD /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /o /win32 "NUL"
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /o "FastHenry2.midl.log" /win32
# ADD BASE RSC /l 0x410 /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x410 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 Htmlhelp.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "FastHenry2 - Win32 Release"
# Name "FastHenry2 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\CntrItem.cpp
# End Source File
# Begin Source File

SOURCE=.\FastHenry2.cpp
# End Source File
# Begin Source File

SOURCE=.\FastHenry2.odl
# End Source File
# Begin Source File

SOURCE=.\FastHenry2.rc
# End Source File
# Begin Source File

SOURCE=.\FastHenryDoc.cpp
# End Source File
# Begin Source File

SOURCE=.\FastHenryView.cpp
# End Source File
# Begin Source File

SOURCE=.\FastHenryWindow.cpp
# End Source File
# Begin Source File

SOURCE=.\MainFrm.cpp
# End Source File
# Begin Source File

SOURCE=.\RunDialog.cpp
# End Source File
# Begin Source File

SOURCE=.\StdAfx.cpp
# ADD CPP /Yc"stdafx.h"
# End Source File
# Begin Source File

SOURCE=.\UnicodeString.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\CntrItem.h
# End Source File
# Begin Source File

SOURCE=.\FastHenry2.h
# End Source File
# Begin Source File

SOURCE=.\FastHenryDoc.h
# End Source File
# Begin Source File

SOURCE=.\FastHenryView.h
# End Source File
# Begin Source File

SOURCE=.\FastHenryWindow.h
# End Source File
# Begin Source File

SOURCE=.\FHStructs.h
# End Source File
# Begin Source File

SOURCE=.\FHWindow.h
# End Source File
# Begin Source File

SOURCE=.\MainFrm.h
# End Source File
# Begin Source File

SOURCE=.\Resource.h
# End Source File
# Begin Source File

SOURCE=.\RunDialog.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# Begin Source File

SOURCE=.\UnicodeString.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;cnt;rtf;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\res\FastHenry2.ico
# End Source File
# Begin Source File

SOURCE=.\res\FastHenry2.rc2
# End Source File
# Begin Source File

SOURCE=.\res\FastHenryDoc.ico
# End Source File
# Begin Source File

SOURCE=.\res\Toolbar.bmp
# End Source File
# End Group
# Begin Group "FastHenry"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\FastHenry\addgroundplane.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\barnoldi.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\breakupSeg.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\calcp.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\capsolve.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\CMPLX.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\contact.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\cx_ludecomp.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\default_opts.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\deg_mutual.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\direct.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\dist_betw_fils.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\fillM.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\find_nonuni_path.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\findpaths.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\gmres.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\GP.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\hole.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\induct.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\INDUCT.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\joelself.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\math2.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\MATH2.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulDo.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulGlobal.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulGlobal.h
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulLocal.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulMats.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulMulti.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulSetup.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mulStruct.h
# End Source File
# Begin Source File

SOURCE=.\FastHenry\mutual.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\newPrecond.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\parse_command_line.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\PATRAN.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\prec_cost.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\precond.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\read_tree.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\readGeom.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\regurgitate.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\RESUSAGE.H
# End Source File
# Begin Source File

SOURCE=.\FastHenry\savemat_mod.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\setupComputePsi.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\setupMulti.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\uglieralloc.c
# End Source File
# Begin Source File

SOURCE=.\FastHenry\writefastcap.c
# End Source File
# End Group
# Begin Group "Sparse"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\sparse\spAllocate.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spBuild.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spCompat.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spConfig.h
# End Source File
# Begin Source File

SOURCE=.\sparse\spDefs.h
# End Source File
# Begin Source File

SOURCE=.\sparse\spFactor.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spFortran.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spMatrix.h
# End Source File
# Begin Source File

SOURCE=.\sparse\spOutput.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spSolve.c
# End Source File
# Begin Source File

SOURCE=.\sparse\spUtils.c
# End Source File
# End Group
# Begin Source File

SOURCE=.\hlp\fasthenry2.hhp

!IF  "$(CFG)" == "FastHenry2 - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
OutDir=.\Release
InputPath=.\hlp\fasthenry2.hhp
InputName=fasthenry2

"$(OutDir)\$(InputName).chm" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	hhc.exe .\hlp\$(InputName).hhp 
	if exist "hlp\$(InputName).chm" copy "hlp\$(InputName).chm" $(OutDir) 
	echo off 
	
# End Custom Build

!ELSEIF  "$(CFG)" == "FastHenry2 - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
OutDir=.\Debug
InputPath=.\hlp\fasthenry2.hhp
InputName=fasthenry2

"$(OutDir)\$(InputName).chm" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	hhc.exe .\hlp\$(InputName).hhp 
	if exist "hlp\$(InputName).chm" copy "hlp\$(InputName).chm" $(OutDir) 
	echo off 
	
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FastHenry2.reg
# End Source File
# Begin Source File

SOURCE=.\HISTORY\History.txt
# End Source File
# Begin Source File

SOURCE=.\HISTORY\History_detailed.txt
# End Source File
# Begin Source File

SOURCE=.\License_and_history\license.txt
# End Source File
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project
