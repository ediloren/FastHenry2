; CLW file contains information for the MFC ClassWizard

[General Info]
Version=1
LastClass=CFastHenryDoc
LastTemplate=CDialog
NewFileInclude1=#include "stdafx.h"
NewFileInclude2=#include "FastHenry2.h"
LastPage=0

ClassCount=6
Class1=CFastHenryApp
Class2=CFastHenryDoc
Class3=CFastHenryView
Class4=CMainFrame

ResourceCount=8
Resource1=IDR_MAINFRAME (English (U.S.))
Resource2=IDR_MAINFRAME
Resource4=IDR_CNTR_INPLACE (English (U.S.))
Class5=CAboutDlg
Resource3=IDD_ABOUTBOX
Resource5=IDD_ABOUTBOX (English (U.S.))
Resource6=IDD_DIALOG_RUN
Class6=CRunDialog
Resource7=IDR_CNTR_INPLACE
Resource8=IDD_DIALOG_RUN (English (U.S.))

[CLS:CFastHenryApp]
Type=0
HeaderFile=FastHenry2.h
ImplementationFile=FastHenry2.cpp
Filter=W
LastObject=ID_FASTHENRY_RUN
BaseClass=CWinApp
VirtualFilter=AC

[CLS:CFastHenryDoc]
Type=0
HeaderFile=FastHenryDoc.h
ImplementationFile=FastHenryDoc.cpp
Filter=N
BaseClass=CRichEditDoc
VirtualFilter=DC
LastObject=CFastHenryDoc

[CLS:CFastHenryView]
Type=0
HeaderFile=FastHenryView.h
ImplementationFile=FastHenryView.cpp
Filter=W
LastObject=CFastHenryView
BaseClass=CRichEditView
VirtualFilter=VWC

[CLS:CMainFrame]
Type=0
HeaderFile=MainFrm.h
ImplementationFile=MainFrm.cpp
Filter=T
BaseClass=CFrameWnd
VirtualFilter=fWC
LastObject=CMainFrame



[CLS:CAboutDlg]
Type=0
HeaderFile=FastHenry2.cpp
ImplementationFile=FastHenry2.cpp
Filter=D

[DLG:IDD_ABOUTBOX]
Type=1
Class=CAboutDlg
ControlCount=8
Control1=IDC_STATIC,static,1342177283
Control2=IDC_STATIC,static,1342308480
Control3=IDC_STATIC,static,1342308352
Control4=IDOK,button,1342373889
Control5=IDC_STATIC,static,1342308352
Control6=IDC_STATIC,static,1342308352
Control7=IDC_STATIC,static,1342308352
Control8=IDC_STATIC,static,1342308352

[MNU:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FASTHENRY_RUN
Command2=ID_FILE_PRINT
Command3=ID_FILE_PRINT_PREVIEW
Command4=ID_FILE_PRINT_SETUP
Command5=ID_APP_EXIT
Command6=ID_EDIT_COPY
Command7=ID_EDIT_CLEARALL
Command8=ID_EDIT_SELECT_ALL
Command9=ID_EDIT_FIND
Command10=ID_EDIT_REPEAT
Command11=ID_VIEW_TOOLBAR
Command12=ID_VIEW_STATUS_BAR
Command13=ID_FASTHENRY_RUN
Command14=ID_FASTHENRY_STOP
Command15=ID_HELP
Command16=ID_APP_ABOUT
CommandCount=16

[MNU:IDR_CNTR_INPLACE]
Type=1
Class=CFastHenryView
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_FILE_SAVE_AS
Command5=ID_FILE_PRINT
Command6=ID_FILE_PRINT_PREVIEW
Command7=ID_FILE_PRINT_SETUP
Command8=ID_FILE_MRU_FILE1
Command9=ID_APP_EXIT
CommandCount=9

[ACL:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_EDIT_SELECT_ALL
Command2=ID_EDIT_COPY
Command3=ID_EDIT_FIND
Command4=ID_EDIT_REPLACE
Command5=ID_FASTHENRY_RUN
Command6=ID_FILE_PRINT
Command7=ID_EDIT_PASTE
Command8=ID_EDIT_UNDO
Command9=ID_EDIT_CUT
Command10=ID_HELP
Command11=ID_CONTEXT_HELP
Command12=ID_EDIT_REPEAT
Command13=ID_EDIT_COPY
Command14=ID_EDIT_PASTE
Command15=ID_OLE_EDIT_PROPERTIES
Command16=ID_EDIT_CUT
Command17=ID_EDIT_UNDO
CommandCount=17

[ACL:IDR_CNTR_INPLACE]
Type=1
Class=CFastHenryView
Command1=ID_FASTHENRY_RUN
Command2=ID_FILE_PRINT
Command3=ID_HELP
Command4=ID_CONTEXT_HELP
CommandCount=4

[DLG:IDD_ABOUTBOX (English (U.S.))]
Type=1
Class=CAboutDlg
ControlCount=7
Control1=IDC_STATIC,static,1342177283
Control2=IDC_STATIC,static,1342308480
Control3=IDC_STATIC,static,1342308352
Control4=IDOK,button,1342373889
Control5=IDC_STATIC,static,1342308352
Control6=IDC_STATIC,static,1342308352
Control7=IDC_STATIC,static,1342308352

[TB:IDR_MAINFRAME (English (U.S.))]
Type=1
Class=?
Command1=ID_FASTHENRY_RUN
Command2=ID_EDIT_COPY
Command3=ID_APP_ABOUT
Command4=ID_CONTEXT_HELP
CommandCount=4

[MNU:IDR_MAINFRAME (English (U.S.))]
Type=1
Class=CMainFrame
Command1=ID_FASTHENRY_RUN
Command2=ID_FILE_PRINT
Command3=ID_FILE_PRINT_PREVIEW
Command4=ID_FILE_PRINT_SETUP
Command5=ID_APP_EXIT
Command6=ID_EDIT_COPY
Command7=ID_EDIT_CLEARALL
Command8=ID_EDIT_SELECT_ALL
Command9=ID_EDIT_FIND
Command10=ID_EDIT_REPEAT
Command11=ID_VIEW_TOOLBAR
Command12=ID_VIEW_STATUS_BAR
Command13=ID_FASTHENRY_RUN
Command14=ID_FASTHENRY_STOP
Command15=ID_HELP
Command16=ID_APP_ABOUT
CommandCount=16

[MNU:IDR_CNTR_INPLACE (English (U.S.))]
Type=1
Class=?
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_FILE_SAVE_AS
Command5=ID_FILE_PRINT
Command6=ID_FILE_PRINT_PREVIEW
Command7=ID_FILE_PRINT_SETUP
Command8=ID_FILE_MRU_FILE1
Command9=ID_APP_EXIT
CommandCount=9

[ACL:IDR_MAINFRAME (English (U.S.))]
Type=1
Class=?
Command1=ID_EDIT_SELECT_ALL
Command2=ID_EDIT_COPY
Command3=ID_EDIT_FIND
Command4=ID_EDIT_REPLACE
Command5=ID_FASTHENRY_RUN
Command6=ID_FILE_PRINT
Command7=ID_EDIT_PASTE
Command8=ID_EDIT_UNDO
Command9=ID_EDIT_CUT
Command10=ID_HELP
Command11=ID_CONTEXT_HELP
Command12=ID_EDIT_REPEAT
Command13=ID_EDIT_COPY
Command14=ID_EDIT_PASTE
Command15=ID_OLE_EDIT_PROPERTIES
Command16=ID_EDIT_CUT
Command17=ID_EDIT_UNDO
CommandCount=17

[ACL:IDR_CNTR_INPLACE (English (U.S.))]
Type=1
Class=?
Command1=ID_FASTHENRY_RUN
Command2=ID_FILE_PRINT
Command3=ID_HELP
Command4=ID_CONTEXT_HELP
CommandCount=4

[DLG:IDD_DIALOG_RUN (English (U.S.))]
Type=1
Class=CRunDialog
ControlCount=45
Control1=IDC_INPUT_FILE_NAME,edit,1350631552
Control2=IDC_BROWSE_INPUT_FILE,button,1342242816
Control3=IDC_STATIC4,static,1342308352
Control4=IDC_STATIC6,static,1342308352
Control5=IDC_STATIC1,static,1342308352
Control6=IDC_STATIC16,static,1342308352
Control7=IDC_STATIC15,static,1342308352
Control8=IDC_STATIC2,static,1342308352
Control9=IDC_STATIC3,static,1342308352
Control10=IDC_STATIC5,static,1342308352
Control11=IDC_ROM_ORDER,edit,1350631552
Control12=IDC_FILENAME_SUFFIX,edit,1350631552
Control13=IDC_MULTIPOLE_ORDER,edit,1350631552
Control14=IDC_TOLERANCE_RTOL,edit,1350631552
Control15=IDC_TOLERANCE_ATOL,edit,1350631552
Control16=IDC_SOLVE_ITERATIONS,edit,1350631552
Control17=IDC_INITIAL_REFINEMENT,edit,1350631552
Control18=IDC_SHELLS_RADIUS,edit,1350631552
Control19=IDC_STATIC7,static,1342308352
Control20=IDC_STATIC8,static,1342308352
Control21=IDC_STATIC9,static,1342308352
Control22=IDC_STATIC10,static,1342308352
Control23=IDC_STATIC11,static,1342308352
Control24=IDC_STATIC12,static,1342308352
Control25=IDC_STATIC13,static,1342308352
Control26=IDC_STATIC14,static,1342308352
Control27=IDC_STATIC20,static,1342308352
Control28=IDC_MATRIX_SOLUTION_METHOD,combobox,1344339971
Control29=IDC_MATRIX_VECTOR_PRODUCT,combobox,1344339971
Control30=IDC_PRECONDITIONER,combobox,1344339971
Control31=IDC_PARTITIONING_LEVELS,combobox,1342242818
Control32=IDC_VISUALIZATION_MODE,combobox,1344339971
Control33=IDC_GROUND_PLANE,combobox,1344339971
Control34=IDC_DUMP_MATRICES,combobox,1344339971
Control35=IDC_DUMP_FILE_TYPE,combobox,1344339971
Control36=IDC_COMPUTE_PORTS,edit,1352728772
Control37=IDC_AUTOMATIC_REFINE,button,1342242819
Control38=IDC_EXIT_AFTER_ROM,button,1342242819
Control39=IDC_REGURGITATE,button,1342242819
Control40=IDC_DEBUG_INFO,button,1342242819
Control41=IDCANCEL,button,1342242816
Control42=IDC_RESET,button,1342242816
Control43=IDOK,button,1342242817
Control44=IDC_STATIC17,static,1342308352
Control45=IDC_PRINT_RESULTS_ON_SCREEN,button,1342242819

[CLS:CRunDialog]
Type=0
HeaderFile=RunDialog.h
ImplementationFile=RunDialog.cpp
BaseClass=CDialog
Filter=D
LastObject=CRunDialog
VirtualFilter=dWC

[DLG:IDD_DIALOG_RUN]
Type=1
Class=CRunDialog
ControlCount=45
Control1=IDC_INPUT_FILE_NAME,edit,1350631552
Control2=IDC_BROWSE_INPUT_FILE,button,1342242816
Control3=IDC_STATIC4,static,1342308352
Control4=IDC_STATIC6,static,1342308352
Control5=IDC_STATIC1,static,1342308352
Control6=IDC_STATIC16,static,1342308352
Control7=IDC_STATIC15,static,1342308352
Control8=IDC_STATIC2,static,1342308352
Control9=IDC_STATIC3,static,1342308352
Control10=IDC_STATIC5,static,1342308352
Control11=IDC_ROM_ORDER,edit,1350631552
Control12=IDC_FILENAME_SUFFIX,edit,1350631552
Control13=IDC_MULTIPOLE_ORDER,edit,1350631552
Control14=IDC_TOLERANCE_RTOL,edit,1350631552
Control15=IDC_TOLERANCE_ATOL,edit,1350631552
Control16=IDC_SOLVE_ITERATIONS,edit,1350631552
Control17=IDC_INITIAL_REFINEMENT,edit,1350631552
Control18=IDC_SHELLS_RADIUS,edit,1350631552
Control19=IDC_STATIC7,static,1342308352
Control20=IDC_STATIC8,static,1342308352
Control21=IDC_STATIC9,static,1342308352
Control22=IDC_STATIC10,static,1342308352
Control23=IDC_STATIC11,static,1342308352
Control24=IDC_STATIC12,static,1342308352
Control25=IDC_STATIC13,static,1342308352
Control26=IDC_STATIC14,static,1342308352
Control27=IDC_STATIC20,static,1342308352
Control28=IDC_MATRIX_SOLUTION_METHOD,combobox,1344339971
Control29=IDC_MATRIX_VECTOR_PRODUCT,combobox,1344339971
Control30=IDC_PRECONDITIONER,combobox,1344339971
Control31=IDC_PARTITIONING_LEVELS,combobox,1342242818
Control32=IDC_VISUALIZATION_MODE,combobox,1344339971
Control33=IDC_GROUND_PLANE,combobox,1344339971
Control34=IDC_DUMP_MATRICES,combobox,1344339971
Control35=IDC_DUMP_FILE_TYPE,combobox,1344339971
Control36=IDC_COMPUTE_PORTS,edit,1352728772
Control37=IDC_AUTOMATIC_REFINE,button,1342242819
Control38=IDC_EXIT_AFTER_ROM,button,1342242819
Control39=IDC_REGURGITATE,button,1342242819
Control40=IDC_DEBUG_INFO,button,1342242819
Control41=IDCANCEL,button,1342242816
Control42=IDC_RESET,button,1342242816
Control43=IDOK,button,1342242817
Control44=IDC_STATIC17,static,1342308352
Control45=IDC_PRINT_RESULTS_ON_SCREEN,button,1342242819

[TB:IDR_MAINFRAME]
Type=1
Class=?
Command1=ID_FASTHENRY_RUN
Command2=ID_EDIT_COPY
Command3=ID_APP_ABOUT
Command4=ID_CONTEXT_HELP
CommandCount=4

