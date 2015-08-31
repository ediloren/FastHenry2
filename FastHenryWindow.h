// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// CPP header file for module containing functions linking
// FastHenry with Windows I/O window


#ifndef FASTHENRYWINDOW_H
#define FASTHENRYWINDOW_H

#include <winuser.h>

#include "FastHenry2.h"
#include "FHStructs.h"

#define FH_NORMAL_END 0
#define FH_USER_BREAK 1
#define FH_GENERIC_ERROR 1000

#define MAX_TITLE_LENGHT 64
#define MAX_OUTPUT_TEXT_LEN 2048

// exported variables & functions

extern volatile char bFHContinue;
extern volatile BOOL g_bIsFCRunning;
extern HWND FHHwnd;
extern char g_sTitle[64];
extern char g_sOutputText[MAX_OUTPUT_TEXT_LEN];

extern "C" void FHSetName(char *name);
extern "C" void FHOnClosing(int cause);
extern "C" void FHExit(int cause);
extern "C" int viewprintf(FILE *out, const char *fmt,...);

void OnFastHenryExit();

#endif

