// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// C Header file with definitions for linking FastHenry with Windows
// I/O window


#ifndef FHWINDOW_H
#define FHWINDOW_H

#include "FHStructs.h"

// defines used with FHOnClosing function defined in FastHenryWindow.h/.cpp
#define FH_NORMAL_END 0
#define FH_USER_BREAK 1
#define FH_GENERIC_ERROR 1000

// maximum length of a row/port name
#define FH_MAX_NAME_LEN		256

extern volatile char bFHContinue;
extern volatile struct impMatrix strctImpMatrix;

extern void FHSetName(char *name);

#endif

