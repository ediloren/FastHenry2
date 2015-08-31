// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// Functions linking FastHenry with Windows I/O window


#include "stdafx.h"

#include "FastHenryWindow.h"
#include <stdio.h>

extern CFastHenryApp *mainApp = NULL;

void FHSetName(char *name)
{
	LRESULT val;

	static const UINT UWM_SET_TITLE = RegisterWindowMessage(_T("UWM_SET_TITLE-FastHenry-Enrico_Di_Lorenzo"));

	// safe copy of FastHenry input file name
	strncpy(g_sTitle, name, MAX_TITLE_LENGHT);
	g_sTitle[MAX_TITLE_LENGHT-1] = 0;
	
	// ask main process to update window title
	if( IsWindow(FHHwnd) ) 
		val = ::SendMessage(FHHwnd, UWM_SET_TITLE, 0, 0);
}

// function called on FastHenry closing (a return() will follow 
// the call to this routine)
void FHOnClosing(int cause)
{
	if(cause == FH_USER_BREAK) {
		viewprintf(stderr, "\n\nWARNING: Program execution stopped on user request!");
	}
}

// function called to exit FastHenry (no any other function call can 
// follow: thread exits here)
void FHExit(int cause)
{
	if(cause == FH_USER_BREAK) {
		viewprintf(stderr, "\n\nWARNING: Program execution stopped on user request!");
	}

	// if exiting thread directly, won't pass through RunFHThread() in FastHenry2
	// any more. Therefore must clean up here
	OnFastHenryExit();

	// close all open files
	_fcloseall();

	AfxEndThread(0);
}

int viewprintf(FILE *out, const char *fmt,...)
{
	int ret;
	unsigned long color;
	LRESULT val;

	static const UINT UWM_OUTPUT_TEXT = RegisterWindowMessage(_T("UWM_OUTPUT_TEXT-FastHenry-Enrico_Di_Lorenzo"));
	static const UINT UWM_LOG_TEXT = RegisterWindowMessage(_T("UWM_LOG_TEXT-FastHenry-Enrico_Di_Lorenzo"));

	va_list arg_ptr;

	va_start(arg_ptr, fmt);

	ret = vsprintf(g_sOutputText, fmt, arg_ptr);

	va_end(arg_ptr);

	if (out == stdout) {
		color = FHV_BLACK;
	}
	else
		color = FHV_RED;

	// ask main process to output text in given color
	if(IsWindow(FHHwnd)) {
		// SendMessage() does not return until the window procedure 
		// has processed the message, so no synchronization issues
		val = ::SendMessage(FHHwnd, UWM_OUTPUT_TEXT, color, 0);
		// However, automation cannot make calls during input-synchronous calls,
		// like SendMessage(); so post a message to the main process to ask
		// for Automation log. The log message is stored during the processing
		// of the previous SendMessage()
		val = ::PostMessage(FHHwnd, UWM_LOG_TEXT, color, 0);
	}

	return ret;
}
