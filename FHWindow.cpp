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

#include "FHWindow.h"
#include "MainFrm.h"
#include "FastHenryDoc.h"
#include "FastHenryView.h"
#include <stdio.h>

extern CMainFrame *pFHMainFrame;

void FHSetName( char *name)
{

	CFastHenryDoc *doc;
	CFastHenryView *view;

	// retrieve the document
	view = (CFastHenryView *) pFHMainFrame->GetActiveView(); 
	doc = (CFastHenryDoc *) view->GetDocument();
	
	if( name == NULL || strcmp(name, "-") == 0 ) {
		doc->SetTitle(_T(""));
	}
	else {
		doc->SetTitle( (LPCTSTR) name);
	}
}