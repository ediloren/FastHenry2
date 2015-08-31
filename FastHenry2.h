// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// FastHenry2.h : main header file for the FASTHENRY2 application
//

#if !defined(AFX_FASTHENRY2_H__E1385E85_3037_11D5_9282_74D314C10000__INCLUDED_)
#define AFX_FASTHENRY2_H__E1385E85_3037_11D5_9282_74D314C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols
#include "RunDialog.h"
#include "FastHenryDoc.h"

/////////////////////////////////////////////////////////////////////////////
// CFastHenryApp:
// See FastHenry2.cpp for the implementation of this class
//

class CFastHenryApp : public CWinApp
{
public:
	CFastHenryApp();
	~CFastHenryApp();
	void freeCmdLine();
	BOOL OnCopyData(CWnd* pWnd, COPYDATASTRUCT *pCopyDataStruct);
	BOOL parseCmdLine(int *argc, char ***argv, const char *commandStr);
	void LaunchFastHenry();

	int argc;
	char **argv;
	// result matrices
	COleSafeArray m_clsResMatrix;
	COleSafeArray m_clsIndMatrix;
	COleSafeArray m_clsFrequencies;
	COleSafeArray m_clsRowPortNames;
	COleSafeArray m_clsColPortNames;
	// install path
	CString m_strAppPath;
	// default run path
	CString m_strSamplePath;

protected:
	int getSubstring(const char *buffer, char *substr, int *skip);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFastHenryApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

	void RunDialogLaunch();

// Implementation
private:

	CRunDialog m_clsRunDlg;
	// Server object for document creation
	COleTemplateServer m_server;
	// document template
	CSingleDocTemplate* m_pDocTemplate;

	//{{AFX_MSG(CFastHenryApp)
	afx_msg void OnAppAbout();
	afx_msg void OnFastHenryRun();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FASTHENRY2_H__E1385E85_3037_11D5_9282_74D314C10000__INCLUDED_)
