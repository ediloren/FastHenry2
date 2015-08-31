// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// Copyright (c) 2012 Enrico Di Lorenzo, www.fastfieldsolvers.com
// All rights reserved

#if !defined(AFX_RUNDIALOG_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_RUNDIALOG_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// RunDialog.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CRunDialog dialog

class CRunDialog : public CDialog
{
// Construction
public:
	CRunDialog(CWnd* pParent = NULL);   // standard constructor
	char CRunDialog::generateCmdLine(int *argc, char ***argv);

	CString m_strSamplePath;

// Dialog Data
	//{{AFX_DATA(CRunDialog)
	enum { IDD = IDD_DIALOG_RUN };
	CEdit	m_ctrlInputFileName;
	BOOL	m_bAutoRefinement;
	CString	m_sComputePorts;
	BOOL	m_bDebugInfo;
	int		m_iDumpFileType;
	int		m_iDumpMatrices;
	BOOL	m_bExitAfterROM;
	CString	m_sFilenameSuffix;
	int		m_iGndPlane;
	int		m_iInitialRefinement;
	int		m_iMatrixSolution;
	int		m_iMtxVctProduct;
	int		m_iMultipoleOrder;
	CString	m_sPartitioningLevels;
	int		m_iPrecondMethod;
	BOOL	m_bRegurgitate;
	CString	m_sRestrictToPorts;
	int		m_iROMOrder;
	int		m_iShellRadius;
	int		m_iSolveIterations;
	float	m_fATol;
	float	m_fRTol;
	int		m_iVisualizationMode;
	CString	m_sInputFileName;
	BOOL	m_bPrintResults;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CRunDialog)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CRunDialog)
	virtual BOOL OnInitDialog();
	afx_msg void OnReset();
	afx_msg void OnBrowseInputFile();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_RUNDIALOG_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
