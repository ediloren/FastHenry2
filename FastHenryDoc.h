// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// FastHenryDoc.h : interface of the CFastHenryDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_FASTHENRYDOC_H__E1385E8B_3037_11D5_9282_74D314C10000__INCLUDED_)
#define AFX_FASTHENRYDOC_H__E1385E8B_3037_11D5_9282_74D314C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

// Unicode string class 
#include "UnicodeString.h"

class CFastHenryDoc : public CRichEditDoc
{
protected: // create from serialization only
	CFastHenryDoc();
	DECLARE_DYNCREATE(CFastHenryDoc)

// Attributes
public:
	void EndCallback();
	void LogCallback(CUnicodeString &text, unsigned long color);

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFastHenryDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	protected:
	virtual BOOL SaveModified();
	//}}AFX_VIRTUAL
	virtual CRichEditCntrItem* CreateClientItem(REOBJECT* preo) const;

// Implementation
public:
	virtual ~CFastHenryDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	HRESULT AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...);

	IDispatch *m_pDispEndCallback;
	CUnicodeString m_clsEndCallbackName;
	IDispatch *m_pDispLogCallback;
	CUnicodeString m_clsLogCallbackName;

// Generated message map functions
protected:
	//{{AFX_MSG(CFastHenryDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

	// Generated OLE dispatch map functions
	//{{AFX_DISPATCH(CFastHenryDoc)
	afx_msg BOOL IsRunning();
	afx_msg void ShowWindow();
	afx_msg BOOL Quit();
	afx_msg VARIANT GetResistance();
	afx_msg VARIANT GetInductance();
	afx_msg VARIANT GetRowPortNames();
	afx_msg VARIANT GetColPortNames();
	afx_msg VARIANT GetFrequencies();
	afx_msg void Stop();
	afx_msg BOOL SetEndCallback(LPDISPATCH callback, LPCTSTR cbName);
	afx_msg BOOL SetLogCallback(LPDISPATCH callback, LPCTSTR cbName);
	afx_msg BOOL Run(LPCTSTR commandLine);
	//}}AFX_DISPATCH
	DECLARE_DISPATCH_MAP()
	DECLARE_INTERFACE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FASTHENRYDOC_H__E1385E8B_3037_11D5_9282_74D314C10000__INCLUDED_)
