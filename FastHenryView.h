// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// FastHenryView.h : interface of the CFastHenryView class
//
/////////////////////////////////////////////////////////////////////////////


#if !defined(AFX_FASTHENRYVIEW_H__E1385E8D_3037_11D5_9282_74D314C10000__INCLUDED_)
#define AFX_FASTHENRYVIEW_H__E1385E8D_3037_11D5_9282_74D314C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "FHStructs.h"

class CFastHenry2CntrItem;

class CFastHenryView : public CRichEditView
{
protected: // create from serialization only
	CFastHenryView();
	DECLARE_DYNCREATE(CFastHenryView)

// Attributes
public:
	CFastHenryDoc* GetDocument();

// Operations
public:
	void SetColor(long color);
	void OutputText(char *text, unsigned long color);

private:
	void SetDefaultCharFormat();
	void SetDefaultParaFormat();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFastHenryView)
	public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual void OnInitialUpdate(); // called first time after construct
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

	// override to select what can be pasted in the view
	BOOL CanPaste() const;

// Implementation
public:
	virtual ~CFastHenryView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CFastHenryView)
	afx_msg void OnDestroy();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnFasthenryStop();
	afx_msg void OnEditClearall();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in FastHenryView.cpp
inline CFastHenryDoc* CFastHenryView::GetDocument()
   { return (CFastHenryDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FASTHENRYVIEW_H__E1385E8D_3037_11D5_9282_74D314C10000__INCLUDED_)
