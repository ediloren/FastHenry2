// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// FastHenryView.cpp : implementation of the CFastHenryView class
//


#include "stdafx.h"
#include "FastHenry2.h"

#include "FastHenryDoc.h"
#include "FastHenryView.h"

#include <stdio.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// globals visible to C functions
extern "C" volatile char bFHContinue;

extern "C" int viewprintf(FILE *out, const char *fmt,...);
extern "C" void FHSetName( char *name);

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView

IMPLEMENT_DYNCREATE(CFastHenryView, CRichEditView)

BEGIN_MESSAGE_MAP(CFastHenryView, CRichEditView)
	//{{AFX_MSG_MAP(CFastHenryView)
	ON_WM_DESTROY()
	ON_WM_KEYDOWN()
	ON_WM_CHAR()
	ON_COMMAND(ID_FASTHENRY_STOP, OnFasthenryStop)
	ON_COMMAND(ID_EDIT_CLEARALL, OnEditClearall)
	//}}AFX_MSG_MAP

	// Standard printing commands
	//ON_COMMAND(ID_FILE_PRINT, CRichEditView::OnFilePrint)
	//ON_COMMAND(ID_FILE_PRINT_DIRECT, CRichEditView::OnFilePrint)
	//ON_COMMAND(ID_FILE_PRINT_PREVIEW, CRichEditView::OnFilePrintPreview)

END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView construction/destruction

CFastHenryView::CFastHenryView()
{
	// add init code
}

CFastHenryView::~CFastHenryView()
{
}

BOOL CFastHenryView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CRichEditView::PreCreateWindow(cs);
}

void CFastHenryView::OnInitialUpdate() 
{
	long len;

	CRichEditView::OnInitialUpdate();
	
	// set the desired char and paragraph formats
	GetRichEditCtrl().SetSel(0,-1);
	SetDefaultCharFormat();
	SetDefaultParaFormat();	
	len = GetTextLength();
	GetRichEditCtrl().SetSel(len,len);

	// this option are needed to set auto vertical scroll (and word wrap
	// at end of line)
	GetRichEditCtrl().SetOptions(ECOOP_SET, ECO_AUTOVSCROLL);
	m_nWordWrap = WrapToWindow;
	WrapChanged();

	// Set the printing margins (720 twips = 1/2 inch).
	SetMargins(CRect(720, 720, 720, 720));
}

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView printing

BOOL CFastHenryView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}


void CFastHenryView::OnDestroy()
{
	// Deactivate the item on destruction; this is important
	// when a splitter view is being used.
   COleClientItem* pActiveItem = GetDocument()->GetInPlaceActiveItem(this);
   if (pActiveItem != NULL && pActiveItem->GetActiveView() == this)
   {
      pActiveItem->Deactivate();
      ASSERT(GetDocument()->GetInPlaceActiveItem(this) == NULL);
   }
   CRichEditView::OnDestroy();
}


/////////////////////////////////////////////////////////////////////////////
// CFastHenryView diagnostics

#ifdef _DEBUG
void CFastHenryView::AssertValid() const
{
	CRichEditView::AssertValid();
}

void CFastHenryView::Dump(CDumpContext& dc) const
{
	CRichEditView::Dump(dc);
}

CFastHenryDoc* CFastHenryView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CFastHenryDoc)));
	return (CFastHenryDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView message handlers

void CFastHenryView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	// TODO: Add your message handler code here and/or call default
	
	//view must be read only
	//CRichEditView::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CFastHenryView::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	// TODO: Add your message handler code here and/or call default
	
	//view must be read only
	//CRichEditView::OnChar(nChar, nRepCnt, nFlags);
}

void CFastHenryView::OnFasthenryStop() 
{
	// signal to the working thread to stop execution
	bFHContinue	= FALSE;
}

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView overrides

BOOL CFastHenryView::CanPaste() const
{
	// no pasting
	return FALSE;
}

/////////////////////////////////////////////////////////////////////////////
// CFastHenryView operations

void CFastHenryView::SetDefaultCharFormat()
{
	CHARFORMAT2 cf;
	
	cf = GetCharFormatSelection();

	// reset all possible char formats to the only ones desired
	cf.dwMask = CFM_BOLD | CFM_CHARSET | CFM_COLOR | CFM_FACE | CFM_ITALIC
		| CFM_OFFSET | CFM_PROTECTED | CFM_STRIKEOUT | CFM_UNDERLINE;
	cf.dwEffects = 0;
	cf.yOffset = 0;
	cf.crTextColor = RGB(0, 0, 0); // color black
	cf.bCharSet = ANSI_CHARSET;
	cf.bPitchAndFamily = DEFAULT_PITCH | FF_MODERN;
	strcpy(cf.szFaceName, "Courier New");

	// set the char format on current selection
	SetCharFormat(cf);
	// set the char format as default for the view (any new text will 
	// have this format)
	GetRichEditCtrl().SetDefaultCharFormat(cf);
}

void CFastHenryView::SetDefaultParaFormat()
{
	PARAFORMAT2 pf;

	pf = GetParaFormatSelection();
	
	pf.dwMask = PFM_ALIGNMENT | PFM_NUMBERING |  PFM_OFFSET | PFM_OFFSETINDENT
		| PFM_RIGHTINDENT | PFM_STARTINDENT |  PFM_TABSTOPS;

	pf.wNumbering = 0;
	pf.dxStartIndent = 0;
	pf.dxRightIndent = 10;
	pf.dxOffset = 0;
	pf.wAlignment = PFA_LEFT;
	pf.cTabCount = 0;
	pf.cTabCount = 0;

	SetParaFormat(pf);
}

void CFastHenryView::SetColor(long color)
{
	CHARFORMAT2 cf;
	
	cf = GetCharFormatSelection();

	// reset all possible char formats to the only ones desired
	cf.dwMask = CFM_COLOR;
	cf.crTextColor = color; // desired color 

	// set the char format on current selection
	SetCharFormat(cf);
}

void CFastHenryView::OnEditClearall() 
{
	GetRichEditCtrl().SetSel(0, -1);
	GetRichEditCtrl().ReplaceSel(_T(""));
	HideCaret();
}

// Remark: no synchronization object is needed, since OutputText()
// calls are issued only by the main process, as an answer to a user
// action or upon request by the working thread through a message,
// not directly calling the function
void CFastHenryView::OutputText(char *text, unsigned long color)
{
	long lastPos;

	SetColor(color);
	lastPos = GetTextLength();
	GetRichEditCtrl().SetSel(lastPos, lastPos);
	GetRichEditCtrl().ReplaceSel(text);
	HideCaret();
}

/////////////////////////////////////////////////////////////////////////////
// globals


