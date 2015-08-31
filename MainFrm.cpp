// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// MainFrm.cpp : implementation of the CMainFrame class
//


#include "stdafx.h"
#include "FastHenry2.h"
#include "FastHenryDoc.h"
#include "FastHenryView.h"

#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern volatile BOOL g_bIsFCRunning;
extern CMainFrame *pFHMainFrame = NULL;
extern HWND FHHwnd;

// globals to communicate with working thread
char g_sTitle[MAX_TITLE_LENGHT];
char g_sOutputText[MAX_OUTPUT_TEXT_LEN];

/////////////////////////////////////////////////////////////////////////////
// CMainFrame

IMPLEMENT_DYNCREATE(CMainFrame, CFrameWnd)

static const UINT UWM_SET_TITLE = RegisterWindowMessage(_T("UWM_SET_TITLE-FastHenry-Enrico_Di_Lorenzo"));
static const UINT UWM_OUTPUT_TEXT = RegisterWindowMessage(_T("UWM_OUTPUT_TEXT-FastHenry-Enrico_Di_Lorenzo"));
static const UINT UWM_LOG_TEXT = RegisterWindowMessage(_T("UWM_LOG_TEXT-FastHenry-Enrico_Di_Lorenzo"));
static const UINT UWM_THREAD_END = RegisterWindowMessage(_T("UWM_THREAD_END-FastHenry-Enrico_Di_Lorenzo"));

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
	//{{AFX_MSG_MAP(CMainFrame)
	ON_WM_CREATE()
	ON_WM_CLOSE()

	// catch messages from 'remote' control console 
	ON_WM_COPYDATA()

	//}}AFX_MSG_MAP
	// Global help commands
	ON_COMMAND(ID_HELP, OnHelp)

	ON_REGISTERED_MESSAGE(UWM_SET_TITLE, OnSetTitle)
	ON_REGISTERED_MESSAGE(UWM_OUTPUT_TEXT, OnOutputText)
	ON_REGISTERED_MESSAGE(UWM_LOG_TEXT, OnLogText)
	ON_REGISTERED_MESSAGE(UWM_THREAD_END, OnThreadEnd)

END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

/////////////////////////////////////////////////////////////////////////////
// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	// TODO: add member initialization code here

	// global object visible also to C modules
	pFHMainFrame = this;
	
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	if (!m_wndToolBar.Create(this) ||
		!m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

	if (!m_wndStatusBar.Create(this) ||
		!m_wndStatusBar.SetIndicators(indicators,
		  sizeof(indicators)/sizeof(UINT)))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}

	// TODO: Remove this if you don't want tool tips or a resizeable toolbar
	m_wndToolBar.SetBarStyle(m_wndToolBar.GetBarStyle() |
		CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);

	// TODO: Delete these three lines if you don't want the toolbar to
	//  be dockable
	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);

	FHHwnd = GetSafeHwnd();
	ASSERT(IsWindow(FHHwnd));

	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	//  Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	// make program name appear before file name
	cs.style &= ~FWS_PREFIXTITLE;

	return CFrameWnd::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWnd::Dump(dc);
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMainFrame message handlers

void CMainFrame::OnClose() 
{
	CFastHenryView *view;

	if( g_bIsFCRunning == TRUE ) {
		view = (CFastHenryView *) GetActiveView();
		view->OutputText("Error on user 'Close' request: cannot close window while FastHenry is running\n", FHV_RED);
	}
	else {

		CFastHenryDoc *doc;

		// set the document as unmodified, not to have the prompt
		// for changes saving
		doc = (CFastHenryDoc *) GetActiveDocument();
		doc->SetModifiedFlag(FALSE);

		CFrameWnd::OnClose();
	}
}

// process working thread request to change main window title
afx_msg LRESULT CMainFrame::OnSetTitle(WPARAM, LPARAM)
{
	CFastHenryDoc *doc;

	doc = (CFastHenryDoc *) GetActiveDocument();

	if( g_sTitle == NULL || strcmp(g_sTitle, "-") == 0 ) {
		doc->SetTitle(_T(""));
	}
	else {
		doc->SetTitle( (LPCTSTR) g_sTitle);
	}


	return 0;
}

afx_msg LRESULT CMainFrame::OnOutputText(WPARAM color, LPARAM)
{
	CFastHenryView *view;

	view = (CFastHenryView *) GetActiveView();

	view->OutputText(g_sOutputText, color);

	// store string in log list
	m_clsLogList.push_front(COUPLE(g_sOutputText, color));

	return strlen(g_sOutputText);
}

afx_msg LRESULT CMainFrame::OnLogText(WPARAM color, LPARAM)
{
	CFastHenryDoc *pDoc;
	CString logString;

	pDoc = (CFastHenryDoc *) GetActiveDocument();

	// empty log list
	while(m_clsLogList.size() > 0) {
		// call callback routine, passing the oldest element of the
		// log list as argument
		pDoc->LogCallback(m_clsLogList.back().first, m_clsLogList.back().second);
		// remove element from the list
		m_clsLogList.pop_back();
	}
	
	return 0;
}

afx_msg LRESULT CMainFrame::OnThreadEnd(WPARAM, LPARAM)
{
	CFastHenryDoc *pDoc;
	
	pDoc = (CFastHenryDoc *) GetActiveDocument();
	// call callback routine
	pDoc->EndCallback();
	
	return 0;
}

BOOL CMainFrame::OnCopyData(CWnd* pWnd, COPYDATASTRUCT* pCopyDataStruct) 
{

	// call main function in CFastHenryApp
	((CFastHenryApp *)AfxGetApp())->OnCopyData(pWnd, pCopyDataStruct);

	return CFrameWnd::OnCopyData(pWnd, pCopyDataStruct);
}

void CMainFrame::OnHelp() 
{
	CString path;

	// get install directory path
	path = ((CFastHenryApp*)AfxGetApp())->m_strAppPath;
	// and form full application path name
	path += _T("\\FastHenry2.chm::/WelcometoFastHenry.htm");


	HtmlHelp(NULL, path, HH_DISPLAY_TOPIC, 0);
}
