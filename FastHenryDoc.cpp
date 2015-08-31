// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// FastHenryDoc.cpp : implementation of the CFastHenryDoc class
//

#include "stdafx.h"

// include for _bstrt_t class
#include <comdef.h>

#include "FastHenry2.h"
#include "FastHenryDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern volatile BOOL g_bIsFCRunning;
extern "C" volatile char bFHContinue;
extern HWND FHHwnd;

/////////////////////////////////////////////////////////////////////////////
// CFastHenryDoc

IMPLEMENT_DYNCREATE(CFastHenryDoc, CRichEditDoc)

BEGIN_MESSAGE_MAP(CFastHenryDoc, CRichEditDoc)
	//{{AFX_MSG_MAP(CFastHenryDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
	// Enable default OLE container implementation
	ON_UPDATE_COMMAND_UI(ID_OLE_EDIT_LINKS, CRichEditDoc::OnUpdateEditLinksMenu)
	ON_COMMAND(ID_OLE_EDIT_LINKS, CRichEditDoc::OnEditLinks)
	ON_UPDATE_COMMAND_UI(ID_OLE_VERB_FIRST, CRichEditDoc::OnUpdateObjectVerbMenu)
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(CFastHenryDoc, CRichEditDoc)
	//{{AFX_DISPATCH_MAP(CFastHenryDoc)
	DISP_FUNCTION(CFastHenryDoc, "IsRunning", IsRunning, VT_BOOL, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "ShowWindow", ShowWindow, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "Quit", Quit, VT_BOOL, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "GetResistance", GetResistance, VT_VARIANT, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "GetInductance", GetInductance, VT_VARIANT, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "GetRowPortNames", GetRowPortNames, VT_VARIANT, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "GetColPortNames", GetColPortNames, VT_VARIANT, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "GetFrequencies", GetFrequencies, VT_VARIANT, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "Stop", Stop, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(CFastHenryDoc, "SetEndCallback", SetEndCallback, VT_BOOL, VTS_DISPATCH VTS_BSTR)
	DISP_FUNCTION(CFastHenryDoc, "SetLogCallback", SetLogCallback, VT_BOOL, VTS_DISPATCH VTS_BSTR)
	DISP_FUNCTION(CFastHenryDoc, "Run", Run, VT_BOOL, VTS_BSTR)
	//}}AFX_DISPATCH_MAP
END_DISPATCH_MAP()

// Note: we add support for IID_IFastHenry to support typesafe binding
//  from VBA.  This IID must match the GUID that is attached to the 
//  dispinterface in the .ODL file.

// {3BAD3C03-3868-4C3C-B265-BDAE2A4FA4F5}
static const IID IID_IFastHenry =
{ 0x3bad3c03, 0x3868, 0x4c3c, { 0xb2, 0x65, 0xbd, 0xae, 0x2a, 0x4f, 0xa4, 0xf5 } };

BEGIN_INTERFACE_MAP(CFastHenryDoc, CRichEditDoc)
	INTERFACE_PART(CFastHenryDoc, IID_IFastHenry, Dispatch)
END_INTERFACE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFastHenryDoc construction/destruction

CFastHenryDoc::CFastHenryDoc()
{
	// OLE init
	EnableAutomation();
	// OLE
	AfxOleLockApp();
	// IDispatch pointers
	m_pDispEndCallback = NULL;
	m_pDispLogCallback = NULL;
}

CFastHenryDoc::~CFastHenryDoc()
{
	if(m_pDispEndCallback != NULL) {
		// tell Automation to remove the reference to this object,
		// so memory can be released
		m_pDispEndCallback->Release();
	}
	if(m_pDispLogCallback != NULL) {
		// tell Automation to remove the reference to this object,
		// so memory can be released
		m_pDispLogCallback->Release();
	}

	// Close the OLE Library
	AfxOleUnlockApp();
}

BOOL CFastHenryDoc::OnNewDocument()
{
	if (!CRichEditDoc::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}

CRichEditCntrItem* CFastHenryDoc::CreateClientItem(REOBJECT* preo) const
{
	// cast away constness of this
	return NULL; //return new CFastHenry2CntrItem(preo, (CFastHenryDoc*) this);
}



/////////////////////////////////////////////////////////////////////////////
// CFastHenryDoc serialization

void CFastHenryDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}

	// Calling the base class CRichEditDoc enables serialization
	//  of the container document's COleClientItem objects.
	CRichEditDoc::Serialize(ar);
}

/////////////////////////////////////////////////////////////////////////////
// CFastHenryDoc diagnostics

#ifdef _DEBUG
void CFastHenryDoc::AssertValid() const
{
	CRichEditDoc::AssertValid();
}

void CFastHenryDoc::Dump(CDumpContext& dc) const
{
	CRichEditDoc::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFastHenryDoc commands

// overridden to inhibit document to prompt user for saving changes
BOOL CFastHenryDoc::SaveModified() 
{
	// returning TRUE means that the document can be closed 
	return TRUE;
}

// AutoWrap helper function to implement a client Automation interface
//
// It is used here to call callback functions by our Automation server
// 
// For information about callbacks, see online help under:
//
// - HOWTO: Automate Excel From VC++ Without Using MFC, KB216686
// - b2c.exe VisualBasic to VC++ Automation client converter by Microsoft (web site)
// - Building OLE Automation Service Components in Visual Basic and Visual C++
//   (Hotel Manager Example) and related example project in VBasic
// - Locating Resources To Study Automation, KB152023
// - Office Automation Using VC++, KB196776
// - Platform SDK -> COM and ActiveX -> Automation -> Overview of Automation ->
//   How Do Clients and Objects Interact?
// - Platform SDK -> COM and ActiveX -> Automation -> Accessing ActiveX Objects ->
//   Creating Applications and Tools That Access Objects -> 
//   Accessing Members Through IDispatch / Accessing Members Through VTBLs
// - Inside OLE 2dn Edition
HRESULT CFastHenryDoc::AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...)
{
    // Begin variable-argument list...
    va_list marker;
    va_start(marker, cArgs);

    if(!pDisp) {
        MessageBox(NULL, "NULL IDispatch passed to AutoWrap()", "Error", 0x10010);
        _exit(0);
    }

    // Variables used...
    DISPPARAMS dp = { NULL, NULL, 0, 0 };
    DISPID dispidNamed = DISPID_PROPERTYPUT;
    DISPID dispID;
    HRESULT hr;
    char buf[200];
    char szName[200];
    
    // Convert down to ANSI
    WideCharToMultiByte(CP_ACP, 0, ptName, -1, szName, 256, NULL, NULL);
    
    // Get DISPID for name passed...
    hr = pDisp->GetIDsOfNames(IID_NULL, &ptName, 1, LOCALE_USER_DEFAULT, &dispID);
    if(FAILED(hr)) {
        sprintf(buf, "IDispatch::GetIDsOfNames(\"%s\") failed w/err 0x%08lx", szName, hr);
        MessageBox(NULL, buf, "AutoWrap()", 0x10010);
        _exit(0);
        return hr;
    }
    
    // Allocate memory for arguments...
    VARIANT *pArgs = new VARIANT[cArgs+1];
    // Extract arguments...
    for(int i=0; i<cArgs; i++) {
        pArgs[i] = va_arg(marker, VARIANT);
    }
    
    // Build DISPPARAMS
    dp.cArgs = cArgs;
    dp.rgvarg = pArgs;
    
    // Handle special-case for property-puts!
    if(autoType & DISPATCH_PROPERTYPUT) {
        dp.cNamedArgs = 1;
        dp.rgdispidNamedArgs = &dispidNamed;
    }
    
    // Make the call!
    hr = pDisp->Invoke(dispID, IID_NULL, LOCALE_SYSTEM_DEFAULT, autoType, &dp, pvResult, NULL, NULL);
    if(FAILED(hr)) {
        sprintf(buf, "IDispatch::Invoke(\"%s\"=%08lx) failed w/err 0x%08lx", szName, dispID, hr);
        MessageBox(NULL, buf, "AutoWrap()", 0x10010);
        _exit(0);
        return hr;
    }
    // End variable-argument section...
    va_end(marker);
    
    delete [] pArgs;
    
    return hr;
}

BOOL CFastHenryDoc::IsRunning() 
{
	return g_bIsFCRunning;
}

void CFastHenryDoc::ShowWindow() 
{
	POSITION pos;
	CView* pView;
	CFrameWnd* pFrameWnd;

	pos = GetFirstViewPosition();
	pView = GetNextView(pos);
	if (pView != NULL)
	{
		pFrameWnd = pView->GetParentFrame();
		pFrameWnd->ActivateFrame(SW_SHOW);
		pFrameWnd = pFrameWnd->GetParentFrame();
		if (pFrameWnd != NULL)
			pFrameWnd->ActivateFrame(SW_SHOW);
	}  
	
	// Save the handler of the main window in global variable,
	// so working threads are able to send back messages
    // (threads are not allowed to access any other structure
	// of the calling process except via its window handler,
	// or an ASSERT will fail)
	//
	// No more used here; now assignment is done in CMainFrame::OnCreate(),
	// see comments in CFastHenryApp::InitInstance()
	//FHHwnd = AfxGetApp()->m_pMainWnd->GetSafeHwnd();
	//ASSERT(FHHwnd);

	// perform view initial update
	pView->OnInitialUpdate();
}

BOOL CFastHenryDoc::Quit() 
{
	CWnd *pWnd;

	// cannot quit while FastHenry is running, must stop
	// thread before
	if(g_bIsFCRunning == TRUE)
		return FALSE;

	pWnd = AfxGetApp()->GetMainWnd();

	pWnd->PostMessage(WM_CLOSE);

	return TRUE;
}

VARIANT CFastHenryDoc::GetResistance() 
{
	CFastHenryApp *theApp;
	
	theApp = (CFastHenryApp *)AfxGetApp();

	// build safe-array from capacitance matrix safe array
	COleSafeArray resMatrix(theApp->m_clsResMatrix);

	// return the safe-array encapsulated in a VARIANT
	return resMatrix.Detach();
}

VARIANT CFastHenryDoc::GetInductance() 
{
	CFastHenryApp *theApp;
	
	theApp = (CFastHenryApp *)AfxGetApp();

	// build safe-array from capacitance matrix safe array
	COleSafeArray indMatrix(theApp->m_clsIndMatrix);

	// return the safe-array encapsulated in a VARIANT
	return indMatrix.Detach();
}

VARIANT CFastHenryDoc::GetFrequencies() 
{
	CFastHenryApp *theApp;
	
	theApp = (CFastHenryApp *)AfxGetApp();

	// build safe-array from conductor names safe array
	COleSafeArray nameArray(theApp->m_clsFrequencies);

	// return the safe-array encapsulated in a VARIANT
	return nameArray.Detach();
}

VARIANT CFastHenryDoc::GetRowPortNames() 
{
	CFastHenryApp *theApp;
	
	theApp = (CFastHenryApp *)AfxGetApp();

	// build safe-array from conductor names safe array
	COleSafeArray nameArray(theApp->m_clsRowPortNames);

	// return the safe-array encapsulated in a VARIANT
	return nameArray.Detach();
}

VARIANT CFastHenryDoc::GetColPortNames() 
{
	CFastHenryApp *theApp;
	
	theApp = (CFastHenryApp *)AfxGetApp();

	// build safe-array from conductor names safe array
	COleSafeArray nameArray(theApp->m_clsColPortNames);

	// return the safe-array encapsulated in a VARIANT
	return nameArray.Detach();
}

BOOL CFastHenryDoc::Run(LPCTSTR commandLine) 
{
	CFastHenryApp *theApp;

	theApp = (CFastHenryApp *)AfxGetApp();

	if( g_bIsFCRunning == TRUE) {
		return FALSE;
	} 

	theApp->parseCmdLine(&(theApp->argc), &(theApp->argv), commandLine);

	theApp->LaunchFastHenry();

	return TRUE;
}

void CFastHenryDoc::Stop() 
{
	// signal to the working thread to stop execution
	bFHContinue	= FALSE;
}

BOOL CFastHenryDoc::SetEndCallback(LPDISPATCH callback, LPCTSTR cbName) 
{
	if(callback != NULL) {
		m_pDispEndCallback = callback;
		m_clsEndCallbackName = cbName;
		// tell Automation that there is a new reference to this object
		// so that it does not prematurely release it;
		// otherwise, the reference to the 'callback' object would go
		// out of scope once the SetEndCallback() method is complete,
		// and Automation therefore would dereference and releases the memory. 
		m_pDispEndCallback->AddRef();
		return TRUE;
	}
	else {
		// no callback anymore
		m_pDispEndCallback->Release();
		m_pDispEndCallback = NULL;
		return FALSE;
	}
}

// set an automation callback upon FastHenry working thread end
void CFastHenryDoc::EndCallback()
{
	if(m_pDispEndCallback != NULL) {
		AutoWrap(DISPATCH_METHOD, NULL, m_pDispEndCallback, m_clsEndCallbackName.GetPointer(), 0);
	}
}

BOOL CFastHenryDoc::SetLogCallback(LPDISPATCH callback, LPCTSTR cbName) 
{
	if(callback != NULL) {
		m_pDispLogCallback = callback;
		m_clsLogCallbackName = cbName;
		// tell Automation that there is a new reference to this object
		// so that it does not prematurely release it;
		// otherwise, the reference to the 'callback' object would go
		// out of scope once the SetEndCallback() method is complete,
		// and Automation therefore would dereference and releases the memory. 
		m_pDispLogCallback->AddRef();
		return TRUE;
	}
	else {
		// no callback anymore
		m_pDispLogCallback->Release();
		m_pDispLogCallback = NULL;
		return FALSE;
	}
}

// set an automation callback upon FastHenry working thread log messages (viewprintf()..)
void CFastHenryDoc::LogCallback(CUnicodeString &text, unsigned long color)
{
	VARIANT parm[2];

 	if(m_pDispLogCallback != NULL) {

		parm[0].vt = VT_BSTR;
		parm[0].bstrVal = ::SysAllocString(text.GetPointer());
		parm[1].vt = VT_I4;
		parm[1].lVal = color;
		
		AutoWrap(DISPATCH_METHOD, NULL, m_pDispLogCallback, m_clsLogCallbackName.GetPointer(), 2, parm[1], parm[0]);

		VariantClear(&parm[0]);
		VariantClear(&parm[1]);
	}
}
