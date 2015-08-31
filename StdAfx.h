// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//


#if !defined(AFX_STDAFX_H__E1385E87_3037_11D5_9282_74D314C10000__INCLUDED_)
#define AFX_STDAFX_H__E1385E87_3037_11D5_9282_74D314C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

// to avoid the warning "_WIN32_WINNT not defined. Defaulting to _WIN32_WINNT_MAXVER (see WinSDKVer.h)"
// The list of valid values for _WIN32_WINNT can be found at https://msdn.microsoft.com/en-us/library/aa383745%28VS.85%29.aspx
#define _WIN32_WINNT _WIN32_WINNT_WINXP

#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers

// disable deprecation of standard C functions like strcpy(), where the compiler
// suggest to use strcpy_s() instead
#define _CRT_SECURE_NO_WARNINGS

#include <afxwin.h>         // MFC core and standard components
#include <afxext.h>         // MFC extensions
#include <afxole.h>         // MFC OLE classes
#include <afxodlgs.h>       // MFC OLE dialog classes
#include <afxdisp.h>        // MFC OLE automation classes
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			// MFC support for Windows Common Controls
#endif // _AFX_NO_AFXCMN_SUPPORT

#include <afxrich.h>		// MFC rich edit classes

// Html help
#include <htmlhelp.h>

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__E1385E87_3037_11D5_9282_74D314C10000__INCLUDED_)
