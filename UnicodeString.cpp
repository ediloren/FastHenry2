// UnicodeString.cpp: implementation of the CUnicodeString class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "UnicodeString.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CUnicodeString::CUnicodeString()
{
	m_pString = NULL;
	m_lLen = 0;
}

CUnicodeString::CUnicodeString(const CUnicodeString& string)
{
	m_lLen = string.m_lLen;
	m_pString = new WCHAR[m_lLen];
	memcpy(m_pString, string.m_pString, sizeof(WCHAR)*m_lLen);
}

CUnicodeString::CUnicodeString(LPCTSTR string)
{
	m_pString = NULL;

	*this = string;
}

CUnicodeString::~CUnicodeString()
{
	deleteString();
}

CUnicodeString& CUnicodeString::operator=(LPCTSTR string)
{
	deleteString();

	// calculate wchar length of string
	m_lLen = MultiByteToWideChar(CP_ACP, 0, string, -1, NULL, NULL);
	// allocate memory
	m_pString = new WCHAR[m_lLen];
	// perform actual conversion from ANSI to Unicode
	MultiByteToWideChar(CP_ACP, 0, string, -1, m_pString, m_lLen);

	return *this;
}

LPWSTR CUnicodeString::GetPointer()
{
	return m_pString;
}

void CUnicodeString::deleteString()
{
	if(m_pString != NULL) {
		delete[] m_pString;
		m_pString = NULL;
	}

	m_lLen = 0;
}
