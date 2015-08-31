// UnicodeString.h: interface for the CUnicodeString class
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UNICODESTRING_H__DD097111_A6A2_471E_9FFC_2F3AF397CD95__INCLUDED_)
#define AFX_UNICODESTRING_H__DD097111_A6A2_471E_9FFC_2F3AF397CD95__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CUnicodeString  
{
public:
	CUnicodeString();
	// copy constructor
	CUnicodeString(const CUnicodeString& string);
	CUnicodeString(LPCTSTR string);
	~CUnicodeString();
	CUnicodeString& operator=(LPCTSTR string);
	LPWSTR GetPointer();

protected:
	void deleteString();

	LPWSTR m_pString;
	long m_lLen;
};

#endif // !defined(AFX_UNICODESTRING_H__DD097111_A6A2_471E_9FFC_2F3AF397CD95__INCLUDED_)
