// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// CntrItem.h : interface of the CFastHenry2CntrItem class
//


#if !defined(AFX_CNTRITEM_H__E1385E90_3037_11D5_9282_74D314C10000__INCLUDED_)
#define AFX_CNTRITEM_H__E1385E90_3037_11D5_9282_74D314C10000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

class CFastHenryDoc;
class CFastHenryView;

class CFastHenry2CntrItem : public CRichEditCntrItem
{
	DECLARE_SERIAL(CFastHenry2CntrItem)

// Constructors
public:
	CFastHenry2CntrItem(REOBJECT* preo = NULL, CFastHenryDoc* pContainer = NULL);
		// Note: pContainer is allowed to be NULL to enable IMPLEMENT_SERIALIZE.
		//  IMPLEMENT_SERIALIZE requires the class have a constructor with
		//  zero arguments.  Normally, OLE items are constructed with a
		//  non-NULL document pointer.

// Attributes
public:
	CFastHenryDoc* GetDocument()
		{ return (CFastHenryDoc*)CRichEditCntrItem::GetDocument(); }
	CFastHenryView* GetActiveView()
		{ return (CFastHenryView*)CRichEditCntrItem::GetActiveView(); }

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFastHenry2CntrItem)
	public:
	protected:
	//}}AFX_VIRTUAL

// Implementation
public:
	~CFastHenry2CntrItem();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_CNTRITEM_H__E1385E90_3037_11D5_9282_74D314C10000__INCLUDED_)
