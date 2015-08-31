// Copyright (c) 2013 FastFieldSolvers S.r.l.
// http://www.fastfieldsolvers.com
// All Rights reserved
//
// Usage is subject to the License that you should have received
// as a text file together with these source files. The License text is also available
// in the relevant source code distribution from http://www.fastfieldsolvers.com
// or contacting FastFieldSolvers S.r.l.


// CntrItem.cpp : implementation of the CFastHenry2CntrItem class
//


#include "stdafx.h"
#include "FastHenry2.h"

#include "FastHenryDoc.h"
#include "FastHenryView.h"
#include "CntrItem.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFastHenry2CntrItem implementation

IMPLEMENT_SERIAL(CFastHenry2CntrItem, CRichEditCntrItem, 0)

CFastHenry2CntrItem::CFastHenry2CntrItem(REOBJECT* preo, CFastHenryDoc* pContainer)
	: CRichEditCntrItem(preo, pContainer)
{
	// TODO: add one-time construction code here
	
}

CFastHenry2CntrItem::~CFastHenry2CntrItem()
{
	// TODO: add cleanup code here
	
}

/////////////////////////////////////////////////////////////////////////////
// CFastHenry2CntrItem diagnostics

#ifdef _DEBUG
void CFastHenry2CntrItem::AssertValid() const
{
	CRichEditCntrItem::AssertValid();
}

void CFastHenry2CntrItem::Dump(CDumpContext& dc) const
{
	CRichEditCntrItem::Dump(dc);
}
#endif

/////////////////////////////////////////////////////////////////////////////
