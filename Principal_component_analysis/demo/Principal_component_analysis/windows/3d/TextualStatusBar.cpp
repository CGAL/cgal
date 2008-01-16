// TextualStatusBar.cpp : implementation file
//

#include "stdafx.h"
#include "TextualStatusBar.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CTextualStatusBar

CTextualStatusBar::CTextualStatusBar()
{
}

CTextualStatusBar::~CTextualStatusBar()
{
}

BEGIN_MESSAGE_MAP(CTextualStatusBar, CStatusBar)
	//{{AFX_MSG_MAP(CTextualStatusBar)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CTextualStatusBar message handlers

int CTextualStatusBar::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if(CStatusBar::OnCreate(lpCreateStruct) == -1)
		return -1;

	if(!m_wndToolTip.Create(this, TTS_ALWAYSTIP) ||
		 !m_wndToolTip.AddTool(this, "Tool tip for status bar"))
		AfxMessageBox("Unable to add tool tip for status bar");

	return 0;
}

void CTextualStatusBar::OnDestroy()
{
	m_wndToolTip.DestroyWindow();
	CStatusBar::OnDestroy();
}

BOOL CTextualStatusBar::PreTranslateMessage(MSG* pMsg)
{
	m_wndToolTip.RelayEvent(pMsg);
	return CStatusBar::PreTranslateMessage(pMsg);
}

