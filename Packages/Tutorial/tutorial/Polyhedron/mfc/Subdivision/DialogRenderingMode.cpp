// DialogRenderingMode.cpp : implementation file
//

#include "stdafx.h"
#include "Mesh.h"
#include "DialogRenderingMode.h"
#include ".\dialogrenderingmode.h"


// CDialogRenderingMode dialog

IMPLEMENT_DYNAMIC(CDialogRenderingMode, CDialog)
CDialogRenderingMode::CDialogRenderingMode(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogRenderingMode::IDD, pParent)
{
	m_Mode = 0;
}

CDialogRenderingMode::~CDialogRenderingMode()
{
}

void CDialogRenderingMode::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CDialogRenderingMode, CDialog)
	ON_WM_LBUTTONUP()
	ON_WM_CLOSE()
	ON_WM_DESTROY()
END_MESSAGE_MAP()


// CDialogRenderingMode message handlers

void CDialogRenderingMode::OnLButtonUp(UINT nFlags, 
																			 CPoint point)
{
	int x = (int)(point.x/200);
	int y = (int)(point.y/150);
	m_Mode = 3*y+x;
	this->PostMessage(WM_CLOSE);
}

void CDialogRenderingMode::OnClose()
{
	CDialog::OnClose();
}

void CDialogRenderingMode::OnDestroy()
{
	CDialog::OnDestroy();
}
