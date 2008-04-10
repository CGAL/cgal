// DialogOptions_4.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Gyroviz.h"
#include "DialogOptions_4.h"


// Boîte de dialogue CDialogOptions_4

IMPLEMENT_DYNAMIC(CDialogOptions_4, CDialog)
CDialogOptions_4::CDialogOptions_4(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions_4::IDD, pParent)
	, m_point_size(2)
{
}

CDialogOptions_4::~CDialogOptions_4()
{
}

void CDialogOptions_4::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_EDIT_POINT_SIZE,m_point_size);
}


BEGIN_MESSAGE_MAP(CDialogOptions_4, CDialog)
END_MESSAGE_MAP()
