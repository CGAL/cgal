// DialogOptions_2.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Gyroviz.h"
#include "DialogOptions_2.h"


// Boîte de dialogue CDialogOptions_2

IMPLEMENT_DYNAMIC(CDialogOptions_2, CDialog)
CDialogOptions_2::CDialogOptions_2(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions_2::IDD, pParent)
	, m_point_size(2)
{
}

CDialogOptions_2::~CDialogOptions_2()
{
}

void CDialogOptions_2::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_EDIT_POINT_SIZE,m_point_size);
}


BEGIN_MESSAGE_MAP(CDialogOptions_2, CDialog)
END_MESSAGE_MAP()


