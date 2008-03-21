// DialogOptions_3.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Gyroviz.h"
#include "DialogOptions_3.h"


// Boîte de dialogue CDialogOptions_3

IMPLEMENT_DYNAMIC(CDialogOptions_3, CDialog)
CDialogOptions_3::CDialogOptions_3(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions_3::IDD, pParent)
	, m_point_size(2)
{
}

CDialogOptions_3::~CDialogOptions_3()
{
}

void CDialogOptions_3::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_EDIT_POINT_SIZE,m_point_size);
}


BEGIN_MESSAGE_MAP(CDialogOptions_3, CDialog)
END_MESSAGE_MAP()


