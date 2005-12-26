// DialogOptions.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Mesh.h"
#include "DialogOptions.h"


// Boîte de dialogue CDialogOptions

IMPLEMENT_DYNAMIC(CDialogOptions, CDialog)
CDialogOptions::CDialogOptions(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions::IDD, pParent)
{
}

CDialogOptions::~CDialogOptions()
{
}

void CDialogOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CDialogOptions, CDialog)
END_MESSAGE_MAP()


// Gestionnaires de messages CDialogOptions
