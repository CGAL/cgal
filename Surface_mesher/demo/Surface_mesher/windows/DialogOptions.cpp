// DialogOptions.cpp : fichier d'implémentation
#include "stdafx.h"
#include "Mesh.h"
#include "DialogOptions.h"


IMPLEMENT_DYNAMIC(CDialogOptions, CDialog)
CDialogOptions::CDialogOptions(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions::IDD, pParent)
	, m_criterion_uniform_size(0)
	, m_max_nb_vertices(0)
	, m_refresh_each(0)
	, m_init_nb_vertices(0)
{
}

CDialogOptions::~CDialogOptions()
{
}

void CDialogOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text(pDX,IDC_EDIT_INIT_NBVERTICES,m_init_nb_vertices);
	DDV_MinMaxUInt(pDX, m_init_nb_vertices, 1,10000000);

	DDX_Text(pDX,IDC_EDIT_UNIFORM_SIZE,m_criterion_uniform_size);
	DDV_MinMaxDouble(pDX, m_criterion_uniform_size, 0.00001,100);

	DDX_Text(pDX,IDC_EDIT_NB_MAX_VERTICES,m_max_nb_vertices);
	DDV_MinMaxUInt(pDX, m_max_nb_vertices, 1,10000000);

	DDX_Text(pDX,IDC_EDIT_REFRESH_EACH,m_refresh_each);
	DDV_MinMaxUInt(pDX, m_refresh_each, 1,10000000);
}


BEGIN_MESSAGE_MAP(CDialogOptions, CDialog)
END_MESSAGE_MAP()
