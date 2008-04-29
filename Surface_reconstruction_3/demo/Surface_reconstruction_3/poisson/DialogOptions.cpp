// DialogOptions.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Poisson.h"
#include "DialogOptions.h"


// Boîte de dialogue CDialogOptions

IMPLEMENT_DYNAMIC(CDialogOptions, CDialog)
CDialogOptions::CDialogOptions(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogOptions::IDD, pParent)
	, m_sm_distance(0)
	, m_sm_radius(0)
	, m_sm_angle(0)
	, m_dr_max_vertices(0)
	, m_dr_shell_size(0)
	, m_dr_sizing(0)
	, m_contouring_value(0)
  , m_number_of_neighbours(0)
  , m_outlier_percentage(0)
{
}

CDialogOptions::~CDialogOptions()
{
}

void CDialogOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

  // processing
	DDX_Text(pDX,IDC_EDIT_OUTLIER_PERCENTAGE,m_outlier_percentage);

	DDX_Text(pDX,IDC_EDIT_SM_ANGLE,m_sm_angle);
	DDX_Text(pDX,IDC_EDIT_SM_RADIUS,m_sm_radius);
	DDX_Text(pDX,IDC_EDIT_SM_DISTANCE,m_sm_distance);

	DDX_Text(pDX,IDC_EDIT_DR_SHELL_SIZE,m_dr_shell_size);
	DDX_Text(pDX,IDC_EDIT_DR_SIZING,m_dr_sizing);
	DDX_Text(pDX,IDC_EDIT_DR_MAXV,m_dr_max_vertices);

	DDX_Text(pDX,IDC_EDIT_CONTOURING_VALUE,m_contouring_value);

	DDX_Text(pDX,IDC_EDIT_NB_OF_NEIGHBOURS,m_number_of_neighbours);
}


BEGIN_MESSAGE_MAP(CDialogOptions, CDialog)
END_MESSAGE_MAP()


// Gestionnaires de messages CDialogOptions
