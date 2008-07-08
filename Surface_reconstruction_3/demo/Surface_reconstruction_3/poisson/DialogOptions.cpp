// DialogOptions.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "Poisson.h"
#include "DialogOptions.h"


// Boîte de dialogue CDialogOptions

IMPLEMENT_DYNAMIC(CDialogOptions, CDialog)
CDialogOptions::CDialogOptions(CWnd* pParent /*=NULL*/)
    : CDialog(CDialogOptions::IDD, pParent)
    , m_sm_distance_poisson(0)
    , m_sm_error_bound(0)
    , m_sm_radius(0)
    , m_sm_angle(0)
    , m_dr_max_vertices(0)
    , m_dr_shell_size(0)
    , m_dr_sizing(0)
    , m_contouring_value(0)
    , m_lambda(0)
    , m_sm_distance_apss(0)
    , m_number_of_neighbours(0)
    , m_threshold_percent_avg_knn_sq_dst(0)
    , m_clustering_step(0)
    , m_random_simplification_percentage(0)
{
}

CDialogOptions::~CDialogOptions()
{
}

void CDialogOptions::DoDataExchange(CDataExchange* pDX)
{
    CDialog::DoDataExchange(pDX);

    // processing
    DDX_Text(pDX,IDC_EDIT_MIN_CAMERAS_CONE_ANGLE, m_min_cameras_cone_angle);
    DDX_Text(pDX,IDC_EDIT_AVG_KNN_SQ_DST_PERCENTAGE, m_threshold_percent_avg_knn_sq_dst);
    DDX_Text(pDX,IDC_EDIT_SM_ANGLE,m_sm_angle);
    DDX_Text(pDX,IDC_EDIT_SM_RADIUS,m_sm_radius);
    DDX_Text(pDX,IDC_EDIT_SM_DISTANCE_POISSON,m_sm_distance_poisson);
    DDX_Text(pDX,IDC_EDIT_SM_ERROR_BOUND,m_sm_error_bound);
    DDX_Text(pDX,IDC_EDIT_DR_SHELL_SIZE,m_dr_shell_size);
    DDX_Text(pDX,IDC_EDIT_DR_SIZING,m_dr_sizing);
    DDX_Text(pDX,IDC_EDIT_DR_MAXV,m_dr_max_vertices);
    DDX_Text(pDX,IDC_EDIT_CONTOURING_VALUE,m_contouring_value);
    DDX_Text(pDX,IDC_EDIT_LAMBDA,m_lambda);
    DDX_Text(pDX,IDC_EDIT_SM_DISTANCE_APSS,m_sm_distance_apss);
    DDX_Text(pDX,IDC_EDIT_NB_OF_NEIGHBOURS,m_number_of_neighbours);
    DDX_Text(pDX,IDC_EDIT_CLUSTERING_STEP, m_clustering_step);
    DDX_Text(pDX,IDC_EDIT_RANDOM_SIMPLIFICATION_PERCENTAGE, m_random_simplification_percentage);
}


BEGIN_MESSAGE_MAP(CDialogOptions, CDialog)
END_MESSAGE_MAP()


// Gestionnaires de messages CDialogOptions
