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
    , m_sm_radius_poisson(0)
    , m_sm_angle_poisson(0)
    , m_sm_distance_apss(0)
    , m_sm_radius_apss(0)
    , m_sm_angle_apss(0)
    , m_nb_neighbors_avg_spacing(0)
    , m_nb_neighbors_remove_outliers(0)
    , m_nb_neighbors_smooth_jet_fitting(0)
    , m_nb_neighbors_pca_normals(0)
    , m_nb_neighbors_jet_fitting_normals(0)
    , m_nb_neighbors_mst(0)
    , m_smoothness_apss(0)
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
    DDX_Text(pDX,IDC_EDIT_SM_ANGLE_POISSON,m_sm_angle_poisson);
    DDX_Text(pDX,IDC_EDIT_SM_RADIUS_POISSON,m_sm_radius_poisson);
    DDX_Text(pDX,IDC_EDIT_SM_DISTANCE_POISSON,m_sm_distance_poisson);
    DDX_Text(pDX,IDC_EDIT_SM_ANGLE_APSS,m_sm_angle_apss);
    DDX_Text(pDX,IDC_EDIT_SM_RADIUS_APSS,m_sm_radius_apss);
    DDX_Text(pDX,IDC_EDIT_SM_DISTANCE_APSS,m_sm_distance_apss);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_AVG_SPACING,m_nb_neighbors_avg_spacing);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_OUTLIER_REMOVAL,m_nb_neighbors_remove_outliers);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_SMOOTH_JET_FITTING,m_nb_neighbors_smooth_jet_fitting);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_PCA_NORMALS,m_nb_neighbors_pca_normals);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_JET_FITTING_NORMALS,m_nb_neighbors_jet_fitting_normals);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_MST,m_nb_neighbors_mst);
    DDX_Text(pDX,IDC_EDIT_NB_NEIGHBORS_APSS,m_smoothness_apss);
    DDX_Text(pDX,IDC_EDIT_CLUSTERING_STEP, m_clustering_step);
    DDX_Text(pDX,IDC_EDIT_RANDOM_SIMPLIFICATION_PERCENTAGE, m_random_simplification_percentage);
}


BEGIN_MESSAGE_MAP(CDialogOptions, CDialog)
END_MESSAGE_MAP()


// Gestionnaires de messages CDialogOptions
