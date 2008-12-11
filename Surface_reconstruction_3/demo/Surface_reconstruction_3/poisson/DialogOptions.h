#pragma once

#include "resource.h"


// Class CDialogOptions = Options dialog box
class CDialogOptions : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions)

public:
	CDialogOptions(CWnd* pParent = NULL);
	virtual ~CDialogOptions();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_sm_distance_poisson; // upper bound of distance to surface (Poisson)
	double m_sm_error_bound_poisson; // error bound to stop dichotomy (Poisson)
	double m_sm_radius_poisson; // upper bound of Delaunay balls radii (Poisson)
	double m_sm_angle_poisson; // lower bound of facets angles (degrees) (Poisson)
	double m_sm_distance_apss; // upper bound of distance to surface (APSS)
	double m_sm_error_bound_apss; // error bound to stop dichotomy (APSS)
	double m_sm_radius_apss; // upper bound of Delaunay balls radii (APSS)
	double m_sm_angle_apss; // lower bound of facets angles (degrees) (APSS)
	unsigned int m_dr_max_vertices;
	double m_dr_shell_size;
	double m_dr_sizing;
	double m_contouring_value; // Poisson contouring value (TEST)
	double m_lambda; // laplacian smoothing
	unsigned int m_nb_neighbors_avg_spacing; // K-nearest neighbors (average spacing)
	double m_nb_neighbors_outliers_removal; // K-nearest neighbors (outliers_removal)
	double m_nb_neighbors_smooth_jet_fitting; // K-nearest neighbors (smooth points by Jet Fitting)
	double m_nb_neighbors_pca_normals; // K-nearest neighbors (estimate normals by PCA)
	double m_nb_neighbors_jet_fitting_normals; // K-nearest neighbors (estimate normals by Jet Fitting)
	unsigned int m_nb_neighbors_mst; // K-nearest neighbors (orient normals by MST)
	unsigned int m_nb_neighbors_apss; // K-nearest neighbors (APSS)
	double m_min_cameras_cone_angle; // min angle of camera's cone (degrees) 
    double m_threshold_percent_avg_knn_sq_dst; // percentage of outliers to remove 
    double m_clustering_step; // Grid's step for simplification by clustering 
    double m_random_simplification_percentage; // percentage of random points to remove 
};
