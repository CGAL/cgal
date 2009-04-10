#pragma once

#include "resource.h"


// Class CDialogOptions = Options dialog box
class CDialogOptions : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions)

public:
	CDialogOptions(CWnd* pParent = NULL);
	virtual ~CDialogOptions();

	enum { IDD = IDD_DIALOG_OPTIONS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_sm_distance_poisson; // Approximation error w.r.t. point set radius (Poisson)
	double m_sm_radius_poisson; // Max triangle radius w.r.t. point set radius (Poisson)
	double m_sm_angle_poisson; // Min triangle angle (degrees) (Poisson)
	double m_sm_distance_apss; // Approximation error w.r.t. point set radius (APSS)
	double m_sm_radius_apss; // Max triangle radius w.r.t. point set radius (APSS)
	double m_sm_angle_apss; // Min triangle angle (degrees) (APSS)
	unsigned int m_dr_max_vertices;
	double m_dr_shell_size;
	double m_dr_sizing;
	double m_contouring_value; // Poisson contouring value (TEST)
	double m_lambda; // laplacian smoothing
	unsigned int m_nb_neighbors_avg_spacing; // K-nearest neighbors (average spacing)
	double m_nb_neighbors_outlier_removal; // K-nearest neighbors (outlier_removal)
	double m_nb_neighbors_smooth_jet_fitting; // K-nearest neighbors (smooth points by Jet Fitting)
	double m_nb_neighbors_pca_normals; // K-nearest neighbors (estimate normals by PCA)
	double m_nb_neighbors_jet_fitting_normals; // K-nearest neighbors (estimate normals by Jet Fitting)
	unsigned int m_nb_neighbors_mst; // K-nearest neighbors (orient normals by MST)
	unsigned int m_nb_neighbors_apss; // #neighbors to compute APPS sphere fitting
	double m_min_cameras_cone_angle; // min angle of camera's cone (degrees)
    double m_threshold_percent_avg_knn_sq_dst; // percentage of outliers to remove
    double m_clustering_step; // Grid's step for simplification by clustering
    double m_random_simplification_percentage; // percentage of random points to remove
};
