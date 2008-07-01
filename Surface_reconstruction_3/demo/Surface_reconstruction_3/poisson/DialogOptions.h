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
	double m_sm_distance; // upper bound of distance to surface
	double m_sm_error_bound; // error bound to stop dichotomy
	double m_sm_radius; // upper bound of Delaunay balls radii
	double m_sm_angle; // lower bound of facets angles (degrees)
	unsigned int m_dr_max_vertices;
	double m_dr_shell_size;
	double m_dr_sizing;
	double m_contouring_value; // Poisson contouring value (TEST)
	double m_lambda;
	double m_projection_error; // APSS projection error
	unsigned int m_number_of_neighbours; // K-nearest neighbours
	double m_min_cameras_cone_angle; // min angle of camera's cone (degrees) 
    int m_threshold_percent_avg_knn_sq_dst; // percentage of outliers to remove 
};
