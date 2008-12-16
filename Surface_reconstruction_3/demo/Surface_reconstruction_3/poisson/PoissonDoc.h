// PoissonDoc.h : interface of the CPoissonDoc class
//

#ifndef _DOC_
#define _DOC_
#pragma once

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>

// This package
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Poisson_implicit_function.h>
#include <CGAL/APSS_implicit_function.h>

// This demo
#include "poisson_dt3.h"
#include "UI_point_3.h"
#include "Point_set_3.h"
#include "Triangular_surface_3.h"

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Tetrahedron_3 Tetrahedron;

// Point set
typedef Point_set_3<Kernel> Point_set;

/// Type of points in Point_set_3
typedef UI_point_3<Kernel> UI_point; ///< Position + normal + cameras + selection flag
// Its superclasses:
typedef Gyroviz_point_3<Kernel> Gyroviz_point; ///< Position + normal + cameras
typedef UI_point::Point_with_normal Point_with_normal; ///< Position + normal

// Type of points normal
typedef UI_point::Normal Normal; ///< Model of OrientableNormal_3 concept.

// Poisson's 3D Delaunay triangulation and implicit function
typedef Poisson_dt3<Kernel> Dt3;
typedef CGAL::Poisson_implicit_function<Kernel, Dt3> Poisson_implicit_function;

// APSS implicit function
typedef CGAL::APSS_implicit_function<Kernel> APSS_implicit_function;

// Surface mesh generator 
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> SVb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel> SCb;
typedef CGAL::Triangulation_data_structure_3<SVb, SCb> STds;
typedef CGAL::Delaunay_triangulation_3<Kernel, STds> STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;

// Reconstructed surface
typedef Triangular_surface_3<Kernel> Triangular_surface;


// MFC document.
//
// This class owns the edited points. They exist under 2 forms:
// - m_points[] array of points + normals.
// - the m_poisson_dt 3D triangulation used by the m_poisson_function implicit function
//   (in fact, the points in m_poisson_dt are a non editable copy of m_points[]).
//
// Only 1 form is visible on screen and editable at a given time. 
// This is controlled by m_edit_mode.
//
class CPoissonDoc : public CDocument
{
protected: // create from serialization only
	CPoissonDoc();
	DECLARE_DYNCREATE(CPoissonDoc)

// Public types
public:

    // Available edit modes
    enum Edit_mode { NO_EDIT_MODE, POINT_SET, POISSON, APSS };

// Data members
private:

    // Current edit mode
    Edit_mode m_edit_mode;

    // Input point set
    Point_set m_points;

    // Poisson implicit function
    Poisson_implicit_function* m_poisson_function; // Poisson implicit function
    Dt3* m_poisson_dt; // The Poisson equation is solved on the vertices of m_poisson_dt
    bool m_triangulation_refined; // Is Delaunay refinement applied?
    bool m_poisson_solved; // Is the Poisson equation solved?

    // - APSS implicit function
    APSS_implicit_function* m_apss_function;

    // Surface mesher 
    STr m_surface_mesher_dt; // 3D-Delaunay triangulation
    C2t3 m_surface_mesher_c2t3; // 2D-complex in m_surface_mesher_dt
    Triangular_surface m_surface; // Surface reconstructed by Surface mesher

    // Surface reconstructed by marching tet contouring
    Triangular_surface m_contour;

    // Poisson options
	double m_sm_distance_poisson; // upper bound of distance to surface (Poisson)
	double m_sm_error_bound_poisson; // error bound to stop dichotomy (Poisson)
	double m_sm_radius_poisson; // upper bound of Delaunay balls radii (Poisson)
	double m_sm_angle_poisson; // lower bound of facets angles (degrees) (Poisson)
    double m_dr_sizing; // 3 Delaunay refinements options
    double m_dr_shell_size;
    unsigned int m_dr_max_vertices;
    double m_contouring_value; // Poisson contouring value (TEST)
    double m_lambda;  // laplacian smoothing

    // APSS options
	double m_sm_distance_apss; // upper bound of distance to surface (APSS)
	double m_sm_error_bound_apss; // error bound to stop dichotomy (APSS)
	double m_sm_radius_apss; // upper bound of Delaunay balls radii (APSS)
	double m_sm_angle_apss; // lower bound of facets angles (degrees) (APSS)
	unsigned int m_nb_neighbors_apss; // K-nearest neighbors (APSS)

    // Average Spacing options
	unsigned int m_nb_neighbors_avg_spacing; // K-nearest neighbors (average spacing)

    // Smoothing options
	double m_nb_neighbors_smooth_jet_fitting; // K-nearest neighbors (smooth points by Jet Fitting)

    // Normals Computing options
	double m_nb_neighbors_pca_normals; // K-nearest neighbors (estimate normals by PCA)
	double m_nb_neighbors_jet_fitting_normals; // K-nearest neighbors (estimate normals by Jet Fitting)
	unsigned int m_nb_neighbors_mst; // K-nearest neighbors (orient normals by MST)

    // Outlier removal
    double m_min_cameras_cone_angle; // min angle of camera's cone (degrees)
    double m_threshold_percent_avg_knn_sq_dst; // percentage of outliers to remove
	double m_nb_neighbors_outliers_removal; // K-nearest neighbors (outliers_removal)

    // Point set simplification
    double m_clustering_step; // Grid's step for simplification by clustering 
    double m_random_simplification_percentage; // percentage of random points to remove

// Public methods
public:

    // Get input point set. It is exported as const objects as
    // the rules to modify it are complex. See comment above.
    //
    // - as array of points + normals
    const Point_set* points() const { return &m_points; }
    //
    // - as Poisson implicit function
    const Poisson_implicit_function* poisson_function() const 
    { return m_poisson_function; }
    //
    // - as APSS implicit function
    const APSS_implicit_function* apss_function() const 
    { return m_apss_function; }
    //
    // The current edit mode indicates which form is valid.
    Edit_mode edit_mode() const { return m_edit_mode; }

    // Surface reconstructed by Surface mesher
    const Triangular_surface* surface_mesher_surface() const { return &m_surface; }

    // Surface reconstructed by marching tet contouring
    const Triangular_surface* marching_tet_countour() const { return &m_contour; }

// Private methods
private:

    // Update the number of vertices and tetrahedra in the status bar
    // and write them to cerr.
    void update_status();
    // Write user message in status bar and cerr	
    void status_message(char* fmt,...);
    // Write user message in message box and cerr
    void prompt_message(char* fmt,...);

    // Clean up current edit mode
    void CloseMode();

    // Check the accuracy of normals direction estimation.
    // If original normals are available, compare with them and select normals with large deviation.
    // @return true on success.
    bool verify_normals_direction();
    // Check the accuracy of normals orientation.
    // Count and select non-oriented normals.
    // If original normals are available, compare with them and select flipped normals.
    // @return true on success.
    bool verify_normals_orientation();

// MFC generated
public:
	virtual ~CPoissonDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
    virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
    afx_msg void OnFileSaveSurface();
    afx_msg void OnUpdateFileSaveSurface(CCmdUI *pCmdUI);
    afx_msg void OnFileSaveAs();
    afx_msg void OnUpdateFileSaveAs(CCmdUI *pCmdUI);
    afx_msg void OnEditOptions();
    afx_msg void OnEditDelete();
    afx_msg void OnUpdateEditDelete(CCmdUI *pCmdUI);
    afx_msg void OnEditResetSelection();
    afx_msg void OnUpdateEditResetSelection(CCmdUI *pCmdUI);
	afx_msg void OnOneStepPoissonReconstructionWithNormalizedDivergence();
    afx_msg void OnReconstructionDelaunayRefinement();
    afx_msg void OnReconstructionPoisson();
	afx_msg void OnReconstructionPoissonNormalized();
    afx_msg void OnAlgorithmsRefineInShell();
    afx_msg void OnReconstructionPoissonSurfaceMeshing();
    afx_msg void OnUpdateReconstructionPoisson(CCmdUI *pCmdUI);
    afx_msg void OnUpdateReconstructionPoissonSurfaceMeshing(CCmdUI *pCmdUI);
    afx_msg void OnAlgorithmsMarchingTetContouring();
    afx_msg void OnUpdateAlgorithmsMarchingTetContouring(CCmdUI *pCmdUI);
    afx_msg void OnAlgorithmsExtrapolatenormals();
    afx_msg void OnAlgorithmsPoissonStatistics();
    afx_msg void OnUpdateAlgorithmsPoissonStatistics(CCmdUI *pCmdUI);
    afx_msg void OnAlgorithmsEstimateNormalsByPCA();
    afx_msg void OnAlgorithmsEstimateNormalsByJetFitting();
    afx_msg void OnAlgorithmsOrientNormalsWrtCameras();
    afx_msg void OnAlgorithmsOrientNormalsWithMST();
    afx_msg void OnAlgorithmsSmoothUsingJetFitting();
    afx_msg void OnModePointSet();
    afx_msg void OnUpdateModePointSet(CCmdUI *pCmdUI);
    afx_msg void OnModePoisson();
    afx_msg void OnUpdateModePoisson(CCmdUI *pCmdUI);
    afx_msg void OnCreatePoissonTriangulation();
    afx_msg void OnUpdateCreatePoissonTriangulation(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsSmoothUsingJetFitting(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsEstimateNormalsByPCA(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsEstimateNormalByJetFitting(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsOrientNormalsWithMST(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsOrientNormalsWrtCameras(CCmdUI *pCmdUI);
    afx_msg void OnUpdateReconstructionDelaunayRefinement(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsRefineInShell(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAlgorithmsExtrapolateNormals(CCmdUI *pCmdUI);
    afx_msg void OnAlgorithmsOutliersRemovalWrtCamerasConeAngle();
    afx_msg void OnUpdateAlgorithmsOutliersRemovalWrtCamerasConeAngle(CCmdUI *pCmdUI);
    afx_msg void OnOutliersRemovalWrtAvgKnnSqDist();
    afx_msg void OnUpdateOutliersRemovalWrtAvgKnnSqDist(CCmdUI *pCmdUI);
    afx_msg void OnAnalysisAverageSpacing();
    afx_msg void OnOneStepPoissonReconstruction();
    afx_msg void OnUpdateOneStepPoissonReconstruction(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAnalysisAverageSpacing(CCmdUI *pCmdUI);
    afx_msg void OnReconstructionApssReconstruction();
    afx_msg void OnUpdateReconstructionApssReconstruction(CCmdUI *pCmdUI);
    afx_msg void OnModeAPSS();
    afx_msg void OnUpdateModeAPSS(CCmdUI *pCmdUI);
	afx_msg void OnCalculateAverageSpacing();
	afx_msg void OnExtrapolateNormalsUsingGaussianKernel();
	afx_msg void OnReconstructionSaveas();
    afx_msg void OnPointCloudSimplificationByClustering();
    afx_msg void OnPointCloudSimplificationRandom();
    afx_msg void OnUpdatePointCloudSimplificationByClustering(CCmdUI *pCmdUI);
    afx_msg void OnUpdatePointCloudSimplificationRandom(CCmdUI *pCmdUI);
    afx_msg void OnRadialNormalsOrientation();
    afx_msg void OnUpdateRadialNormalsOrientation(CCmdUI *pCmdUI);
    afx_msg void OnFlipNormals();
    afx_msg void OnUpdateFlipNormals(CCmdUI *pCmdUI);
    afx_msg void OnUpdateStepExtrapolatenormalsusinggaussiankernels(CCmdUI *pCmdUI);
    afx_msg void OnUpdateReconstructionPoissonNormalized(CCmdUI *pCmdUI);
    afx_msg void OnUpdateReconstructionOneStepNormalized(CCmdUI *pCmdUI);
    afx_msg void OnUpdateStepCalculateaveragespacing(CCmdUI *pCmdUI);
};


#endif // _DOC_