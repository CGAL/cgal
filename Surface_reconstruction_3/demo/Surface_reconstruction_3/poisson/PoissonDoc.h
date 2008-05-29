// PoissonDoc.h : interface of the CPoissonDoc class
//

#ifndef _DOC_
#define _DOC_
#pragma once

// CGAL
#include <CGAL/basic.h>
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
#include "Gyroviz_point_3.h"
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

// Point + normal
typedef Gyroviz_point_3<Kernel> Point_with_normal; // Model of the PointWithNormal_3 concept
typedef Point_with_normal::Normal Normal; // Model of OrientedNormal_3 concept

// Point set
typedef Point_set_3<Kernel> Point_set;

// Poisson's Delaunay triangulation 3 and implicit function
typedef Poisson_dt3<Kernel> Dt3;
typedef CGAL::Poisson_implicit_function<Kernel, Dt3> Poisson_implicit_function;

// APSS implicit function
typedef CGAL::APSS_implicit_function<Kernel,Point_with_normal> APSS_implicit_function;

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

    // Surface mesher options
    double m_sm_angle;
    double m_sm_radius;
    double m_sm_distance;

    // Delaunay refinement options
    double m_dr_sizing;
    double m_dr_shell_size;
    unsigned int m_dr_max_vertices;

    // Surface mesher and marching tet common options
    double m_contouring_value;

    // Normal estimation options
    unsigned int m_number_of_neighbours;

    // Outlier removal
    double m_outlier_percentage;

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
    // Utility: compute elapsed time
    double duration(const double time_init);

    // Clean up current edit mode
    void CloseMode();

// MFC generated
public:
	virtual ~CPoissonDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
    virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
    afx_msg void OnReconstructionDelaunayRefinement();
    afx_msg void OnReconstructionPoisson();
    afx_msg void OnAlgorithmsRefineInShell();
    afx_msg void OnReconstructionPoissonSurfaceMeshing();
    afx_msg void OnEditOptions();
    afx_msg void OnUpdateReconstructionPoisson(CCmdUI *pCmdUI);
    afx_msg void OnUpdateReconstructionPoissonSurfaceMeshing(CCmdUI *pCmdUI);
    afx_msg void OnAlgorithmsMarchingTetContouring();
    afx_msg void OnUpdateAlgorithmsMarchingTetContouring(CCmdUI *pCmdUI);
    afx_msg void OnFileSaveSurface();
    afx_msg void OnUpdateFileSaveSurface(CCmdUI *pCmdUI);
    afx_msg void OnFileSaveAs();
    afx_msg void OnUpdateFileSaveAs(CCmdUI *pCmdUI);
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
    afx_msg void OnRemoveOutliers();
    afx_msg void OnUpdateRemoveOutliers(CCmdUI *pCmdUI);
    afx_msg void OnAnalysisAverageSpacing();
    afx_msg void OnOneStepPoissonReconstruction();
    afx_msg void OnUpdateOneStepPoissonReconstruction(CCmdUI *pCmdUI);
    afx_msg void OnUpdateAnalysisAverageSpacing(CCmdUI *pCmdUI);
    afx_msg void OnReconstructionApssReconstruction();
    afx_msg void OnUpdateReconstructionApssReconstruction(CCmdUI *pCmdUI);
    afx_msg void OnModeAPSS();
    afx_msg void OnUpdateModeAPSS(CCmdUI *pCmdUI);
};


#endif // _DOC_