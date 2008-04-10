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
#include <CGAL/Poisson_implicit_function.h>

// This demo
#include "poisson_dt3.h"

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Tetrahedron_3 Tetrahedron;

// Poisson's Delaunay triangulation 3
typedef Poisson_dt3<Kernel> Dt3;
typedef CGAL::Poisson_implicit_function<Kernel, Dt3> Poisson_implicit_function;

// Surface mesher 
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> SVb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel> SCb;
typedef CGAL::Triangulation_data_structure_3<SVb, SCb> STds;
typedef CGAL::Delaunay_triangulation_3<Kernel, STds> STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;


class CPoissonDoc : public CDocument
{
protected: // create from serialization only
	CPoissonDoc();
	DECLARE_DYNCREATE(CPoissonDoc)

// Data members
private:
  // Poisson implicit
  Poisson_implicit_function m_poisson_function; // Poisson implicit function
  Dt3 m_poisson_dt; // The Poisson equation is solved on the vertices of m_poisson_dt
  bool m_triangulation_refined; // Is Delaunay refinement applied?
  bool m_poisson_solved; // Is the Poisson equation solved?

	// Surface mesher 
  STr m_surface_mesher_dt; // 3D-Delaunay triangulation
  C2t3 m_surface_mesher_c2t3; // 2D-complex in m_surface_mesher_dt

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

// Public methods
public:

  /// Get Poisson implicit function.
  Poisson_implicit_function& poisson_function()
  {
    return m_poisson_function;
  }
  const Poisson_implicit_function& poisson_function() const
  {
    return m_poisson_function;
  }

// Private methods
private:

	// misc status stuff
	void update_status();
	void status_message(char* fmt,...);
	double duration(const double time_init);

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
	afx_msg void OnReconstructionDelaunayrefinement();
	afx_msg void OnReconstructionPoisson();
	afx_msg void OnAlgorithmsRefineinshell();
	afx_msg void OnReconstructionSurfacemeshing();
	afx_msg void OnEditOptions();
  afx_msg void OnUpdateReconstructionPoisson(CCmdUI *pCmdUI);
  afx_msg void OnUpdateReconstructionSurfacemeshing(CCmdUI *pCmdUI);
  afx_msg void OnAlgorithmsMarchingtetcontouring();
  afx_msg void OnUpdateAlgorithmsMarchingtetcontouring(CCmdUI *pCmdUI);
  afx_msg void OnFileSaveSurface();
  afx_msg void OnUpdateFileSaveSurface(CCmdUI *pCmdUI);
  afx_msg void OnFileSaveAs();
  afx_msg void OnUpdateFileSaveAs(CCmdUI *pCmdUI);
	afx_msg void OnAlgorithmsExtrapolatenormals();
	afx_msg void OnAlgorithmsPoissonStatistics();
  afx_msg void OnUpdateAlgorithmsPoissonstatistics(CCmdUI *pCmdUI);
  afx_msg void OnAlgorithmsEstimateNormalsByPCA();
  afx_msg void OnAlgorithmsEstimateNormalsByJetFitting();
  afx_msg void OnAlgorithmsOrientNormalsWrtCameras();
  afx_msg void OnAlgorithmsOrientNormalsWithMST();
};


#endif // _DOC_