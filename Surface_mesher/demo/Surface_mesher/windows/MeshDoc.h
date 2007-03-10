// MeshDoc.h : interface of the CMeshDoc class
#pragma once


#include <CGAL/basic.h>
#include "lib/dt3.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges_without_boundary.h>
#include <CGAL/Surface_mesher/Surface_mesher_manifold.h>
#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Surface_mesher/Polyhedral.h>

#include <fstream>

struct K2 : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Robust_circumcenter_traits_3<K2>  K;
typedef K::Point_3 Point;
typedef K::FT FT;

typedef CGAL::Triangulation_vertex_base_3<K> Vb;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<K, Vb> Vb2;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb2, Cb> Tds;
typedef DT3<K, Tds> Del;
// typedef CGAL::Delaunay_triangulation_3<K, Tds> Del;
typedef CGAL::Complex_2_in_triangulation_3_surface_mesh<Del> C2t3;

// Oracle
typedef CGAL::Surface_mesher::Polyhedral <Del> Oracle;
typedef CGAL::Surface_mesher::Refine_criterion<Del> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;

typedef CGAL::Surface_mesher::Surface_mesher<Del, Oracle, Criteria> SM;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges<Del, Oracle, Criteria> SMRE;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary<Del, Oracle, Criteria> SMREWB;
typedef CGAL::Surface_mesher::Surface_mesher_manifold<Del, Oracle, Criteria> SMM;
typedef CGAL::Surface_mesher::Surface_mesher_regular_edges_without_boundary_base<Del, Oracle, Criteria> SMREWBB;
typedef CGAL::Surface_mesher::Surface_mesher_manifold<Del, Oracle, Criteria,
          CGAL::Surface_mesher::Surface_mesher_manifold_base <Del, Oracle, Criteria, SMREWBB> > SMMWB;

typedef SM Surface_mesher;  // basic mesher
// typedef SMRE Surface_mesher;  // regular edges
// typedef SMM Surface_mesher;  // manifold with boundary
//typedef SMMWB Surface_mesher;  // manifold without boundary

// polyhedron
#include "lib/enriched_polyhedron.h"
#include <CGAL/IO/Polyhedron_iostream.h>
typedef Enriched_polyhedron<K2,Enriched_items> Mesh;

class CMeshDoc : public CDocument
{
protected: // create from serialization only
	CMeshDoc();
	DECLARE_DYNCREATE(CMeshDoc)

// Attributes
public:

	Del m_del;
  Mesh m_mesh;
  Oracle m_oracle;

	// options
	FT m_criterion_uniform_size;
	unsigned int m_max_nb_vertices;
	unsigned int m_refresh_each;
	unsigned int m_init_nb_vertices;

// Operations
public:

	void generate_initial_point_sample(std::list<Point>& points,
		                                 const unsigned int nb);
	void refresh(void);

	// status message
	void StatusMessage( char* fmt, ... );

// Overrides
	public:
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CMeshDoc();
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
	virtual BOOL OnSaveDocument(LPCTSTR lpszPathName);
	afx_msg void OnEditOptions();
	afx_msg void OnMeshingRun();
};
