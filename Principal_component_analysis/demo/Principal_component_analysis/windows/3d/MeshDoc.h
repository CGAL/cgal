// MeshDoc.h : interface of the CMeshDoc class
#pragma once

// CGAL
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include "lib/enriched_polyhedron.h"

#include <list>

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Enriched_polyhedron<Kernel,Enriched_items> Mesh;

class CMeshDoc : public CDocument
{
protected: // create from serialization only
	CMeshDoc();
	DECLARE_DYNCREATE(CMeshDoc)

// Attributes
public:

  // data sets
	std::list<Point> m_points;
	std::list<Triangle> m_triangles;

  // triangle mesh represented as a CGAL polyhedron
  Mesh m_mesh;

  // fitting primitives
  FT m_quality;
  Point m_centroid;
  Line m_fitting_line;
  Plane m_fitting_plane;

  void gl_draw_fitting_primitives();

// Operations
public:

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
  afx_msg void OnFitFitpointset();
  afx_msg void OnFitTriangleset32899();
};


