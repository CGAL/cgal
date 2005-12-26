// MeshDoc.h : interface of the CMeshDoc class
#pragma once

// CGAL
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "lib/enriched_polyhedron.h"

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Enriched_polyhedron<Kernel,Enriched_items> Mesh;

class CMeshDoc : public CDocument
{
protected: // create from serialization only
	CMeshDoc();
	DECLARE_DYNCREATE(CMeshDoc)

// Attributes
public:

  // point set
  std::list<Point> m_points;

  // triangle mesh represented as a polyhedron
  Mesh m_mesh;

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
};


