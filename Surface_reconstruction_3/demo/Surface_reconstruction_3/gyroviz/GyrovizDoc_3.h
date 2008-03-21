// GyrovizDoc_3.h : interface of the CGyrovizDoc_3 class
//
#ifndef _DOC_
#define _DOC_
#pragma once


// CGAL
#include "GyrovizKernel.h"
#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

// This demo
#include "Gyroviz_dt2.h"
#include "Gyroviz_dt3.h"


// Gyroviz's Delaunay triangulation 2-3
typedef CGAL::Triangulation_vertex_base_with_info_2<Gyroviz_info_for_dt2,K> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef Gyroviz_dt2<K,Tds>         Dt2;

typedef Dt2::Face_handle              Face_handle;
typedef Dt2::Finite_faces_iterator    Finite_faces_iterator;
typedef Dt2::Finite_edges_iterator    Finite_edges_iterator;
typedef Dt2::Finite_vertices_iterator Finite_vertices_iterator;



// Gyroviz's Delaunay triangulation 3
typedef CGAL::Triangulation_vertex_base_with_info_3<Gyroviz_info_for_dt3,K> Vb3;
typedef CGAL::Triangulation_cell_base_with_info_3<int,K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb3,Cb> Tds3;
typedef Gyroviz_dt3<K,Tds3>         Dt3;

//typedef Dt3::Finite_facets_iterator   Finite_facets_iterator;
//typedef Dt3::Finite_edges_iterator    Finite_edges_iterator;
//typedef Dt3::Finite_vertices_iterator Finite_vertices_iterator;
//typedef Dt3::Finite_cells_iterator    Finite_cells_iterator;


class CGyrovizDoc_3 : public CDocument
{
protected: // create from serialization only
	CGyrovizDoc_3();
	DECLARE_DYNCREATE(CGyrovizDoc_3)

// Data members
private:
  // Triangulation
  Dt2 m_gyroviz_dt; // The Gyroviz equation is solved on the vertices of m_gyroviz_dt
  Dt3 m_gyroviz_dt3;

	// options
	double m_point_size; // OpenGL point size
	
// Public methods
public:

  // Get triangulation.
  Dt2& get_dt2()
  {
    return m_gyroviz_dt;
  }
  const Dt2& get_dt2() const
  {
    return m_gyroviz_dt;
  }

  // Get triangulation.
  Dt3& get_dt3()
  {
    return m_gyroviz_dt3;
  }
  const Dt3& get_dt3() const
  {
    return m_gyroviz_dt3;
  }
  
  // cuurent objet is a dt2 
  bool is_dt2(){
    return (m_gyroviz_dt.number_of_vertices() > 0);
  }


  // Get OpenGL point size.
  double get_point_size() const
  {
    return m_point_size;
  }

// Private methods
private:

	// misc status stuff
	void update_status();
	void status_message(char* fmt,...);
	double duration(const double time_init);

// MFC generated
public:
	virtual ~CGyrovizDoc_3();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
	afx_msg void OnEditOptions();
  afx_msg void OnFileSaveSurface();
  afx_msg void OnUpdateFileSaveSurface(CCmdUI *pCmdUI);
  afx_msg void OnFileSaveAs();
  afx_msg void OnUpdateFileSaveAs(CCmdUI *pCmdUI);
};


#endif // _DOC_