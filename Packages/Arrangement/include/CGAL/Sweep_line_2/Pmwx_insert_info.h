// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Pmwx_insert_info.h
// package       : arr (1.87)
// maintainer    : Tali Zvi <talizvi@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PMWX_INSERT_INFO_H
#define CGAL_PMWX_INSERT_INFO_H

CGAL_BEGIN_NAMESPACE

template<class Vertex_handle, class Halfedge_handle>
class Pmwx_insert_info
{
public:

  /*! Constructor */
  Pmwx_insert_info() : 
    m_vertex(Vertex_handle(NULL)), 
    m_halfedge(Halfedge_handle(NULL))
  {
  }

  void setVertexHandle(Vertex_handle vh) {
    m_vertex = vh;
  }

  Vertex_handle getVertexHandle() const {
    return m_vertex;
  }

  void setHalfedgeHandle(Halfedge_handle h) {
    m_halfedge = h;
  }

  Halfedge_handle getHalfedgeHandle() const {
    return m_halfedge;
  }

  void Print()
  {
    if ( m_vertex == Vertex_handle(NULL))
      std::cout << "vertex: NULL\n";
    else
      std::cout << "vertex: " << m_vertex->point() << "\n";

    if ( m_halfedge == Halfedge_handle(NULL))
      std::cout << "halfedge: NULL\n";
    else 
      std::cout << "halfedge: " << m_halfedge->source()->point() 
		<< " " << m_halfedge->target()->point() << "\n";
  }
 
private:
  Vertex_handle m_vertex;
  Halfedge_handle m_halfedge;
  
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_INSERT_INFO_H

