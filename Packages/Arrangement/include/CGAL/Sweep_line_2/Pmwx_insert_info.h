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

template<class _Halfedge_handle>
class Pmwx_insert_info
{
public:

  typedef _Halfedge_handle Halfedge_handle;

  /*! Constructor */
  Pmwx_insert_info() : 
    m_halfedge(Halfedge_handle(NULL))
  {
  }

  void setHalfedgeHandle(Halfedge_handle h) {
    m_halfedge = h;
  }

  Halfedge_handle getHalfedgeHandle() const {
    return m_halfedge;
  }

  void Print()
  {
    if ( m_halfedge == Halfedge_handle(NULL))
      std::cout << "halfedge: NULL\n";
    else 
      std::cout << "halfedge: " << m_halfedge->source()->point() 
		<< " " << m_halfedge->target()->point() << "\n";
  }
 
private:
  Halfedge_handle m_halfedge;
  
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_INSERT_INFO_H

