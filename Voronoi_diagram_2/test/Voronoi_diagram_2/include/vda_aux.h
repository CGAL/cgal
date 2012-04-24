// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_AUX_H
#define VDA_AUX_H 1

#include <CGAL/basic.h>
#include <cassert>
#include <iostream>
#include <CGAL/Triangulation_utils_2.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

template<class Vertex_handle, class Site_t>
struct Project_site
{
  typedef Site_t Site_2;
  Site_2 operator()(const Vertex_handle& v) const
  {
    return v->site();
  }
};


template<class Vertex_handle, class Site_t>
struct Project_point
{
  typedef Site_t Site_2;
  Site_2 operator()(const Vertex_handle& v) const
  {
    return v->point();
  }
};


template<class VDA, class Point_t>
struct Project_dual
{
  typedef Point_t Point_2;
  typedef typename VDA::Delaunay_graph::Face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    return vda.dual().dual( f );
  }
};

template<class VDA, class Point_t>
struct Project_primal
{
  typedef Point_t Point_2;
  typedef typename VDA::Delaunay_graph::Face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    return vda.dual().primal( f );
  }
};

template<class VDA, class Site_t>
struct Project_ag_dual
{
  typedef Site_t Site_2;
  typedef typename Site_2::Point_2 Point_2;
  typedef typename VDA::Delaunay_graph::Face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    CGAL::Object o = vda.dual().dual(f);
    Site_2 s;
    if ( CGAL::assign(s, o) ) {
      return s.point();
    } else{
      assert( false );
      return Point_2();
    }
  }
};



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


template<class VDA, class Projector>
void print_halfedge(const VDA& vda,
		    const typename VDA::Halfedge_handle& he,
		    const Projector& project,
		    std::ostream& os = std::cout) {
  typename VDA::Delaunay_graph::Edge e = he->dual_edge();
  print_dual_edge(vda, e, project, os);
}

template<class VDA, class Projector>
void print_halfedge(const VDA& vda,
		    const typename VDA::Halfedge& he,
		    const Projector& project,
		    std::ostream& os = std::cout) {
  typename VDA::Delaunay_graph::Edge e = he.dual_edge();
  print_dual_edge(vda, e, project, os);
}


template<class VDA, class Projector>
void print_dual_edge(const VDA& vda,
		     const typename VDA::Delaunay_graph::Edge& e,
		     const Projector& project,
		     std::ostream& os = std::cout) {
  typedef CGAL::Triangulation_cw_ccw_2  CW_CCW_2;

  print_dual_edge(vda,
		  e.first->vertex( CW_CCW_2::ccw(e.second) ),
		  e.first->vertex( CW_CCW_2::cw(e.second) ),
		  project,
		  os);
}


template<class VDA, class Projector>
void print_dual_edge(const VDA& vda,
		     const typename VDA::Delaunay_graph::Vertex_handle& v_src,
		     const typename VDA::Delaunay_graph::Vertex_handle& v_trg,
		     const Projector& project,
		     std::ostream& os = std::cout)
{
  if ( vda.dual().is_infinite(v_src) ) {
    os << "inf";
  } else {
    os << project(v_src);
  }
  os << " - ";
  if ( vda.dual().is_infinite(v_trg) ) {
    os << "inf";
  } else {
    os << project(v_trg);
  }
  os << std::endl;
}




template<class VDA, class Projector>
void print_dual_vertex(const VDA& vda,
		       const typename VDA::Delaunay_graph::Vertex_handle& v,
		       const Projector& project,
		       std::ostream& os = std::cout) {
   if ( vda.dual().is_infinite(v) ) {
     os << "inf" << std::endl;
   } else {
     os << project(v) << std::endl;
   }
}


void print_separator(std::ostream& os) {
  os << std::endl << std::endl;
  os << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << std::endl;
  os << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
     << std::endl;
  os << std::endl << std::endl;
}


#endif // VDA_AUX_H
