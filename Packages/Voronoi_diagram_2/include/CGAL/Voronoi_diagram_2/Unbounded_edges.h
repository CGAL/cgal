// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H
#define CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
//#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================

template<class VDA, class Base_it>
class Bounded_edge_tester
{
 private:
  // this class returns true if the edge is bounded

  // this is essentially VDA::Non_degenerate_edges_iterator
  typedef Base_it                                 Base_iterator;
  typedef typename VDA::Dual_graph::Edge          Dual_edge;
  typedef typename VDA::Dual_graph::Face_handle   Dual_face_handle;

  //  typedef Triangulation_cw_ccw_2    CW_CCW_2;

 public:
  Bounded_edge_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return false; }

    Dual_edge e = it->dual_edge();

    Dual_face_handle df1 = e.first;
    Dual_face_handle df2 = e.first->neighbor(e.second);

    return !vda_->dual().is_infinite(df1) && !vda_->dual().is_infinite(df2);
  }

 private:
  const VDA* vda_;
};

//=========================================================================

template<class VDA, class Base_it>
class Unbounded_edge_tester
{
 private:
  // this class returns true if the edge is unbounded

  // this is essentially VDA::Non_degenerate_edges_iterator
  typedef Base_it                                 Base_iterator;
  typedef typename VDA::Dual_graph::Edge          Dual_edge;
  typedef typename VDA::Dual_graph::Face_handle   Dual_face_handle;

 public:
  Unbounded_edge_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return true; }

    Dual_edge e = it->dual_edge();

    Dual_face_handle df1 = e.first;
    Dual_face_handle df2 = e.first->neighbor(e.second);

    return vda_->dual().is_infinite(df1) || vda_->dual().is_infinite(df2);
  }

 private:
  const VDA* vda_;
};

//=========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H
