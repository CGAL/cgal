// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H
#define CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
//#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================

template<class VDA, class Base_it>
class Bounded_edge_tester
{
 private:
  // this class returns true if the edge is bounded

  // this is essentially VDA::Non_degenerate_edges_iterator
  typedef Base_it                                     Base_iterator;
  typedef typename VDA::Delaunay_graph::Edge          Delaunay_edge;
  typedef typename VDA::Delaunay_graph::Face_handle   Delaunay_face_handle;

 public:
  Bounded_edge_tester(const VDA* vda = nullptr) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return false; }

    Delaunay_edge e = it->dual();

    Delaunay_face_handle df1 = e.first;
    Delaunay_face_handle df2 = e.first->neighbor(e.second);

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
  typedef Base_it                                     Base_iterator;
  typedef typename VDA::Delaunay_graph::Edge          Delaunay_edge;
  typedef typename VDA::Delaunay_graph::Face_handle   Delaunay_face_handle;

 public:
  Unbounded_edge_tester(const VDA* vda = nullptr) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return true; }

    Delaunay_edge e = it->dual();

    Delaunay_face_handle df1 = e.first;
    Delaunay_face_handle df2 = e.first->neighbor(e.second);

    return vda_->dual().is_infinite(df1) || vda_->dual().is_infinite(df2);
  }

 private:
  const VDA* vda_;
};

//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_EDGES_H
