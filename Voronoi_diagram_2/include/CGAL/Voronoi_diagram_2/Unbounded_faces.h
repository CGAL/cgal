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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
#define CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================

template<class VDA, class Base_it>
class Bounded_face_tester
{
 private:
  // this class returns true if the face is bounded

  // this is essentially VDA::Non_degenerate_faces_iterator
  typedef Base_it                                      Base_iterator;
  typedef typename VDA::Delaunay_graph                 Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle       Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Vertex_circulator
  Delaunay_vertex_circulator;

 public:
  Bounded_face_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return false; }

    Delaunay_vertex_handle v = it.base();

    Delaunay_vertex_circulator vc = vda_->dual().incident_vertices(v);
    Delaunay_vertex_circulator vc_start = vc;
    do {
      if ( vda_->dual().is_infinite(vc) ) { return false; }
      ++vc;
    } while ( vc != vc_start );
    return true;
  }

 private:
  const VDA* vda_;
};

//=========================================================================

template<class VDA, class Base_it>
class Unbounded_face_tester
{
 private:
  // this class returns true if the face is unbounded

  // this is essentially VDA::Non_degenerate_faces_iterator
  typedef Base_it                                      Base_iterator;
  typedef typename VDA::Delaunay_graph                 Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle       Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Vertex_circulator
  Delaunay_vertex_circulator;

 public:
  Unbounded_face_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return true; }

    Delaunay_vertex_handle v = it.base();

    Delaunay_vertex_circulator vc = vda_->dual().incident_vertices(v);
    Delaunay_vertex_circulator vc_start = vc;
    do {
      if ( vda_->dual().is_infinite(vc) ) { return true; }
      ++vc;
    } while ( vc != vc_start );
    return false;
  }

 private:
  const VDA* vda_;
};

//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
