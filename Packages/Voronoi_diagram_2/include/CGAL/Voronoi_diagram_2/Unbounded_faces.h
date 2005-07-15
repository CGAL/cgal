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

#ifndef CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
#define CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================

template<class VDA, class Base_it>
class Bounded_face_tester
{
 private:
  // this class returns true if the face is bounded

  // this is essentially VDA::Non_degenerate_faces_iterator
  typedef Base_it                                      Base_iterator;
  typedef typename VDA::Delaunay_graph                 Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_circulator   Dual_vertex_circulator;
  typedef typename Delaunay_graph::Vertex_handle       Dual_vertex_handle;

 public:
  Bounded_face_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return false; }

    Dual_vertex_handle v = it.base();

    Dual_vertex_circulator vc = vda_->dual().incident_vertices(v);
    Dual_vertex_circulator vc_start = vc;
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
  typedef typename Delaunay_graph::Vertex_circulator   Dual_vertex_circulator;
  typedef typename Delaunay_graph::Vertex_handle       Dual_vertex_handle;

 public:
  Unbounded_face_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Base_iterator& it) const {
    if ( vda_->dual().dimension() < 2 ) { return true; }

    Dual_vertex_handle v = it.base();

    Dual_vertex_circulator vc = vda_->dual().incident_vertices(v);
    Dual_vertex_circulator vc_start = vc;
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

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_UNBOUNDED_FACES_H
