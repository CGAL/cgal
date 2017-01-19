// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_REFINE_EDGES_VISITOR_H
#define CGAL_MESH_2_REFINE_EDGES_VISITOR_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Mesher_level.h>

namespace CGAL {
namespace Mesh_2 {

/**
 * This class is the visitor needed when Refine_edges<Tr> if called from 
 * Refine_faces<Tr>.
 * \param Faces_mesher should be instanciated with Refine_face_base<Tr>.
 */
template <typename Faces_mesher>
class Refine_edges_visitor : public ::CGAL::Null_mesh_visitor
{
public:
  typedef typename Faces_mesher::Triangulation Tr;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Point Point;
  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

  typedef typename details::Refine_edges_base_types<Tr>::Constrained_edge
                               Constrained_edge;

  typedef typename Faces_mesher::Previous_level Edges_mesher;
private:
  Faces_mesher& faces_mesher;
  Edges_mesher& edges_mesher;
  Vertex_handle &va, &vb;
  bool &mark_at_left, &mark_at_right;
  Null_mesh_visitor& null_mesh_visitor;

public:
  Refine_edges_visitor(Faces_mesher& faces_mesher_,
		       Edges_mesher& edges_mesher_,
                       Null_mesh_visitor& null)
    : faces_mesher(faces_mesher_),
      edges_mesher(edges_mesher_),
      va(edges_mesher.visitor_va),
      vb(edges_mesher.visitor_vb),
      mark_at_left(edges_mesher.visitor_mark_at_left),
      mark_at_right(edges_mesher.visitor_mark_at_right),
      null_mesh_visitor(null)
  {
  }

  Null_mesh_visitor previous_level() const { return null_mesh_visitor; }
  
  /** 
   * Store vertex handles and markers at left and right of the edge \c e.
   */
  void before_conflicts(const Edge& e, const Point&)
  {
    const Face_handle& fh = e.first;
    const int edge_index = e.second;

    va = fh->vertex(Tr::cw (edge_index));
    vb = fh->vertex(Tr::ccw(edge_index));
    
    mark_at_right = fh->is_in_domain();
    mark_at_left = fh->neighbor(edge_index)->is_in_domain();
  }

  void before_insertion(const Edge&, const Point& p, Zone& z)
  {
    faces_mesher.before_insertion_impl(Face_handle(), p, z);
  }

  /** Restore markers in the star of \c v. */
  void after_insertion(const Vertex_handle& v)
  {
    Tr& tr = faces_mesher.triangulation_ref_impl();

    int dummy;
    // if we put edge_index instead of dummy, Intel C++ does not find
    // a matching function for is_edge
    Face_handle fh;

    tr.is_edge(va, v, fh, dummy);
    // set fh to the face at the right of [va,v]

    typename Tr::Face_circulator fc = tr.incident_faces(v, fh), fcbegin(fc);
    // circulators are counter-clockwise, so we start at the right of
    // [va,v]
    do {
      if( !tr.is_infinite(fc) )
        fc->set_in_domain(mark_at_right);
      ++fc;
    } while ( fc->vertex(tr.ccw(fc->index(v))) != vb );
    // we are now at the left of [va,vb]
    do {
      if( !tr.is_infinite(fc) )
        fc->set_in_domain(mark_at_left);
      ++fc;
    } while ( fc != fcbegin );

    // then let's update bad faces
    faces_mesher.compute_new_bad_faces(v);

    CGAL_expensive_assertion(faces_mesher.check_bad_faces());
  }

  template <typename E, typename P, typename Z>
  void after_no_insertion(E, P, Z) const {}
}; // end class Refine_edges_visitor

/**
 * This class is the visitor for Refine_faces<Tr>.
 */
template <typename Faces_mesher>
class Refine_edges_visitor_from_faces
{
public:
  typedef Refine_edges_visitor<Faces_mesher> Previous_level;
  typedef typename Faces_mesher::Previous_level Edges_mesher;

private:
  Previous_level previous;

public:

  Refine_edges_visitor_from_faces(Faces_mesher& faces_mesher,
				  Edges_mesher& edges_mesher,
                                  Null_mesh_visitor& null)
    : previous(faces_mesher, edges_mesher, null)
  {
  }

  Previous_level& previous_level() { return previous; }

  template <typename E, typename P>
  void before_conflicts(E, P) const {}

  template <typename E, typename P, typename Z>
  void before_insertion(E, P, Z) const {}

  template <typename V>
  void after_insertion(V) const {}

  template <typename E, typename P, typename Z>
  void after_no_insertion(E, P, Z) const {}
}; // end class Refine_edges_visitor_from_faces

} // end namespace Mesh_2
} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_EDGES_VISITOR_H
