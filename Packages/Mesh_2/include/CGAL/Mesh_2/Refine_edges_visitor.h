// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef MESH_2_REFINE_EDGES_VISITOR_H
#define MESH_2_REFINE_EDGES_VISITOR_H

#include <CGAL/Mesher_level.h>

namespace CGAL {
namespace Mesh_2 {

/**
 * This class is the visitor needed when Refine_edges<Tr> if called from 
 * Refine_faces<Tr>.
 *
 * \param Mesher_base should be instanciated with Refine_face_base<Tr>.
 */
template <typename Mesher_base>
class Refine_edges_visitor : public ::CGAL::Null_mesh_visitor
{
public:
  typedef typename Mesher_base::Triangulation Tr;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Point Point;
  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

  typedef typename details::Refine_edges_base_types<Tr>::Constrained_edge
                               Constrained_edge;
private:
  Mesher_base& mesher_base;
  Vertex_handle va, vb;
  bool mark_at_right, mark_at_left;
  Null_mesh_visitor& null_mesh_visitor;

public:
  Refine_edges_visitor(Mesher_base& mesher_base_, Null_mesh_visitor& null)
    : mesher_base(mesher_base_), null_mesh_visitor(null)
  {
  }

  Null_mesh_visitor previous_level() const { return null_mesh_visitor; }
  
  /** 
   * Store vertex handles and markers at left and right of the edge \c e.
   */
  void before_conflicts(const Constrained_edge& e, const Point&)
  {
    Tr& tr = mesher_base.get_triangulation_ref();

    va = e.first;
    vb = e.second;

    Face_handle fh;
    int edge_index;
    CGAL_assertion_code( bool should_be_true = )
    tr.is_edge(va, vb, fh, edge_index);
    CGAL_assertion( should_be_true );
    
    mark_at_right = fh->is_marked();
    mark_at_left = fh->neighbor(edge_index)->is_marked();
  }

  /** Restore markers in the star of \c v. */
  void after_insertion(const Vertex_handle& v)
  {
    Tr& tr = mesher_base.get_triangulation_ref();

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
        fc->set_marked(mark_at_right);
      ++fc;
    } while ( fc->vertex(tr.ccw(fc->index(v))) != vb );
    // we are now at the left of [va,vb]
    do {
      if( !tr.is_infinite(fc) )
        fc->set_marked(mark_at_left);
      ++fc;
    } while ( fc != fcbegin );

    // then let's update bad faces
    mesher_base.compute_new_bad_faces(v);
  }

  template <typename E, typename P, typename Z>
  void before_insertion(E, P, Z) const {}

  template <typename E, typename P, typename Z>
  void after_no_insertion(E, P, Z) const {}
}; // end class Refine_edges_visitor

/**
 * This class is the visitor for Refine_faces<Tr>.
 */
template <typename Mesher_base>
class Refine_edges_visitor_from_faces
{
public:
  typedef Refine_edges_visitor<Mesher_base> Previous_level;

private:
  Previous_level previous;

public:

  Refine_edges_visitor_from_faces(Mesher_base& mesher_base,
                                  Null_mesh_visitor& null)
    : previous(mesher_base, null)
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

#endif // MESH_2_REFINE_EDGES_VISITOR_H
