// ======================================================================
//
// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_POLYHEDRON_DECORATOR_H_01282002
#define CGAL_POLYHEDRON_DECORATOR_H_01282002

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>

namespace CGAL {

template <class Poly>
class Polyhedron_decorator_3 {
  typedef Poly                                                      Polyhedron;

  typedef typename boost::property_map<Polyhedron, vertex_point_t>::type Vertex_pmap;
  typedef typename boost::property_traits<Vertex_pmap>::value_type  Point;

  typedef typename Kernel_traits<Point>::Kernel                     Kernel;
  typedef typename Kernel::FT                                       FT;

  typedef typename boost::graph_traits<Poly>::vertex_descriptor     Vertex_handle;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor   Halfedge_handle;
  typedef typename boost::graph_traits<Poly>::face_descriptor       Facet_handle;

  typedef typename boost::graph_traits<Poly>::vertex_iterator       Vertex_iterator;
  typedef typename boost::graph_traits<Poly>::edge_iterator         Edge_iterator;
  typedef typename boost::graph_traits<Poly>::face_iterator         Facet_iterator;

  typedef Halfedge_around_face_circulator<Poly>     Halfedge_around_facet_circulator;
  typedef Halfedge_around_target_circulator<Poly>   Halfedge_around_vertex_circulator;

public:
  /** Insert a new vertex into a helfedge h (a--b)

      Precondition:
           h
      a <-----> b
          -h
      h is the halfedge connecting vertex a to b
      -h is the opposite halfedge connecting b to a

      Postcondition:
           h         r
      a <-----> V <-----> b
          -h         -r
      V is the return vertex whose geometry is UNDEFINED.
      -r is the returned halfedge that is pointing to V
  */

#if 0
  static Vertex* insert_vertex(Polyhedron& p, Halfedge* h) {
    return insert_vertex(p, Halfedge_handle(h)).ptr();
  }
  //
  static Vertex* insert_vertex(Polyhedron& p, Vertex* a, Vertex* b) {
    return insert_vertex(p, Vertex_handle(a), Vertex_handle(b)).ptr();
  }
#endif
  //
  static Vertex_handle insert_vertex(Polyhedron& p, Halfedge_handle h) {
    return target(insert_vertex_return_edge(p, h),p);
  }
  //
  static Vertex_handle insert_vertex(Polyhedron& p,
                                     Vertex_handle a, Vertex_handle b) {
    return target(insert_vertex_return_edge(p, a, b),p);
  }

#if 0
  //
  static Halfedge* insert_vertex_return_edge(Polyhedron& p, Halfedge* h) {
    return insert_vertex_return_edge(p, Halfedge_handle(h)).ptr();
  }
  //
  static Halfedge* insert_vertex_return_edge(Polyhedron& p,
                                             Vertex* a, Vertex* b) {
    return insert_vertex_return_edge(p, Vertex_handle(a), Vertex_handle(b)).ptr();
  }

#endif
  //
  static inline Halfedge_handle insert_vertex_return_edge(Polyhedron& p,
                                                          Halfedge_handle h);
  //
  static inline Halfedge_handle insert_vertex_return_edge(Polyhedron& p,
                                                          Vertex_handle a,
                                                          Vertex_handle b);

  /** Insert a new edge (two halfedges) between the two vertices

      Precondition:
      vertex a and b are in the SAME facet and do NOT connect to each other

      Postcondition:
             H
      a <----------> b
      H is the return halfedge connecting vertex a to b.
  */
#if 0
  static Halfedge* insert_edge(Polyhedron& /*p*/, Vertex* a, Vertex* b) {
    return insert_edge(Vertex_handle(a), Vertex_handle(b)).ptr();
  }
#endif
  //
  static Halfedge_handle insert_edge(Polyhedron& p,
                                     Halfedge_handle a,
                                     Halfedge_handle b) {
    return Euler::split_face(a, b, p);
  }
  //
  static inline Halfedge_handle insert_edge(Polyhedron& p,
                                            Vertex_handle a,
                                            Vertex_handle b);

  /** Set the vertex of index vidx to the new position to the new
      position (x,y,z)

      Precondition:
      0 <= vidx < p.size_of_vertices()

      Postcondition:
      If failed procondition, do nothing.
      The new vertex of index vidx has the new position (x,y,z)
  */
  static void set_vertex_position(Polyhedron& p, int vidx, FT x, FT y, FT z) {
    if (vidx >= 0 && vidx < p.size_of_vertices()) {
      Vertex_iterator vitr = p.vertices_begin();
      for (int i = 0; i < vidx; i++) ++vitr;
      vitr->point() = Point(x,y,z);
    }
  }
};

// ======================================================================
//
template <class Poly>
typename Polyhedron_decorator_3<Poly>::Halfedge_handle
Polyhedron_decorator_3<Poly>::insert_vertex_return_edge(Polyhedron& p,
                                                        Halfedge_handle h) {
  Halfedge_handle hopp = opposite(h,p);
  Halfedge_handle r = Euler::split_vertex(prev(hopp,p), h,p);
  if (! is_border(h,p))
    set_halfedge(face(h,p), r, p);
  if (! is_border(hopp,p))
    set_halfedge(face(hopp, p), hopp, p);
  return opposite(r,p);
}

// ======================================================================
//
template <class Poly>
typename Polyhedron_decorator_3<Poly>::Halfedge_handle
Polyhedron_decorator_3<Poly>::insert_vertex_return_edge(Polyhedron& p,
                                                        Vertex_handle a,
                                                        Vertex_handle b) {
  Halfedge_around_vertex_circulator a_cir_begin = a->vertex_begin();
  Halfedge_around_vertex_circulator a_cir = a_cir_begin;
  do {
    if (a_cir->opposite()->vertex() == b)
      return insert_vertex_return_edge(p, a_cir->opposite());
  } while (++a_cir != a_cir_begin);
  CGAL_precondition_msg(0, "vertex a and b must be incident to each other");
  return Halfedge_handle(NULL);
}

// ======================================================================
//
template <class Poly>
typename Polyhedron_decorator_3<Poly>::Halfedge_handle
Polyhedron_decorator_3<Poly>::insert_edge(Polyhedron& p,
                                          Vertex_handle a,
                                          Vertex_handle b) {
  Halfedge_around_vertex_circulator a_cir_begin = a->vertex_begin();
  Halfedge_around_vertex_circulator a_cir = a_cir_begin;
  Halfedge_around_vertex_circulator b_cir_begin = b->vertex_begin();
  Halfedge_around_vertex_circulator b_cir = b_cir_begin;
  do {
    do {
      if (a_cir->facet() == b_cir->facet())
        return p.split_facet(a_cir, b_cir);
    } while (++b_cir != b_cir_begin);
  } while (++a_cir != a_cir_begin);
  CGAL_precondition_msg(0, "vertex a and b must share the same face");
  return Halfedge_handle(NULL);
}

} //namespace CGAL

#endif //CGAL_POLYHEDRON_DECORATOR_H_01282002
