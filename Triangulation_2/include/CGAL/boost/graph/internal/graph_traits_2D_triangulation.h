// Copyright (c) 2007, 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Andreas Fabri,
//                 Fernando Cacciola

#ifndef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
  #error CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION
  #error CGAL_2D_TRIANGULATION is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION_TEMPLATES
  #error CGAL_2D_TRIANGULATION_TEMPLATES is not defined
#endif

#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation_helper.h>

#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>
#include <CGAL/use.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>

#include <utility>

namespace boost {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_traits< CGAL_2D_TRIANGULATION >
{
  typedef CGAL_2D_TRIANGULATION                                               Triangulation;

  struct T2_graph_traversal_category :
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag
  { };

  typedef typename Triangulation::Vertex_handle                               vertex_descriptor;
  typedef CGAL::internal::T2_halfedge_descriptor<Triangulation>               halfedge_descriptor;
  typedef CGAL::internal::T2_edge_descriptor<Triangulation>                   edge_descriptor;
  typedef typename Triangulation::Face_handle                                 face_descriptor;

    // We need to go from 'Finite_vertex_iterator' to 'Vertex_handle' (and even more
    // in the case of RT2, since it has also a hidden filter)
  typedef CGAL::internal::Dereference_to_handle_enforcer<
            Triangulation,
            typename Triangulation::Finite_vertices_iterator,
            vertex_descriptor>                                                vertex_iterator;
  typedef CGAL::internal::Dereference_to_handle_enforcer<
            Triangulation,
            typename Triangulation::Finite_faces_iterator,
            face_descriptor>                                                  face_iterator;
  typedef CGAL::internal::T2_halfedge_iterator<Triangulation,
            typename Triangulation::Finite_edges_iterator>                    halfedge_iterator;
  typedef CGAL::internal::T2_edge_iterator<Triangulation,
            typename Triangulation::Finite_edges_iterator>                    edge_iterator;

  typedef CGAL::internal::T2_vertex_circulator<Triangulation>                 Vertex_circulator;
  typedef CGAL::Counting_iterator<Vertex_circulator>                          adjacency_iterator;
  typedef CGAL::internal::In_edge_circulator<Triangulation>                   In_edge_circ;
  typedef CGAL::Counting_iterator<In_edge_circ>                               in_edge_iterator;
  typedef CGAL::internal::Out_edge_circulator<Triangulation>                  Out_edge_circ;
  typedef CGAL::Counting_iterator<Out_edge_circ>                              out_edge_iterator;

  typedef undirected_tag                                                      directed_category;
  typedef disallow_parallel_edge_tag                                          edge_parallel_category;
  typedef T2_graph_traversal_category                                         traversal_category;

  typedef typename Triangulation::size_type                                   size_type;
  typedef size_type                                                           vertices_size_type;
  typedef size_type                                                           halfedges_size_type;
  typedef size_type                                                           edges_size_type;
  typedef size_type                                                           faces_size_type;
  typedef size_type                                                           degree_size_type;

  // nulls
  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }
  static face_descriptor null_face() { return face_descriptor(); }
};

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_traits<const CGAL_2D_TRIANGULATION >
  : public graph_traits< CGAL_2D_TRIANGULATION >
{ };

} // namespace boost

namespace CGAL {

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  CGAL_precondition(!g.is_infinite(e));
  return e.first->vertex(g.ccw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
source(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
       const CGAL_2D_TRIANGULATION& g)
{
  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));
  return h.first->vertex(g.ccw(h.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
target(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
       const CGAL_2D_TRIANGULATION& g)
{
  CGAL_precondition(!g.is_infinite(e));
  return e.first->vertex(g.cw(e.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor
target(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
       const CGAL_2D_TRIANGULATION& g)
{
  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));
  return h.first->vertex(g.cw(h.second));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
next(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor     face_descriptor;

  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));

  face_descriptor f = h.first;
  int i = h.second;
  if(!g.is_infinite(f))
    return halfedge_descriptor(f, g.ccw(i));

  // Now, 'h' is border halfedge, move on to the adjacent infinite face
  face_descriptor neigh_f = f->neighbor(g.ccw(i));
  return halfedge_descriptor(neigh_f, g.ccw(neigh_f->index(f)));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
prev(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor     face_descriptor;

  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));

  face_descriptor f = h.first;
  int i = h.second;
  if(!g.is_infinite(f))
    return halfedge_descriptor(f, g.cw(i));

  // Now, 'h' is border halfedge, move on to the adjacent infinite face
  face_descriptor neigh_f = f->neighbor(g.cw(i));
  return halfedge_descriptor(neigh_f, g.cw(neigh_f->index(f)));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
opposite(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));
  return halfedge_descriptor(g.mirror_edge(typename CGAL_2D_TRIANGULATION::Edge(h.first, h.second)));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor
face(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  if(g.is_infinite(h.first))
    return boost::graph_traits< CGAL_2D_TRIANGULATION >::null_face();
  else
    return h.first;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor f,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;

  CGAL_USE(g);
  CGAL_precondition(!g.is_infinite(f));
  return halfedge_descriptor(f, 0);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor face_descriptor;

  CGAL_precondition(!g.is_infinite(v));

  face_descriptor f = v->face();
  int i = f->index(v);

  while(g.is_infinite(f))
  {
    f = f->neighbor(g.cw(i));
    i = f->index(v);
  }

  return halfedge_descriptor(f, g.ccw(i));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor e,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;

  CGAL_USE(g);
  CGAL_precondition(!g.is_infinite(e));
  return halfedge_descriptor(e.first, e.second);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor
edge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor h,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;

  CGAL_USE(g);
  CGAL_precondition(!g.is_infinite(std::make_pair(h.first, h.second)));
  return edge_descriptor(h.first, h.second);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
std::pair<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor,
          bool>
edge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
     typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
     const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;

  CGAL_precondition(!g.is_infinite(u));
  CGAL_precondition(!g.is_infinite(v));

  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(u), done(c);
  if(c != 0)
  {
    do
    {
      if(!g.is_infinite(*c))
      {
        // find the index of the other vertex of *c
        int indv = 3 - c->first->index(u) - c->second;
        if(c->first->vertex(indv) == v)
          return std::make_pair(edge_descriptor(c->first, c->second), true);
      }
    } while (++c != done);
  }

  return std::make_pair(edge_descriptor(), false);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
std::pair<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor,
          bool>
halfedge(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor u,
         typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor face_descriptor;

  CGAL_precondition(!g.is_infinite(u));
  CGAL_precondition(!g.is_infinite(v));

  std::pair<edge_descriptor, bool> eb = edge(u, v, g);

  if(!eb.second)
    return std::make_pair(halfedge_descriptor(), false);

  const edge_descriptor& e = eb.first;

  if(e.first->vertex(g.ccw(e.first->index(u))) == v)
  {
    return std::make_pair(halfedge_descriptor(e.first, e.second), true);
  }
  else
  {
    face_descriptor nf = e.first->neighbor(e.second);
    int idx = nf->index(e.first);
    return std::make_pair(halfedge_descriptor(nf, idx), true);
  }
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_iterator>
vertices(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_iterator Iter;
  return make_range( Iter(g.finite_vertices_begin()), Iter(g.finite_vertices_end()) );
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_iterator >
halfedges(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_iterator Iter;
  return make_range(Iter(g.finite_edges_begin()), Iter(g.finite_edges_end()));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_iterator>
edges(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_iterator Iter;
  return make_range(Iter(g.finite_edges_begin()), Iter(g.finite_edges_end()));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator >
faces(const CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_iterator Iter;
  return make_range(Iter(g.finite_faces_begin()), Iter(g.finite_faces_end()));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
       const CGAL_2D_TRIANGULATION& g)
{
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = 0;
  typename CGAL_2D_TRIANGULATION::Edge_circulator c = g.incident_edges(v), done(c);
  if(c != 0)
  {
    do {
      if(!g.is_infinite(*c))
        ++deg;
    } while (++c != done);
  }

  return deg;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
in_degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
          const CGAL_2D_TRIANGULATION& g)
{
  return degree(v, g);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type
out_degree(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
           const CGAL_2D_TRIANGULATION& g)
{
  return degree(v, g);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::in_edge_iterator >
in_edges(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
         const CGAL_2D_TRIANGULATION& g)
{
  typedef CGAL::internal::In_edge_circulator< CGAL_2D_TRIANGULATION >                   Circ;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::in_edge_iterator     Iter;

  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(v, v->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = degree(v, g);

  return make_range(Iter(Circ(ec, g)), Iter(Circ(ec, g), deg));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::out_edge_iterator >
out_edges(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
          const CGAL_2D_TRIANGULATION& g)
{
  typedef CGAL::internal::Out_edge_circulator< CGAL_2D_TRIANGULATION >                  Circ;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::out_edge_iterator    Iter;

  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(v, v->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = degree(v, g);

  return make_range(Iter(Circ(ec, g)), Iter(Circ(ec, g), deg));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline Iterator_range<typename boost::graph_traits< CGAL_2D_TRIANGULATION >::adjacency_iterator>
adjacent_vertices(typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor v,
                  const CGAL_2D_TRIANGULATION& g)
{
  typedef CGAL::internal::T2_vertex_circulator< CGAL_2D_TRIANGULATION >                  Circ;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::adjacency_iterator    Iter;

  typename CGAL_2D_TRIANGULATION::Edge_circulator ec(v, v->face());
  typename boost::graph_traits< CGAL_2D_TRIANGULATION >::degree_size_type deg = degree(v, g);

  return make_range(Iter(Circ(ec, g)), Iter(Circ(ec, g), deg));
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertices_size_type
num_vertices(const CGAL_2D_TRIANGULATION& g)
{
  return g.number_of_vertices();
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::faces_size_type
num_faces(const CGAL_2D_TRIANGULATION& g)
{
  return g.number_of_faces();
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edges_size_type
num_edges(const CGAL_2D_TRIANGULATION& g)
{
  // Euler characteristic for a triangulated topological disk is V-E+F=1
  return num_vertices(g) + num_faces(g) - 1;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedges_size_type
num_halfedges(const CGAL_2D_TRIANGULATION& g)
{
  return num_edges(g) * 2;
}

} // namespace CGAL

#undef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
#undef CGAL_2D_TRIANGULATION
#undef CGAL_2D_TRIANGULATION_TEMPLATES
