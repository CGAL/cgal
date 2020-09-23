// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H

#include <functional>

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation_helper.h>
#include <CGAL/boost/graph/properties_Triangulation_data_structure_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/iterator.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

// The functions and classes in this file allows the user to
// treat a CGAL Triangulation_data_structure_2 object as a boost graph "as is". No
// wrapper is needed for the Triangulation_data_structure_2 object.

namespace boost {

template <class VB, class FB>
struct graph_traits<CGAL::Triangulation_data_structure_2<VB, FB> >
{
  struct TDS2_graph_traversal_category :
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

  typedef CGAL::Triangulation_data_structure_2<VB,FB> Triangulation_data_structure;

  typedef typename Triangulation_data_structure::Vertex_handle                  vertex_descriptor;
  typedef CGAL::internal::T2_halfedge_descriptor<Triangulation_data_structure>  halfedge_descriptor;
  typedef CGAL::internal::T2_edge_descriptor<Triangulation_data_structure>      edge_descriptor;
  typedef typename Triangulation_data_structure::Face_handle                    face_descriptor;

  typedef CGAL::Prevent_deref<typename Triangulation_data_structure::Vertex_iterator> vertex_iterator;
  typedef CGAL::internal::T2_halfedge_iterator<Triangulation_data_structure,
            typename Triangulation_data_structure::Edge_iterator>                     halfedge_iterator;
  typedef CGAL::internal::T2_edge_iterator<Triangulation_data_structure,
            typename Triangulation_data_structure::Edge_iterator>                     edge_iterator;
  typedef CGAL::Prevent_deref<typename Triangulation_data_structure::Face_iterator>   face_iterator;

  typedef CGAL::Counting_iterator<CGAL::internal::TDS2_Out_edge_circulator<typename Triangulation_data_structure::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
  typedef CGAL::Counting_iterator<CGAL::internal::TDS2_In_edge_circulator<typename Triangulation_data_structure::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
  typedef CGAL::Counting_iterator<typename Triangulation_data_structure::Vertex_circulator> Incident_vertices_iterator;
  typedef Incident_vertices_iterator adjacency_iterator;

  typedef undirected_tag directed_category;
  typedef disallow_parallel_edge_tag edge_parallel_category;
  typedef TDS2_graph_traversal_category traversal_category;

  typedef typename Triangulation_data_structure::size_type                          size_type;
  typedef size_type                                                                 vertices_size_type;
  typedef size_type                                                                 edges_size_type;
  typedef size_type                                                                 halfedges_size_type;
  typedef size_type                                                                 faces_size_type;
  typedef size_type                                                                 degree_size_type;

  // nulls
  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
  static face_descriptor null_face()   { return face_descriptor(); }
  static halfedge_descriptor null_halfedge()   { return halfedge_descriptor(); }
};

template <class VB, class FB>
struct graph_traits<const CGAL::Triangulation_data_structure_2<VB, FB> >
  : public graph_traits< CGAL::Triangulation_data_structure_2<VB, FB> >
{ };

} // namespace boost

namespace CGAL {

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
next(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
prev(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
opposite(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
         const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename Triangulation_data_structure_2<VB,FB>::Edge Edge;
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(g.mirror_edge(Edge(e.first, e.second)));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
source(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
target(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
source(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::ccw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor
target(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
       const Triangulation_data_structure_2<VB,FB>& )
{
  return e.first->vertex(Triangulation_data_structure_2<VB,FB>::cw(e.second));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor
face(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>&)
{
  return e.first;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor f,
         const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(f,0);
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor v,
         const Triangulation_data_structure_2<VB,FB>& )
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_descriptor face_descriptor;
  face_descriptor fd = v->face();
  int i = fd->index(v);
  return halfedge_descriptor(fd,Triangulation_data_structure_2<VB,FB>::ccw(i));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor e,
         const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(e.first, e.second);
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor
edge(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_descriptor e,
     const Triangulation_data_structure_2<VB,FB>&)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_descriptor edge_descriptor;
  return edge_descriptor(e.first,e.second);
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_iterator>
vertices(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_iterator Iter;
  return make_range(Iter(g.vertices_begin()), Iter(g.vertices_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edge_iterator>
edges(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits<Triangulation_data_structure_2<VB,FB> >::edge_iterator Iter;
  return make_range(Iter(g.edges_begin()), Iter(g.edges_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_iterator >
halfedges(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedge_iterator Iter;
  return make_range(Iter(g.edges_begin()), Iter(g.edges_end()));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_iterator >
faces(const Triangulation_data_structure_2<VB,FB>& g)
{
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::face_iterator Iter;
  return make_range(Iter(g.faces_begin()), Iter(g.faces_end()));
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
out_degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
           const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::out_edge_iterator >
out_edges(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
          const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator ec(u,u->face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::out_edge_iterator Iter;

  return make_range(Iter(ec), Iter(ec,out_deg));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::in_edge_iterator >
in_edges(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
         const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator ec(u,u->face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::in_edge_iterator Iter;
  return make_range(Iter(ec), Iter(ec,out_deg));
}

template <class VB, class FB>
inline Iterator_range<typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::adjacency_iterator>
adjacent_vertices(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
                  const Triangulation_data_structure_2<VB,FB>& g)
{
  typename Triangulation_data_structure_2<VB,FB>::Vertex_circulator vc = out_edge_iterator(u,u.face());
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type out_deg = out_degree(u,g);
  typedef typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::adjacency_iterator Iter;
  return make_range( Iter(vc), Iter(vc,out_deg) );
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertices_size_type
num_vertices(const Triangulation_data_structure_2<VB,FB>& g)
{
  return g.number_of_vertices();
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::edges_size_type
num_edges(const Triangulation_data_structure_2<VB,FB>& g)
{
  return  g.number_of_vertices() + g.number_of_faces() - 2;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::halfedges_size_type
num_halfedges(const Triangulation_data_structure_2<VB,FB>& g)
{
  return num_edges(g) * 2;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::faces_size_type
num_faces(const Triangulation_data_structure_2<VB,FB>& g)
{
  return g.number_of_faces();
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
in_degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
          const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

template <class VB, class FB>
typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type
degree(typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::vertex_descriptor u,
       const Triangulation_data_structure_2<VB,FB>& g)
{
  typename boost::graph_traits< Triangulation_data_structure_2<VB,FB> >::degree_size_type deg = 0;
  typename Triangulation_data_structure_2<VB,FB>::Edge_circulator c = g.incident_edges(u), done(c);
  if ( c != 0) {
    do {
      ++deg;
    } while (++c != done);
  }
  return deg;
}

} // namespace CGAL

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_DATA_STRUCTURE_2_H
