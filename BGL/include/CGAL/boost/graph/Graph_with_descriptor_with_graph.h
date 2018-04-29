// Copyright (c) 2016  GeometryFactory (France).  All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_GWDWG_H
#define CGAL_BOOST_GRAPH_GWDWG_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <CGAL/boost/graph/Graph_with_descriptor_with_graph_fwd.h>

namespace CGAL
{


template<typename Graph_, typename Descriptor_>
struct Gwdwg_descriptor
{
public:
  typedef Graph_ Graph;
  typedef Descriptor_ Descriptor;

  Graph* graph;
  Descriptor descriptor;

  Gwdwg_descriptor()
    : graph(NULL), descriptor()
  {}

  Gwdwg_descriptor(Descriptor descriptor)
    : graph(NULL), descriptor(descriptor)
  {}

  Gwdwg_descriptor(Descriptor descriptor, Graph& graph)
    : graph(&graph), descriptor(descriptor)
  {}
};

template<typename Graph,typename Descriptor>
bool operator==(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  CGAL_assertion( lhs.graph == rhs.graph || rhs.graph==NULL || lhs.graph==NULL);
  return lhs.descriptor == rhs.descriptor;
}

template<typename Graph,typename Descriptor>
bool operator!=(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  return ! (lhs == rhs);
}

template<typename Graph,typename Descriptor>
bool operator<(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  CGAL_assertion( lhs.graph == rhs.graph || rhs.graph==NULL || lhs.graph==NULL);
  return lhs.descriptor < rhs.descriptor;
}

template<typename Graph,typename Descriptor>
bool operator>(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  CGAL_assertion( lhs.graph == rhs.graph || rhs.graph==NULL || lhs.graph==NULL);
  return lhs.descriptor > rhs.descriptor;
}

template<typename Graph,typename Descriptor>
bool operator<=(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  CGAL_assertion( lhs.graph == rhs.graph || rhs.graph==NULL || lhs.graph==NULL);
  return lhs.descriptor <= rhs.descriptor;
}

template<typename Graph,typename Descriptor>
bool operator>=(const Gwdwg_descriptor<Graph,Descriptor>& lhs,
                const Gwdwg_descriptor<Graph,Descriptor>& rhs)
{
  CGAL_assertion( lhs.graph == rhs.graph || rhs.graph==NULL || lhs.graph==NULL);
  return lhs.descriptor >= rhs.descriptor;
}


template<typename Graph,typename Descriptor>
std::ostream& operator<<(std::ostream& os, const Gwdwg_descriptor<Graph,Descriptor>& gd)
{
  return os << gd.descriptor << "  in a " << typeid(*gd.graph).name();
}

/*!
\ingroup PkgBGLAdaptors

The class `Graph_with_descriptor_with_graph` wraps a graph into another graph in such a way that its descriptors contain a reference to the graph they come from.

For example, calling `source(edge, graph)` will trigger an assertion if `edge` does not belong to `graph`.
It is mainly used for debugging purposes.

Property forwarding
-------------------
All internal properties of the underlying graph are forwarded.

Property maps can be wrapped with `Graph_with_descriptor_with_graph_property_map`.
\tparam Graph must be a model of a `FaceListGraph` and `HalfedgeListGraph`.

\cgalModels `FaceListGraph`
\cgalModels `HalfedgeListGraph`
\cgalModels `MutableFaceGraph` if `Graph` is a model of `MutableFaceGraph`
*/

template<typename Graph_>
struct Graph_with_descriptor_with_graph
{
  typedef Graph_ Graph;
  Graph* graph;

  typedef boost::graph_traits<Graph> gt;
  typedef Gwdwg_descriptor<Graph, typename gt::vertex_descriptor> vertex_descriptor;
  typedef Gwdwg_descriptor<Graph, typename gt::halfedge_descriptor> halfedge_descriptor;
  typedef Gwdwg_descriptor<Graph, typename gt::edge_descriptor> edge_descriptor;
  typedef Gwdwg_descriptor<Graph, typename gt::face_descriptor> face_descriptor;

  Graph_with_descriptor_with_graph()
    : graph(NULL)
  {}

  Graph_with_descriptor_with_graph(Graph& graph)
    : graph(&graph)
  {}
};


template <typename Graph, typename Graph_descriptor, typename Descriptor>
struct Descriptor2Descriptor: public CGAL::unary_function<Graph_descriptor,Descriptor>
{

  Descriptor2Descriptor()
    : graph(NULL)
  {}

  Descriptor2Descriptor(Graph& graph)
    : graph(&graph)
  {}

  Descriptor
  operator()(Graph_descriptor gd) const
  {
    CGAL_assertion(graph!=NULL);
    return Descriptor(gd,*graph);
  }

  Graph* graph;
};


} // namespace CGAL


namespace boost
{

template<typename Graph>
struct graph_traits< CGAL::Graph_with_descriptor_with_graph<Graph> >
{
  typedef CGAL::Graph_with_descriptor_with_graph<Graph> G;
  typedef typename G::vertex_descriptor vertex_descriptor;
  typedef typename G::halfedge_descriptor halfedge_descriptor;
  typedef typename G::edge_descriptor edge_descriptor;
  typedef typename G::face_descriptor face_descriptor;

  typedef CGAL::Descriptor2Descriptor<Graph, typename boost::graph_traits<Graph>::vertex_descriptor,vertex_descriptor> V2V;
  typedef boost::transform_iterator<V2V,typename boost::graph_traits<Graph>::vertex_iterator> vertex_iterator;

  typedef CGAL::Descriptor2Descriptor<Graph, typename boost::graph_traits<Graph>::halfedge_descriptor,halfedge_descriptor> H2H;
  typedef boost::transform_iterator<H2H,typename boost::graph_traits<Graph>::halfedge_iterator> halfedge_iterator;

  typedef CGAL::Descriptor2Descriptor<Graph, typename boost::graph_traits<Graph>::edge_descriptor,edge_descriptor> E2E;
  typedef boost::transform_iterator<E2E,typename boost::graph_traits<Graph>::edge_iterator> edge_iterator;

  typedef CGAL::Descriptor2Descriptor<Graph, typename boost::graph_traits<Graph>::face_descriptor,face_descriptor> F2F;
  typedef boost::transform_iterator<F2F,typename boost::graph_traits<Graph>::face_iterator> face_iterator;

  typedef boost::transform_iterator<E2E,typename boost::graph_traits<Graph>::out_edge_iterator> out_edge_iterator;

  typedef boost::transform_iterator<E2E,typename boost::graph_traits<Graph>::in_edge_iterator> in_edge_iterator;

  typedef boost::graph_traits<Graph> BGTG;
  typedef typename BGTG::directed_category directed_category;
  typedef typename BGTG::edge_parallel_category edge_parallel_category;
  typedef typename BGTG::traversal_category traversal_category;
  typedef typename BGTG::vertices_size_type vertices_size_type;
  typedef typename BGTG::edges_size_type edges_size_type;
  typedef typename BGTG::halfedges_size_type halfedges_size_type;
  typedef typename BGTG::faces_size_type faces_size_type;
  typedef typename BGTG::degree_size_type degree_size_type;

  static vertex_descriptor null_vertex()
  {
    return vertex_descriptor(BGTG::null_vertex());
  }

  static halfedge_descriptor null_halfedge()
  {
    return halfedge_descriptor(BGTG::null_halfedge());
  }

  static face_descriptor null_face()
  {
    return face_descriptor(BGTG::null_face());
  }
};

template<typename Graph>
struct graph_traits< const CGAL::Graph_with_descriptor_with_graph<Graph> >
  : public graph_traits< CGAL::Graph_with_descriptor_with_graph<Graph> >
{};


} // namespace boost


namespace CGAL {

template <typename T1, typename T2>
bool in_same_graph(const T1& t1, const T2& t2)
{
  return t1.graph == t2.graph;
}

template<typename Graph>
typename boost::graph_traits<Graph>::vertices_size_type
num_vertices(const Graph_with_descriptor_with_graph<Graph>& w)
{
  return num_vertices(*w.graph);
}

template<typename Graph>
typename boost::graph_traits<Graph>::edges_size_type
num_edges(const Graph_with_descriptor_with_graph<Graph>& w)
{
  return num_edges(*w.graph);
}

template<typename Graph>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
       const Graph_with_descriptor_with_graph<Graph>& w)
{
  CGAL_assertion(in_same_graph(v,w));
  return degree(v.descriptor, *w.graph);
}

template <class Graph>
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
           const Graph_with_descriptor_with_graph<Graph>& w)
{
  CGAL_assertion(in_same_graph(v,w));
  return out_degree(v.descriptor, *w.graph);
}

template <class Graph>
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
          const Graph_with_descriptor_with_graph<Graph>& w)
{
  CGAL_assertion(in_same_graph(v,w));
  return in_degree(v.descriptor, *w.graph);
}

template <class Graph>
typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
source(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::edge_descriptor e,
       const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  CGAL_assertion(in_same_graph(e,w));
  return vertex_descriptor(source(e.descriptor, *w.graph), *w.graph);
}

template <class Graph>
typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
target(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::edge_descriptor e,
       const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  CGAL_assertion(in_same_graph(e,w));
  return vertex_descriptor(target(e.descriptor, *w.graph), *w.graph);
}

template <class Graph>
std::pair<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::edge_descriptor, bool>
edge(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor u,
     typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
     const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits<Graph>::edge_descriptor g_edge_descriptor;
  typedef typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::edge_descriptor edge_descriptor;
  CGAL_assertion(in_same_graph(u,w));
  CGAL_assertion(in_same_graph(v,w));
  bool b;
  g_edge_descriptor ed;
  boost::tie(ed,b) = edge(u.descriptor, v.descriptor, *w.graph);
  return std::make_pair(edge_descriptor(ed,*w.graph),b);
}


template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_iterator>
vertices(const Graph_with_descriptor_with_graph<Graph> & w)
{
  typename boost::graph_traits<Graph>::vertex_iterator b,e;
  boost::tie(b,e) = vertices(*w.graph);
  return std::make_pair(boost::make_transform_iterator(b,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::V2V(*w.graph)),
                        boost::make_transform_iterator(e,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::V2V(*w.graph)));
}

template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::edge_iterator>
edges(const Graph_with_descriptor_with_graph<Graph> & w)
{
  typename boost::graph_traits<Graph>::edge_iterator b,e;
  boost::tie(b,e) = edges(*w.graph);
  return std::make_pair(boost::make_transform_iterator(b,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)),
                        boost::make_transform_iterator(e,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)));
}

template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
          const Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(v,w));
  typename boost::graph_traits<Graph>::out_edge_iterator b,e;
  boost::tie(b,e) = out_edges(v.descriptor, *w.graph);
  return std::make_pair(boost::make_transform_iterator(b,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)),
                        boost::make_transform_iterator(e,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)));
}

template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::in_edge_iterator>
in_edges(typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(v,w));
  typename boost::graph_traits<Graph>::in_edge_iterator b,e;
  boost::tie(b,e) = in_edges(v.descriptor, *w.graph);
  return std::make_pair(boost::make_transform_iterator(b,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)),
                        boost::make_transform_iterator(e,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::E2E(*w.graph)));
}


//
// MutableHalfedgeGraph
//
template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
add_vertex(Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  return vertex_descriptor(add_vertex(*w.graph),*w.graph);
}


template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
add_vertex(const typename boost::graph_traits<Graph >::vertex_property_type& p,
           Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  return vertex_descriptor(add_vertex(p,*w.graph),*w.graph);
}


template <class Graph>
void
remove_vertex(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
              Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(v,w));
  remove_vertex(v.descriptor, *w.graph);
}

template<typename Graph>
void
reserve(Graph_with_descriptor_with_graph<Graph>& w,
        typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertices_size_type nv,
        typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edges_size_type ne,
        typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::faces_size_type nf)
{
  reserve(*w.graph, nv, ne, nf);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edge_descriptor
add_edge(Graph_with_descriptor_with_graph<Graph> & w)
{
    typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edge_descriptor edge_descriptor;
    return edge_descriptor(add_edge(*w.graph),*w.graph);
}

template <class Graph>
void
remove_edge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edge_descriptor e,
            Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(e,w));
  remove_edge(e.descriptor, *w.graph);
}

template <class Graph>
void
set_target(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h1,
           typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
           Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(h1,w));
  CGAL_assertion(in_same_graph(v,w));
  set_target(h1.descriptor, v.descriptor, *w.graph);
}

template <class Graph>
void
set_next(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h1,
         typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h2,
         Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(h1,w));
  CGAL_assertion(in_same_graph(h2,w));
  set_next(h1.descriptor, h2.descriptor, *w.graph);
}

//
// MutableFaceGraph
//
template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor
add_face(Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor face_descriptor;
  return face_descriptor(add_face(*w.graph),*w.graph);
}


template <class InputIterator, class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor
add_face(InputIterator begin, InputIterator end,
         Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor G_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor G_face_descriptor;
  std::vector<G_vertex_descriptor> vertices;
  for(; begin != end; ++begin){
    typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vd = *begin;
    CGAL_assertion(in_same_graph(vd,w));
    vertices.push_back(vd.descriptor);
  }
  G_face_descriptor fd = add_face(vertices.begin(), vertices.end(), *w.graph);
  return typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::face_descriptor(fd,*w.graph);
}


template <class Graph>
void
remove_face(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor f,
            Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(f,w));
  remove_face(f.descriptor, *w.graph);
}

template <class Graph>
void
set_face(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
         typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor f,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(h,w));
  CGAL_assertion(f==boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::null_face() || in_same_graph(f,w));
  set_face(h.descriptor, f.descriptor, *w.graph);
}

template <class Graph>
void
set_halfedge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor f,
             typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
             Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(f,w));
  CGAL_assertion(in_same_graph(h,w));
  set_halfedge(f.descriptor, h.descriptor, *w.graph);
}

template <class Graph>
void
set_halfedge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
             typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
             const Graph_with_descriptor_with_graph<Graph> & w)
{
  CGAL_assertion(in_same_graph(v,w));
  CGAL_assertion(in_same_graph(h,w));
  set_halfedge(v.descriptor, h.descriptor, *w.graph);
}

//
// HalfedgeGraph
//
template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edge_descriptor
edge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
     const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::edge_descriptor edge_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return edge_descriptor(edge(h.descriptor, *w.graph), *w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
halfedge(typename boost::graph_traits<  Graph_with_descriptor_with_graph<Graph> >::edge_descriptor e,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(e,w));
  return halfedge_descriptor(halfedge(e.descriptor, *w.graph), *w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(v,w));
  return halfedge_descriptor(halfedge(v.descriptor, *w.graph), *w.graph);
}


template <class Graph>
std::pair<typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor, bool>
halfedge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor u,
         typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor v,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  typename boost::graph_traits< Graph >::halfedge_descriptor hd;
  bool b;
  CGAL_assertion(in_same_graph(u,w));
  CGAL_assertion(in_same_graph(v,w));
  boost::tie(hd,b) = halfedge(u.descriptor, v.descriptor, *w.graph);
  return std::make_pair(halfedge_descriptor(hd,*w.graph),b);
}


template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
opposite(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return halfedge_descriptor(opposite(h.descriptor,*w.graph),*w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
source(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
       const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return vertex_descriptor(source(h.descriptor,*w.graph),*w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor
target(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
       const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::vertex_descriptor vertex_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return vertex_descriptor(target(h.descriptor,*w.graph),*w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
next(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
     const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return halfedge_descriptor(next(h.descriptor,*w.graph),*w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
prev(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
     const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return halfedge_descriptor(prev(h.descriptor,*w.graph),*w.graph);
}

//
// HalfedgeListGraph
//

template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::halfedge_iterator>
halfedges(const Graph_with_descriptor_with_graph<Graph> & w)
{
  typename boost::graph_traits<Graph>::halfedge_iterator b,e;
  boost::tie(b,e) = halfedges(*w.graph);
  return std::make_pair(boost::make_transform_iterator(b, typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::H2H(*w.graph)),
                        boost::make_transform_iterator(e, typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::H2H(*w.graph)));
}


template <class Graph>
typename boost::graph_traits<Graph>::halfedges_size_type
num_halfedges(const Graph_with_descriptor_with_graph<Graph> & w)
{
  return num_halfedges(*w.graph);
}

// FaceGraph
template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor
face(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor h,
     const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor face_descriptor;
  CGAL_assertion(in_same_graph(h,w));
  return face_descriptor(face(h.descriptor,*w.graph),*w.graph);
}

template <class Graph>
typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::face_descriptor f,
         const Graph_with_descriptor_with_graph<Graph> & w)
{
  typedef typename boost::graph_traits< Graph_with_descriptor_with_graph<Graph> >::halfedge_descriptor halfedge_descriptor;
  CGAL_assertion(in_same_graph(f,w));
  return halfedge_descriptor(halfedge(f.descriptor,*w.graph),*w.graph);
}


template <class Graph>
CGAL::Iterator_range<typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::face_iterator>
faces(const Graph_with_descriptor_with_graph<Graph> & w)
{
  typename boost::graph_traits<Graph>::face_iterator b,e;
  boost::tie(b,e) = faces(*w.graph);
  return std::make_pair(boost::make_transform_iterator(b,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::F2F(*w.graph)),
                        boost::make_transform_iterator(e,typename boost::graph_traits<Graph_with_descriptor_with_graph<Graph> >::F2F(*w.graph)));
}



template <class Graph>
typename boost::graph_traits<Graph>::vertices_size_type
num_faces(const Graph_with_descriptor_with_graph<Graph> & w)
{
  return num_faces(*w.graph);
}

template <class Graph>
bool
is_valid(const Graph_with_descriptor_with_graph<Graph> & w, bool verbose = false)
{
  return is_valid(*w.graph,verbose);
}


/*!
  \ingroup PkgBGLAdaptors
    `Graph_with_descriptor_with_graph_property_map` enables to forward properties from a
     `Graph` to a `Graph_with_descriptor_with_graph`.
    \cgalModels `Graph_with_descriptor_with_graph_property_map` the same property map concept as `PM`
    @tparam Graph a model of the `FaceListGraph` and `HalfedgeListGraph` concepts.
    @tparam PM a property_map of a `Graph`.

*/
template <typename Graph, typename PM, typename Category =
          typename boost::property_traits<PM>::category>
struct Graph_with_descriptor_with_graph_property_map {

  typedef Category category;
  typedef typename boost::property_traits<PM>::value_type value_type;
  typedef typename boost::property_traits<PM>::reference reference;
  typedef Gwdwg_descriptor<Graph, typename boost::property_traits<PM>::key_type > key_type;

  Graph* graph;
  PM pm;

  Graph_with_descriptor_with_graph_property_map()
    : graph(NULL)
  {}

  Graph_with_descriptor_with_graph_property_map(const Graph& graph, const PM& pm)
    : graph(const_cast<Graph*>(&graph)), pm(pm)
  {}

  template <typename Descriptor>
  friend
  reference
  get(const Graph_with_descriptor_with_graph_property_map<Graph,PM>& gpm, const Descriptor& d)
  {
    CGAL_assertion(gpm.graph!=NULL);
    CGAL_assertion(d.graph == gpm.graph);
    return get(gpm.pm, d.descriptor);
  }

  template <typename Descriptor>
  friend
  void
  put(const Graph_with_descriptor_with_graph_property_map<Graph,PM>& gpm, const Descriptor& d,   const value_type& v)
  {
    CGAL_assertion(gpm.graph!=NULL);
    CGAL_assertion(d.graph == gpm.graph);
    put(gpm.pm, d.descriptor, v);
  }
}; // class Graph_with_descriptor_with_graph_property_map

//specialisation for lvaluepropertymaps
template <typename Graph, typename PM>
struct Graph_with_descriptor_with_graph_property_map<Graph, PM, boost::lvalue_property_map_tag> {

  typedef boost::lvalue_property_map_tag category;
  typedef typename boost::property_traits<PM>::value_type value_type;
  typedef typename boost::property_traits<PM>::reference reference;
  typedef Gwdwg_descriptor<Graph, typename boost::property_traits<PM>::key_type > key_type;

  Graph* graph;
  PM pm;

  value_type& operator[](key_type& k) const
  {
      return get(*this, k);
  }

  Graph_with_descriptor_with_graph_property_map()
    : graph(NULL)
  {}

  Graph_with_descriptor_with_graph_property_map(const Graph& graph, const PM& pm)
    : graph(const_cast<Graph*>(&graph)), pm(pm)
  {}

  template <typename Descriptor>
  friend
  reference
  get(const Graph_with_descriptor_with_graph_property_map<Graph,PM>& gpm, const Descriptor& d)
  {
    CGAL_assertion(gpm.graph!=NULL);
    CGAL_assertion(d.graph == gpm.graph);
    return get(gpm.pm, d.descriptor);
  }

  template <typename Descriptor>
  friend
  void
  put(const Graph_with_descriptor_with_graph_property_map<Graph,PM>& gpm, const Descriptor& d,   const value_type& v)
  {
    CGAL_assertion(gpm.graph!=NULL);
    CGAL_assertion(d.graph == gpm.graph);
    put(gpm.pm, d.descriptor, v);
  }
}; // class Graph_with_descriptor_with_graph_property_map



template <class Graph, class PropertyTag>
Graph_with_descriptor_with_graph_property_map<Graph, typename boost::property_map<Graph, PropertyTag >::type>
get(PropertyTag ptag, const Graph_with_descriptor_with_graph<Graph>& w)
{
  typedef typename boost::property_map<Graph,PropertyTag>::type PM;
  typedef Graph_with_descriptor_with_graph_property_map<Graph, PM> GPM;
  return GPM(*w.graph, get(ptag,*w.graph));
}


template <class Graph, class PropertyTag, class Descriptor>
typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type
get(PropertyTag ptag, const Graph_with_descriptor_with_graph<Graph>& w, const Descriptor& d)
{
  CGAL_assertion(d.graph == w.graph);
  return get(ptag, *w.graph, d.descriptor);
}


template <class Graph, class PropertyTag, class Descriptor>
void
put(PropertyTag ptag, const Graph_with_descriptor_with_graph<Graph>& w, const Descriptor& d, typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type& v)
{
  CGAL_assertion(d.graph == w.graph);
  put(ptag, *w.graph, d.descriptor, v);
}



template <class G, class D>
std::size_t hash_value(CGAL::Gwdwg_descriptor<G,D> d)
{
  return hash_value(d.descriptor);
}

template<typename Graph, typename PropertyTag>
struct graph_has_property<CGAL::Graph_with_descriptor_with_graph<Graph>, PropertyTag>
  : graph_has_property<Graph, PropertyTag> {};
}//end namespace CGAL

namespace boost {
  template <typename Graph, typename PropertyTag>
  struct property_map<CGAL::Graph_with_descriptor_with_graph<Graph>,PropertyTag> {
    typedef CGAL::Graph_with_descriptor_with_graph_property_map<Graph, typename boost::property_map<Graph, PropertyTag >::type> type;
    typedef CGAL::Graph_with_descriptor_with_graph_property_map<Graph, typename boost::property_map<Graph, PropertyTag >::const_type> const_type;
  };


}// namespace boost

#endif //CGAL_BOOST_GRAPH_GWDWG_H
