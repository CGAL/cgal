// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_BOOST_GRAPH_VISITOR_H
#define CGAL_BOOST_GRAPH_VISITOR_H

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>

namespace boost
{
  template<typename V, typename Graph>
  struct graph_traits<
    boost::tuple<std::reference_wrapper<V>, std::reference_wrapper<Graph> > >
    : boost::graph_traits< Graph >
  {
    typedef boost::graph_traits<Graph> Base;
    typedef typename Base::vertex_property_type vertex_property_type;
  };

  template<typename V, typename Graph>
  struct graph_traits<
    boost::tuple<std::reference_wrapper<V>, std::reference_wrapper<Graph> > const >
    : boost::graph_traits< Graph >
  {};

  template<typename V, typename Graph, class PropertyTag>
  struct property_map<
      boost::tuple<std::reference_wrapper<V>, std::reference_wrapper<Graph> >,
      PropertyTag>
      : public property_map<Graph, PropertyTag>
  {};

  template<typename V, typename Graph, class PropertyTag>
  struct property_map<
    const boost::tuple<std::reference_wrapper<V>, std::reference_wrapper<Graph> >,
    PropertyTag>
    : public property_map<Graph, PropertyTag>
  {};

} // namespace boost

namespace CGAL
{
template<typename V, typename Graph>
boost::tuple<std::reference_wrapper<V>,
             std::reference_wrapper<Graph> >
make_graph_with_visitor(V& v, Graph& g)
{
  return boost::make_tuple(std::ref(v), std::ref(g));
}

template<typename Graph>
class Visitor_base
{
public:
  typedef typename boost::graph_traits<Graph> gt;
  typedef typename gt::halfedge_descriptor halfedge_descriptor;
  typedef typename gt::edge_descriptor edge_descriptor;
  typedef typename gt::vertex_descriptor vertex_descriptor;
};

//// OVERLOADS FOR Visitor

template<typename Graph>
void num_vertices(Visitor_base<Graph>& w)
{}
template<typename Graph>
void num_edges(Visitor_base<Graph>& w)
{}
template<typename Graph>
void degree(typename boost::graph_traits<Graph>::vertex_descriptor v
  , const Visitor_base<Graph>& w)
{}
template <class Graph>
void out_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
  , const Visitor_base<Graph>& w)
{}
template <class Graph>
void in_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
  , const Visitor_base<Graph>& w)
{}
template <class Graph>
void source(typename boost::graph_traits<Graph>::edge_descriptor e
  , const Visitor_base<Graph> & w)
{}
template <class Graph>
void target(typename boost::graph_traits<Graph>::edge_descriptor e
  , const Visitor_base<Graph> & w)
{}
template <class Graph>
void edge(typename boost::graph_traits<Graph>::vertex_descriptor u
        , typename boost::graph_traits<Graph>::vertex_descriptor v
        , const Visitor_base<Graph> & w)
{}
template <class Graph>
inline void vertices(const Visitor_base<Graph> & w)
{}
template <class Graph>
inline void edges(const Visitor_base<Graph> & w)
{}
template <class Graph>
inline void in_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
                   , const Visitor_base<Graph> & w)
{}
template <class Graph>
inline void out_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
                    , const Visitor_base<Graph> & w)
{}

//
// MutableHalfedgeGraph
//
template <class Graph>
void add_vertex(Visitor_base<Graph> & w)
{}
template <class Graph>
void add_vertex(const typename boost::graph_traits<Graph >::vertex_property_type& p
               , Visitor_base<Graph> & w)
{}
template <class Graph>
void remove_vertex(typename boost::graph_traits< Graph >::vertex_descriptor v
                 , Visitor_base<Graph> & w)
{}
template <class Graph>
void add_edge(Visitor_base<Graph> & w)
{}
template <class Graph>
void remove_edge(typename boost::graph_traits< Graph >::edge_descriptor e
               , Visitor_base<Graph> & w)
{}
template <class Graph>
void set_target(typename boost::graph_traits< Graph >::halfedge_descriptor h1
              , typename boost::graph_traits< Graph >::vertex_descriptor v
              , Visitor_base<Graph> & w)
{}
template <class Graph>
void set_next(typename boost::graph_traits< Graph >::halfedge_descriptor h1
            , typename boost::graph_traits< Graph >::halfedge_descriptor h2
            , Visitor_base<Graph> & w)
{}

//
// MutableFaceGraph
//
template <class Graph>
void add_face(Visitor_base<Graph> & w)
{}
template <class InputIterator, class Graph>
void add_face(InputIterator begin, InputIterator end,
              Visitor_base<Graph> & w)
{}
template <class Graph>
void remove_face(typename boost::graph_traits< Graph >::face_descriptor f
               , Visitor_base<Graph> & w)
{}
template <class Graph>
void set_face(typename boost::graph_traits< Graph >::halfedge_descriptor h
            , typename boost::graph_traits< Graph >::face_descriptor f
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void set_halfedge(typename boost::graph_traits< Graph >::face_descriptor f
                , typename boost::graph_traits< Graph >::halfedge_descriptor h
                , Visitor_base<Graph> & w)
{}
template <class Graph>
void set_halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
                , typename boost::graph_traits< Graph >::halfedge_descriptor h
                , const Visitor_base<Graph> & w)
{}

//
// HalfedgeGraph
//
template <class Graph>
void edge(typename boost::graph_traits< Graph >::halfedge_descriptor h
        , const Visitor_base<Graph> & w)
{}
template <class Graph>
void halfedge(typename boost::graph_traits< Graph >::edge_descriptor e
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void halfedge(typename boost::graph_traits< Graph >::vertex_descriptor u
            , typename boost::graph_traits< Graph >::vertex_descriptor v
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void opposite(typename boost::graph_traits< Graph >::halfedge_descriptor h
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void source(typename boost::graph_traits< Graph >::halfedge_descriptor h
          , const Visitor_base<Graph> & w)
{}
template <class Graph>
void target(typename boost::graph_traits< Graph >::halfedge_descriptor h
          , const Visitor_base<Graph> & w)
{}
template <class Graph>
void next(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
        , const Visitor_base<Graph> & w)
{}
template <class Graph>
void prev(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
        , const Visitor_base<Graph> & w)
{}

//
// HalfedgeListGraph
//
template <class Graph>
void halfedges(const Visitor_base<Graph> & w)
{}
template <class Graph>
void num_halfedges(const Visitor_base<Graph> & w)
{}

// Graph
template <class Graph>
void face(typename boost::graph_traits< Graph >::halfedge_descriptor h
        , const Visitor_base<Graph> & w)
{}
template <class Graph>
void halfedge(typename boost::graph_traits< Graph >::face_descriptor f
            , const Visitor_base<Graph> & w)
{}
template <class Graph>
void faces(const Visitor_base<Graph> & w)
{}
template <class Graph>
void num_faces(const Visitor_base<Graph> & w)
{}
template <class Graph>
void is_valid(const Visitor_base<Graph> & w, bool verbose = false)
{}
template <class Graph, class PropertyTag>
void get(PropertyTag ptag, const Visitor_base<Graph>& w)
{}


//// OVERLOADS FOR TUPLE<Visitor, Graph>

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::vertices_size_type
num_vertices(const boost::tuple<std::reference_wrapper<Visitor>,
                                std::reference_wrapper<Graph> >& w)
{
  num_vertices(get<0>(w).get());
  return num_vertices(get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::edges_size_type
num_edges(const boost::tuple<std::reference_wrapper<Visitor>,
                             std::reference_wrapper<Graph> >& w)
{
  num_edges(get<0>(w).get());
  return num_edges(get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Graph>::vertex_descriptor v
     , const boost::tuple<std::reference_wrapper<Visitor>,
                          std::reference_wrapper<Graph> >& w)
{
  degree(v, get<0>(w).get());
  return degree(v, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
         , const boost::tuple<std::reference_wrapper<Visitor>,
                              std::reference_wrapper<Graph> >& w)
{
  out_degree(v, get<0>(w).get());
  return out_degree(v, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
        , const boost::tuple<std::reference_wrapper<Visitor>,
                             std::reference_wrapper<Graph> >& w)
{
  in_degree(v, get<0>(w).get());
  return in_degree(v, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::vertex_descriptor
source(typename boost::graph_traits<Graph>::edge_descriptor e
     , const boost::tuple<std::reference_wrapper<Visitor>,
                          std::reference_wrapper<Graph> > & w)
{
  source(e, get<0>(w).get());
  return source(e, get<1>(w).get);
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::vertex_descriptor
target(typename boost::graph_traits<Graph>::edge_descriptor e
     , const boost::tuple<std::reference_wrapper<Visitor>,
                          std::reference_wrapper<Graph> > & w)
{
  target(e, get<0>(w).get());
  return target(e, get<1>(w).get());
}

template <class Graph, class Visitor>
std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
edge(typename boost::graph_traits<Graph>::vertex_descriptor u
   , typename boost::graph_traits<Graph>::vertex_descriptor v
   , const boost::tuple<std::reference_wrapper<Visitor>,
                        std::reference_wrapper<Graph> > & w)
{
  edge(u, v, get<0>(w).get());
  return edge(u, v, get<1>(w).get);
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::vertex_iterator>
vertices(const boost::tuple<std::reference_wrapper<Visitor>,
                            std::reference_wrapper<Graph> >& w)
{
  vertices(get<0>(w).get());
  return vertices(get<1>(w).get());
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::edge_iterator>
edges(const boost::tuple<std::reference_wrapper<Visitor>,
                         std::reference_wrapper<Graph> >& w)
{
  edges(get<0>(w).get());
  return edges(get<1>(w).get());
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::in_edge_iterator>
in_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
       , const boost::tuple<std::reference_wrapper<Visitor>,
                            std::reference_wrapper<Graph> >& w)
{
  in_edges(u, get<0>(w).get());
  return in_edges(u, get<1>(w).get());
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::out_edge_iterator>
out_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
        , const boost::tuple<std::reference_wrapper<Visitor>,
                             std::reference_wrapper<Graph> >& w)
{
  out_edges(u, get<0>(w).get());
  return out_edges(u, get<1>(w).get());
}

//
// MutableHalfedgeGraph
//

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
add_vertex(boost::tuple<std::reference_wrapper<Visitor>,
                        std::reference_wrapper<Graph> >& w)
{
  add_vertex(get<0>(w).get());
  return add_vertex(get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
add_vertex(const typename boost::graph_traits<Graph >::vertex_property_type& p
         , boost::tuple<std::reference_wrapper<Visitor>,
                        std::reference_wrapper<Graph> >& w)
{
  add_vertex(p, get<0>(w).get());
  return add_vertex(p, get<1>(w).get());
}

template <class Graph, class Visitor>
void
remove_vertex(typename boost::graph_traits< Graph >::vertex_descriptor v
            , boost::tuple<std::reference_wrapper<Visitor>,
                           std::reference_wrapper<Graph> >& w)
{
  remove_vertex(v, get<0>(w).get());
  remove_vertex(v, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::edge_descriptor
add_edge(boost::tuple<std::reference_wrapper<Visitor>,
                      std::reference_wrapper<Graph> >& w)
{
  add_edge(get<0>(w).get());
  return add_edge(get<1>(w).get());
}

template <class Graph, class Visitor>
void
remove_edge(typename boost::graph_traits< Graph >::edge_descriptor e
, boost::tuple<std::reference_wrapper<Visitor>,
               std::reference_wrapper<Graph> >& w)
{
  remove_edge(e, get<0>(w).get());
  remove_edge(e, get<1>(w).get());
}

template <class Graph, class Visitor>
void
set_target(typename boost::graph_traits< Graph >::halfedge_descriptor h1
, typename boost::graph_traits< Graph >::vertex_descriptor v
, boost::tuple<std::reference_wrapper<Visitor>,
               std::reference_wrapper<Graph> >& w)
{
  set_target(h1, v, get<0>(w).get());
  set_target(h1, v, get<1>(w).get());
}

template <class Graph, class Visitor>
void
set_next(typename boost::graph_traits< Graph >::halfedge_descriptor h1
       , typename boost::graph_traits< Graph >::halfedge_descriptor h2
       , boost::tuple<std::reference_wrapper<Visitor >,
                      std::reference_wrapper<Graph> >& w)
{
  set_next(h1, h2, get<0>(w).get());
  set_next(h1, h2, get<1>(w).get());
}

//
// MutableFaceGraph
//
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
add_face(boost::tuple<std::reference_wrapper<Visitor>,
  std::reference_wrapper<Graph> >& w)
{
  add_face(get<0>(w).get());
  return add_face(get<1>(w).get());
}

template <class InputIterator, class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
add_face(InputIterator begin,
         InputIterator end,
         boost::tuple<std::reference_wrapper<Visitor>,
         std::reference_wrapper<Graph> >& w)
{
  add_face(begin, end, get<0>(w).get());
  return add_face(begin, end, get<1>(w).get());
}

template <class Graph, class Visitor>
void
remove_face(typename boost::graph_traits< Graph >::face_descriptor f
, boost::tuple<std::reference_wrapper<Visitor>,
std::reference_wrapper<Graph> >& w)
{
  remove_face(f, get<0>(w).get());
  return remove_face(f, get<1>(w).get());
}

template <class Graph, class Visitor>
void
set_face(typename boost::graph_traits< Graph >::halfedge_descriptor h
, typename boost::graph_traits< Graph >::face_descriptor f
, const boost::tuple<std::reference_wrapper<Visitor>,
std::reference_wrapper<Graph> >& w)
{
  set_face(h, f, get<0>(w).get());
  set_face(h, f, get<1>(w).get());
}

template <class Graph, class Visitor>
void
set_halfedge(typename boost::graph_traits< Graph >::face_descriptor f
, typename boost::graph_traits< Graph >::halfedge_descriptor h
, boost::tuple<std::reference_wrapper<Visitor>,
std::reference_wrapper<Graph> >& w)
{
  set_halfedge(f, h, get<0>(w).get());
  set_halfedge(f, h, get<1>(w).get());
}

template <class Graph, class Visitor>
void
set_halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
, typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>,
std::reference_wrapper<Graph> >& w)
{
  set_halfedge(v, h, get<0>(w).get());
  set_halfedge(v, h, get<1>(w).get());
}

//
// HalfedgeGraph
//
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::edge_descriptor
edge(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>,
                     std::reference_wrapper<Graph> >& w)
{
  edge(h, get<0>(w).get());
  return edge(h, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::edge_descriptor e
, const boost::tuple<std::reference_wrapper<Visitor>,
                     std::reference_wrapper<Graph> >& w)
{
  halfedge(e, get<0>(w).get());
  return halfedge(e, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
, const boost::tuple<std::reference_wrapper<Visitor>,
                     std::reference_wrapper<Graph> >& w)
{
  halfedge(v, get<0>(w).get());
  return halfedge(v, get<1>(w).get());
}

template <class Graph, class Visitor>
std::pair< typename boost::graph_traits< Graph >::halfedge_descriptor
  , bool>
  halfedge(typename boost::graph_traits< Graph >::vertex_descriptor u
  , typename boost::graph_traits< Graph >::vertex_descriptor v
  , const boost::tuple<std::reference_wrapper<Visitor>,
                       std::reference_wrapper<Graph> >& w)
{
  halfedge(u, v, get<0>(w).get());
  return halfedge(u, v, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
opposite(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>,
                      std::reference_wrapper<Graph> >& w)
{
  opposite(h, get<0>(w).get());
  return opposite(h, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
source(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>,
                     std::reference_wrapper<Graph> >& w)
{
  source(h, get<0>(w).get());
  return source(h, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
target(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  target(h, get<0>(w).get());
  return target(h, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
next(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
, const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  next(outedge, get<0>(w).get());
  return next(outedge, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
prev(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
, const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  prev(outedge, get<0>(w).get());
  return prev(outedge, get<1>(w).get());
}

//
// HalfedgeListGraph
//
template <class Graph, class Visitor>
CGAL::Iterator_range<typename boost::graph_traits< Graph >::halfedge_iterator>
halfedges(const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  halfedges(get<0>(w).get());
  return halfedges(get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedges_size_type
num_halfedges(const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  num_halfedges(get<0>(w).get());
  return num_halfedges(get<1>(w).get());
}

// Graph
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
face(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  face(h, get<0>(w).get());
  return face(h, get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::face_descriptor f
, const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  halfedge(f, get<0>(w).get());
  return halfedge(f, get<1>(w).get());
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::face_iterator >
faces(const boost::tuple<std::reference_wrapper<Visitor>, std::reference_wrapper<Graph> >& w)
{
  faces(get<0>(w).get());
  return faces(get<1>(w).get());
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::faces_size_type
num_faces(const boost::tuple<std::reference_wrapper<Visitor>,
                             std::reference_wrapper<Graph> >& w)
{
  num_faces(get<0>(w).get());
  return num_faces(get<1>(w).get());
}

template <class Graph, class Visitor>
bool is_valid(const boost::tuple<std::reference_wrapper<Visitor>,
                                 std::reference_wrapper<Graph> >& w
              , bool verbose = false)
{
  is_valid(get<0>(w).get(), verbose);
  return is_valid(get<1>(w).get(), verbose);
}

template <class Graph, class PropertyTag, class Visitor>
typename boost::property_map< Graph, PropertyTag >::type
get(PropertyTag ptag,
    const boost::tuple<std::reference_wrapper<Visitor>,
                       std::reference_wrapper<Graph> >& w)
{
  get(ptag, get<0>(w).get());
  return get(ptag, get<1>(w).get());
}

}//end namespace CGAL

#endif //CGAL_BOOST_GRAPH_VISITOR_H
