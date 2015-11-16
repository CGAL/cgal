// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Jane Tournois

#ifndef CGAL_BOOST_GRAPH_VISITOR_H
#define CGAL_BOOST_GRAPH_VISITOR_H

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/ref.hpp>

namespace boost
{
  template<typename V, typename Graph>
  struct graph_traits<
    boost::tuple<boost::reference_wrapper<V>, boost::reference_wrapper<Graph> > >
    : boost::graph_traits< Graph >
  {
    typedef boost::graph_traits<Graph> Base;
    typedef typename Base::vertex_property_type vertex_property_type;
  };

  template<typename V, typename Graph>
  struct graph_traits<
    boost::tuple<boost::reference_wrapper<V>, boost::reference_wrapper<Graph> > const >
    : boost::graph_traits< Graph >
  {};

  template<typename V, typename Graph, class PropertyTag>
  struct property_map<
      boost::tuple<boost::reference_wrapper<V>, boost::reference_wrapper<Graph> >,
      PropertyTag>
      : public property_map<Graph, PropertyTag>
  {};

  template<typename V, typename Graph, class PropertyTag>
  struct property_map<
    const boost::tuple<boost::reference_wrapper<V>, boost::reference_wrapper<Graph> >,
    PropertyTag>
    : public property_map<Graph, PropertyTag>
  {};

} // namespace boost

namespace CGAL
{
template<typename V, typename Graph>
boost::tuple<boost::reference_wrapper<V>,
             boost::reference_wrapper<Graph> >
make_graph_with_visitor(V& v, Graph& g)
{
  return boost::make_tuple(boost::ref(v), boost::ref(g));
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
num_vertices(const boost::tuple<boost::reference_wrapper<Visitor>,
                                boost::reference_wrapper<Graph> >& w)
{
  num_vertices(boost::unwrap_ref(w.get<0>()));
  return num_vertices(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::edges_size_type
num_edges(const boost::tuple<boost::reference_wrapper<Visitor>,
                             boost::reference_wrapper<Graph> >& w)
{
  num_edges(boost::unwrap_ref(w.get<0>()));
  return num_edges(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Graph>::vertex_descriptor v
     , const boost::tuple<boost::reference_wrapper<Visitor>,
                          boost::reference_wrapper<Graph> >& w)
{
  degree(v, boost::unwrap_ref(w.get<0>()));
  return degree(v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
         , const boost::tuple<boost::reference_wrapper<Visitor>,
                              boost::reference_wrapper<Graph> >& w)
{
  out_degree(v, boost::unwrap_ref(w.get<0>()));
  return out_degree(v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Graph>::vertex_descriptor v
        , const boost::tuple<boost::reference_wrapper<Visitor>,
                             boost::reference_wrapper<Graph> >& w)
{
  in_degree(v, boost::unwrap_ref(w.get<0>()));
  return in_degree(v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::vertex_descriptor
source(typename boost::graph_traits<Graph>::edge_descriptor e
     , const boost::tuple<boost::reference_wrapper<Visitor>,
                          boost::reference_wrapper<Graph> > & w)
{
  source(e, boost::unwrap_ref(w.get<0>()));
  return source(e, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::vertex_descriptor
target(typename boost::graph_traits<Graph>::edge_descriptor e
     , const boost::tuple<boost::reference_wrapper<Visitor>,
                          boost::reference_wrapper<Graph> > & w)
{
  target(e, boost::unwrap_ref(w.get<0>()));
  return target(e, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
edge(typename boost::graph_traits<Graph>::vertex_descriptor u
   , typename boost::graph_traits<Graph>::vertex_descriptor v
   , const boost::tuple<boost::reference_wrapper<Visitor>,
                        boost::reference_wrapper<Graph> > & w)
{
  edge(u, v, boost::unwrap_ref(w.get<0>()));
  return edge(u, v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::vertex_iterator>
vertices(const boost::tuple<boost::reference_wrapper<Visitor>,
                            boost::reference_wrapper<Graph> >& w)
{
  vertices(boost::unwrap_ref(w.get<0>()));
  return vertices(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::edge_iterator>
edges(const boost::tuple<boost::reference_wrapper<Visitor>,
                         boost::reference_wrapper<Graph> >& w)
{
  edges(boost::unwrap_ref(w.get<0>()));
  return edges(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::in_edge_iterator>
in_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
       , const boost::tuple<boost::reference_wrapper<Visitor>,
                            boost::reference_wrapper<Graph> >& w)
{
  in_edges(u, boost::unwrap_ref(w.get<0>()));
  return in_edges(u, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::out_edge_iterator>
out_edges(typename boost::graph_traits<Graph>::vertex_descriptor u
        , const boost::tuple<boost::reference_wrapper<Visitor>,
                             boost::reference_wrapper<Graph> >& w)
{
  out_edges(u, boost::unwrap_ref(w.get<0>()));
  return out_edges(u, boost::unwrap_ref(w.get<1>()));
}

//
// MutableHalfedgeGraph
// 

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
add_vertex(boost::tuple<boost::reference_wrapper<Visitor>,
                        boost::reference_wrapper<Graph> >& w)
{
  add_vertex(boost::unwrap_ref(w.get<0>()));
  return add_vertex(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
add_vertex(const typename boost::graph_traits<Graph >::vertex_property_type& p
         , boost::tuple<boost::reference_wrapper<Visitor>,
                        boost::reference_wrapper<Graph> >& w)
{
  add_vertex(p, boost::unwrap_ref(w.get<0>()));
  return add_vertex(p, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
remove_vertex(typename boost::graph_traits< Graph >::vertex_descriptor v
            , boost::tuple<boost::reference_wrapper<Visitor>,
                           boost::reference_wrapper<Graph> >& w)
{
  remove_vertex(v, boost::unwrap_ref(w.get<0>()));
  remove_vertex(v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::edge_descriptor
add_edge(boost::tuple<boost::reference_wrapper<Visitor>,
                      boost::reference_wrapper<Graph> >& w)
{
  add_edge(boost::unwrap_ref(w.get<0>()));
  return add_edge(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
remove_edge(typename boost::graph_traits< Graph >::edge_descriptor e
, boost::tuple<boost::reference_wrapper<Visitor>,
               boost::reference_wrapper<Graph> >& w)
{
  remove_edge(e, boost::unwrap_ref(w.get<0>()));
  remove_edge(e, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
set_target(typename boost::graph_traits< Graph >::halfedge_descriptor h1
, typename boost::graph_traits< Graph >::vertex_descriptor v
, boost::tuple<boost::reference_wrapper<Visitor>,
               boost::reference_wrapper<Graph> >& w)
{
  set_target(h1, v, boost::unwrap_ref(w.get<0>()));
  set_target(h1, v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
set_next(typename boost::graph_traits< Graph >::halfedge_descriptor h1
       , typename boost::graph_traits< Graph >::halfedge_descriptor h2
       , boost::tuple<boost::reference_wrapper<Visitor >,
                      boost::reference_wrapper<Graph> >& w)
{
  set_next(h1, h2, boost::unwrap_ref(w.get<0>()));
  set_next(h1, h2, boost::unwrap_ref(w.get<1>()));
}

//
// MutableFaceGraph 
//
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
add_face(boost::tuple<boost::reference_wrapper<Visitor>,
  boost::reference_wrapper<Graph> >& w)
{
  add_face(boost::unwrap_ref(w.get<0>()));
  return add_face(boost::unwrap_ref(w.get<1>()));
}

template <class InputIterator, class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
add_face(InputIterator begin,
         InputIterator end,
         boost::tuple<boost::reference_wrapper<Visitor>,
         boost::reference_wrapper<Graph> >& w)
{
  add_face(begin, end, boost::unwrap_ref(w.get<0>()));
  return add_face(begin, end, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
remove_face(typename boost::graph_traits< Graph >::face_descriptor f
, boost::tuple<boost::reference_wrapper<Visitor>,
boost::reference_wrapper<Graph> >& w)
{
  remove_face(f, boost::unwrap_ref(w.get<0>()));
  return remove_face(f, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
set_face(typename boost::graph_traits< Graph >::halfedge_descriptor h
, typename boost::graph_traits< Graph >::face_descriptor f
, const boost::tuple<boost::reference_wrapper<Visitor>,
boost::reference_wrapper<Graph> >& w)
{
  set_face(h, f, boost::unwrap_ref(w.get<0>()));
  set_face(h, f, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
set_halfedge(typename boost::graph_traits< Graph >::face_descriptor f
, typename boost::graph_traits< Graph >::halfedge_descriptor h
, boost::tuple<boost::reference_wrapper<Visitor>,
boost::reference_wrapper<Graph> >& w)
{
  set_halfedge(f, h, boost::unwrap_ref(w.get<0>()));
  set_halfedge(f, h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
void
set_halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
, typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>,
boost::reference_wrapper<Graph> >& w)
{
  set_halfedge(v, h, boost::unwrap_ref(w.get<0>()));
  set_halfedge(v, h, boost::unwrap_ref(w.get<1>()));
}

//
// HalfedgeGraph
//
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::edge_descriptor
edge(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>,
                     boost::reference_wrapper<Graph> >& w)
{
  edge(h, boost::unwrap_ref(w.get<0>()));
  return edge(h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::edge_descriptor e
, const boost::tuple<boost::reference_wrapper<Visitor>,
                     boost::reference_wrapper<Graph> >& w)
{
  halfedge(e, boost::unwrap_ref(w.get<0>()));
  return halfedge(e, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::vertex_descriptor v
, const boost::tuple<boost::reference_wrapper<Visitor>,
                     boost::reference_wrapper<Graph> >& w)
{
  halfedge(v, boost::unwrap_ref(w.get<0>()));
  return halfedge(v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
std::pair< typename boost::graph_traits< Graph >::halfedge_descriptor
  , bool>
  halfedge(typename boost::graph_traits< Graph >::vertex_descriptor u
  , typename boost::graph_traits< Graph >::vertex_descriptor v
  , const boost::tuple<boost::reference_wrapper<Visitor>,
                       boost::reference_wrapper<Graph> >& w)
{
  halfedge(u, v, boost::unwrap_ref(w.get<0>()));
  return halfedge(u, v, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
opposite(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>,
                      boost::reference_wrapper<Graph> >& w)
{
  opposite(h, boost::unwrap_ref(w.get<0>()));
  return opposite(h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
source(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>,
                     boost::reference_wrapper<Graph> >& w)
{
  source(h, boost::unwrap_ref(w.get<0>()));
  return source(h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::vertex_descriptor
target(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  target(h, boost::unwrap_ref(w.get<0>()));
  return target(h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
next(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
, const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  next(outedge, boost::unwrap_ref(w.get<0>()));
  return next(outedge, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
prev(typename boost::graph_traits< Graph >::halfedge_descriptor outedge
, const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  prev(outedge, boost::unwrap_ref(w.get<0>()));
  return prev(outedge, boost::unwrap_ref(w.get<1>()));
}

//
// HalfedgeListGraph
//
template <class Graph, class Visitor>
CGAL::Iterator_range<typename boost::graph_traits< Graph >::halfedge_iterator>
halfedges(const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  halfedges(boost::unwrap_ref(w.get<0>()));
  return halfedges(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedges_size_type
num_halfedges(const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  num_halfedges(boost::unwrap_ref(w.get<0>()));
  return num_halfedges(boost::unwrap_ref(w.get<1>()));
}

// Graph
template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::face_descriptor
face(typename boost::graph_traits< Graph >::halfedge_descriptor h
, const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  face(h, boost::unwrap_ref(w.get<0>()));
  return face(h, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits< Graph >::halfedge_descriptor
halfedge(typename boost::graph_traits< Graph >::face_descriptor f
, const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  halfedge(f, boost::unwrap_ref(w.get<0>()));
  return halfedge(f, boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
inline CGAL::Iterator_range<typename boost::graph_traits<Graph>::face_iterator >
faces(const boost::tuple<boost::reference_wrapper<Visitor>, boost::reference_wrapper<Graph> >& w)
{
  faces(boost::unwrap_ref(w.get<0>()));
  return faces(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
typename boost::graph_traits<Graph>::faces_size_type
num_faces(const boost::tuple<boost::reference_wrapper<Visitor>,
                             boost::reference_wrapper<Graph> >& w)
{
  num_faces(boost::unwrap_ref(w.get<0>()));
  return num_faces(boost::unwrap_ref(w.get<1>()));
}

template <class Graph, class Visitor>
bool is_valid(const boost::tuple<boost::reference_wrapper<Visitor>,
                                 boost::reference_wrapper<Graph> >& w
              , bool verbose = false)
{
  is_valid(boost::unwrap_ref(w.get<0>()), verbose);
  return is_valid(boost::unwrap_ref(w.get<1>()), verbose);
}

template <class Graph, class PropertyTag, class Visitor>
typename boost::property_map< Graph, PropertyTag >::type
get(PropertyTag ptag,
    const boost::tuple<boost::reference_wrapper<Visitor>,
                       boost::reference_wrapper<Graph> >& w)
{
  get(ptag, boost::unwrap_ref(w.get<0>()));
  return get(ptag, boost::unwrap_ref(w.get<1>()));
}

}//end namespace CGAL

#endif //CGAL_BOOST_GRAPH_VISITOR_H
