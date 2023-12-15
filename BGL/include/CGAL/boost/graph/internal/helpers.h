// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Andread Fabri

#ifndef CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
#define CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H

#include <CGAL/iterator.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/value_type_traits.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <tuple>

namespace CGAL {

// breaks a dependency loop between <CGAL/boost/graph/helpers.h>
// and <CGAL/boost/graph/iterator.h>
template <typename Graph> class Halfedge_around_target_iterator;

namespace internal {

template <typename Graph>
void
set_border(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  set_face(h, boost::graph_traits<Graph>::null_face(), g);
}

template <typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
copy(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{
  typename boost::graph_traits<Graph>::edge_descriptor e = add_edge(g);
  typename boost::graph_traits<Graph>::halfedge_descriptor res = halfedge(e,g);
  typename boost::graph_traits<Graph>::halfedge_descriptor ropp = opposite(res, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor hopp = opposite(h, g);
  set_target(res, target(h, g), g);
  set_target(hopp, target(hopp, g), g);
  set_face(res, face(h, g), g);
  set_face(ropp, face(hopp, g), g);
  // note that we cannot call set_next as it then would call set_prev on the  original
  return res;
 }


template <typename Graph>
void
set_vertex_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{ set_halfedge(target(h, g), h, g); }


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  // makes `opposite(h,g)' the successor of h.
  set_next( h, opposite(h, g), g);
}


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , typename boost::graph_traits<Graph>::vertex_descriptor const& v
          , Graph& g)
{
  // makes `h->opposite()' the successor of h and sets the incident
  // vertex of h to v.
  set_next(h, opposite(h, g), g);
  set_target(h, v, g);
  set_halfedge(v, h, g);
}

template <typename Graph>
void
insert_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , typename boost::graph_traits<Graph>::halfedge_descriptor const& h2
           , Graph& g)
{
  set_next(h, next(h2,g), g);
  set_next(h2, opposite(h, g), g);
  set_target(h, target(h2, g), g);
}


template <typename Graph>
void
remove_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , Graph& g)
{
  set_next(h, next(opposite(next(h, g), g), g), g);
}


template <typename Graph>
void
set_face_in_face_loop(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                      typename boost::graph_traits<Graph>::face_descriptor f,
                      Graph& g)
{
  typename boost::graph_traits<Graph>::halfedge_descriptor end = h;
  do {
    set_face(h, f, g);
    h = next(h, g);
  } while ( h != end);
}


template <typename Graph>
void insert_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                     , typename boost::graph_traits<Graph>::halfedge_descriptor const& f
                     , Graph& g)
{
  set_next(h, next(f, g), g);
  set_next(f, h, g);
  set_face(h, face(f, g), g);
}

template <typename Graph>
std::size_t
exact_num_vertices(const Graph& g)
{
  typename boost::graph_traits<Graph>::vertex_iterator beg, end;
  std::tie(beg,end) = vertices(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_halfedges(const Graph& g)
{
  typename boost::graph_traits<Graph>::halfedge_iterator beg, end;
  std::tie(beg,end) = halfedges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_edges(const Graph& g)
{
  typename boost::graph_traits<Graph>::edge_iterator beg, end;
  std::tie(beg,end) = edges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_faces(const Graph& g)
{
  typename boost::graph_traits<Graph>::face_iterator beg, end;
  std::tie(beg,end) = faces(g);
  return std::distance(beg,end);
}

template<typename Graph>
bool
is_isolated(typename boost::graph_traits<Graph>::vertex_descriptor v,
            const Graph& g)
{
  return halfedge(v, g) == boost::graph_traits<Graph>::null_halfedge();
}

template<typename Graph>
void
adjust_incoming_halfedge(typename boost::graph_traits<Graph>::vertex_descriptor v,
                         Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(v, g);
  halfedge_descriptor hh = h;
  if (h != boost::graph_traits<Graph>::null_halfedge())
  {
    if (target(h, g) != v)
    {
      // wrong target, flip
      h = opposite(h, g);
      hh = h;
      set_halfedge(v, h, g);
    }
    do
    {
      if(face(h, g)==boost::graph_traits<Graph>::null_face())
      {
        set_halfedge(v, h, g);
        return;
      }
      h = opposite(next(h, g), g);
    } while (h != hh);
  }
}


} // internal

namespace impl
{
  template<typename PMAP>
  struct Output_iterator_functor
  {
    typedef typename boost::property_traits<PMAP>::key_type input_t;
    typedef typename boost::property_traits<PMAP>::value_type output_t;
    PMAP map;
    Output_iterator_functor(PMAP map)
      :map(map)
    {
    }
    void operator()(const typename std::pair<input_t, output_t>& pair)
    {
      put(map, pair.first, pair.second);
    }
  };

  template<typename PMAP>
  boost::function_output_iterator<Output_iterator_functor<PMAP> > make_functor(PMAP map)
  {
    return boost::make_function_output_iterator(Output_iterator_functor<PMAP>(map));
  }

  inline Emptyset_iterator make_functor(const internal_np::Param_not_found&)
  {
    return Emptyset_iterator();
  }
}//end of impl

template <class PMAP>
struct value_type_traits<boost::function_output_iterator<impl::Output_iterator_functor<PMAP>>>
{
  typedef std::pair<typename impl::Output_iterator_functor<PMAP>::input_t,
                    typename impl::Output_iterator_functor<PMAP>::output_t> type;
};

} // CGAL


#endif // CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
