// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Maxime Gimeno
#ifndef CGAL_BOOST_GRAPH_CONNECTED_COMPONENTS_GRAPH_H
#define CGAL_BOOST_GRAPH_CONNECTED_COMPONENTS_GRAPH_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/assertions.h>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/dynamic_bitset.hpp>

namespace CGAL
{

  /*!
   * \ingroup PkgBGLHelper
   *
   * The class `Connected_components_graph` wraps a graph into another graph in such a way that only the specified connected components are seen from the outside.
   *
   * For example, calling `vertices(graph)` will return an iterator range of all but only the vertices that belong to the selected connected components.
   *
   * \tparam Graph must be a model of a `FaceListGraph` and `HalfedgeListGraph`.
   * \tparam FIMap a model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor` as key and `graph_traits<Graph>::%faces_size_type` as value
   * \tparam VIMap a model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor` as key and `graph_traits<Graph>::%vertices_size_type` as value
   * \tparam HIMap a model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%halfedge_descriptor` as key and `graph_traits<Graph>::%halfedges_size_type` as value
   *
   * \cgalModels `FaceListGraph`
   * \cgalModels `HalfedgeListGraph`
   */
template<typename Graph,
         typename FIMap = typename boost::property_map<Graph, CGAL::face_index_t>::type,
         typename VIMap = typename boost::property_map<Graph, boost::vertex_index_t>::type,
         typename HIMap = typename boost::property_map<Graph, CGAL::halfedge_index_t>::type>
struct Connected_components_graph
{
  typedef boost::graph_traits<Graph>                  gt;
  typedef typename gt::vertex_descriptor              vertex_descriptor;
  typedef typename gt::halfedge_descriptor            halfedge_descriptor;
  typedef typename gt::edge_descriptor                edge_descriptor;
  typedef typename gt::face_descriptor                face_descriptor;
  typedef typename boost::property_traits< FIMap >::value_type face_index_type;
  typedef typename boost::property_traits< VIMap >::value_type vertex_index_type;
  typedef typename boost::property_traits< HIMap >::value_type halfedge_index_type;
  typedef Connected_components_graph<Graph, FIMap, VIMap, HIMap>   Self;

  template <typename FaceComponentMap, class IndexRangeIterator>
  void base_iterator_constructor(IndexRangeIterator begin,
                                 IndexRangeIterator end,
                                 FaceComponentMap fcmap)
  {
    face_patch.resize(num_faces(_graph));
    vertex_patch.resize(num_vertices(_graph));
    halfedge_patch.resize(num_halfedges(_graph));
    boost::unordered_set<typename boost::property_traits<FaceComponentMap>::value_type> pids;
    for(IndexRangeIterator it = begin;
        it != end;
        ++it)
    {
      pids.insert(*it);
    }

    BOOST_FOREACH(face_descriptor fd, faces(_graph) )
    {
      if(pids.find(boost::get(fcmap, fd))!= pids.end() )
      {
        face_patch.set(get(fimap, fd));
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, _graph), _graph))
        {
          halfedge_patch.set(get(himap, hd));
          halfedge_patch.set(get(himap, opposite(hd, _graph)));
          vertex_patch.set(get(vimap, target(hd, _graph)));
        }
      }
    }
  }

  template <typename FaceComponentMap>
  void base_constructor(FaceComponentMap fcmap,
      typename boost::property_traits<FaceComponentMap>::value_type pid)
  {
    face_patch.resize(num_faces(_graph));
    vertex_patch.resize(num_vertices(_graph));
    halfedge_patch.resize(num_halfedges(_graph));
    BOOST_FOREACH(face_descriptor fd, faces(_graph) )
    {
      if(boost::get(fcmap, fd) == pid)
      {
        face_patch.set(get(fimap, fd));
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, _graph), _graph))
        {
          halfedge_patch.set(get(himap, hd));
          halfedge_patch.set(get(himap, opposite(hd, _graph)));
          vertex_patch.set(get(vimap, target(hd, _graph)));
        }
      }
    }
  }
  /*!
   * \brief Creates a Connected_components_graph of the connected components of `graph` specified in the range
   * defined by `begin` and `end`.
  * typedef unspecified_type IndexRangeIterator;

   * \tparam FaceComponentMap a model of `ReadablePropertyMap` with
      `boost::graph_traits<Graph>::%face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.

   * \param graph the graph containing the wanted patches.
   * \param fcmap the property_map that assigns a connected component to each face, with
      `boost::graph_traits<Graph>::%face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \param begin an interator to the beginning of a range of connected components indices of interest.
   * \param end an interator to the element past the end of a range of connected components indices of interest.
   * \param fimap the property map that assigns an index to each face,
   * with boost::graph_traits<Graph>::%face_descriptor as key type and boost::graph_traits<Graph>::%faces_size_type as value type
   * \param vimap the property map that assigns an index to each vertex,
   *  with boost::graph_traits<Graph>::%vertex_descriptor as key type and boost::graph_traits<Graph>::%vertices_size_type as value type
   * \param himap the property map that assigns an index to each halfedge,
   *  with boost::graph_traits<Graph>::%halfedge_descriptor as key type and boost::graph_traits<Graph>::%halfedges_size_type as value type
   */

  template <typename FaceComponentMap, class IndexRangeIterator>
  Connected_components_graph(const Graph& graph,
                             FaceComponentMap fcmap,
                             IndexRangeIterator begin,
                             IndexRangeIterator end,
                           #ifdef DOXYGEN_RUNNING
                             FIMap fimap = get(CGAL::face_index, graph),
                             VIMap vimap = get(boost::vertex_index, graph),
                             HIMap himap = get(CGAL::halfedge_index, graph)
                           #else
                             FIMap fimap,
                             VIMap vimap,
                             HIMap himap
                           #endif
                             )
    : _graph(const_cast<Graph&>(graph)), fimap(fimap), vimap(vimap), himap(himap)
  {
    base_iterator_constructor(begin, end, fcmap);
  }

  template <typename FaceComponentMap, class IndexRangeIterator>
  Connected_components_graph(const Graph& graph,
                             FaceComponentMap fcmap,
                             IndexRangeIterator begin,
                             IndexRangeIterator end)
    : _graph(const_cast<Graph&>(graph))
  {
    fimap = get(CGAL::face_index, graph);
    vimap = get(boost::vertex_index, graph);
    himap = get(CGAL::halfedge_index, graph);
    base_iterator_constructor(begin, end, fcmap);
  }
  /*!
   * \brief Creates a Connected_components_graph of the connected component `pid` of `graph`.
   *
   * \tparam FaceComponentMap a model of `ReadablePropertyMap` with
      `boost::graph_traits<Graph>::%face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \param graph the graph containing the wanted patch.
   * \param fcmap the property_map that assigns a connected component's index to each face, with
      `boost::graph_traits<Graph>::%face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \param pid the index of the connected component of interest.
   * \param fimap the property map that assigns an index to each face,
   * with boost::graph_traits<Graph>::%face_descriptor as key type and boost::graph_traits<Graph>::%faces_size_type as value type
   * \param vimap the property map that assigns an index to each vertex,
   *  with boost::graph_traits<Graph>::%vertex_descriptor as key type and boost::graph_traits<Graph>::%vertices_size_type as value type
   * \param himap the property map that assigns an index to each halfedge,
   *  with boost::graph_traits<Graph>::%halfedge_descriptor as key type and boost::graph_traits<Graph>::%halfedges_size_type as value type
   */
  template <typename FaceComponentMap>
  Connected_components_graph(const Graph& graph,
                             FaceComponentMap fcmap,
                             typename boost::property_traits<FaceComponentMap>::value_type pid,
                           #ifdef DOXYGEN_RUNNING
                             FIMap fimap = get(CGAL::face_index, graph),
                             VIMap vimap = get(boost::vertex_index, graph),
                             HIMap himap = get(CGAL::halfedge_index, graph)
                           #else
                             FIMap fimap,
                             VIMap vimap,
                             HIMap himap
                           #endif
                             )
    : _graph(const_cast<Graph&>(graph)), fimap(fimap), vimap(vimap), himap(himap)
  {
    base_constructor(fcmap, pid);
  }

  template <typename FaceComponentMap>
  Connected_components_graph(const Graph& graph,
                             FaceComponentMap fcmap,
                             typename boost::property_traits<FaceComponentMap>::value_type pid)
    : _graph(const_cast<Graph&>(graph))
  {
    fimap = get(CGAL::face_index, graph);
    vimap = get(boost::vertex_index, graph);
    himap = get(CGAL::halfedge_index, graph);
    base_constructor(fcmap, pid);
  }
  ///returns a const reference to the graph of the Connected_components_graph.
  const Graph& graph()const{ return _graph; }
  ///returns a reference to the graph of the Connected_components_graph.
  Graph& graph(){ return _graph; }

  ///change the selected component's id
  template<class FaceComponentMap>
  void set_selected_component(FaceComponentMap fcmap,
                              typename boost::property_traits<FaceComponentMap>::value_type pid)
  {
    face_indices.clear();
    vertex_indices.clear();
    halfedge_indices.clear();
    base_constructor(fcmap, pid);
  }
  /// change the selected components ids
  template<class FaceComponentMap, class IndexRangeIterator>
  void set_selected_components(FaceComponentMap fcmap,
                              IndexRangeIterator begin,
                              IndexRangeIterator end)
  {
    face_indices.clear();
    vertex_indices.clear();
    halfedge_indices.clear();
    base_iterator_constructor(begin, end, fcmap);
  }
  struct Is_simplex_valid
  {
    Is_simplex_valid(const Self* graph)
      :adapter(graph)
    {}

    Is_simplex_valid()
      :adapter(NULL)
    {}
    template<typename Simplex>
    bool operator()(Simplex s)
    {
      CGAL_assertion(adapter!=NULL);
      return (in_CC(s, *adapter));
    }
    const Self* adapter;
  };

  bool is_in_cc(face_descriptor f) const
  {
    return face_patch[get(fimap, f)];
  }

  bool is_in_cc(vertex_descriptor v) const
  {
    return vertex_patch[get(vimap, v)];
  }

  bool is_in_cc(halfedge_descriptor h) const
  {
    return halfedge_patch[get(himap, h)];
  }
  ///returns the number of faces contained in the specified connected components
  typename boost::graph_traits<Graph>::
  faces_size_type number_of_faces()const
  {
    return face_patch.count();
  }
///returns the number of vertices contained in the specified connected components
  typename boost::graph_traits<Graph>::
  vertices_size_type number_of_vertices()const
  {
    return vertex_patch.count();
  }
///returns the number of halfedges contained in the specified connected components
  typename boost::graph_traits<Graph>::
  halfedges_size_type number_of_halfedges()const
  {
    return halfedge_patch.count();
  }

  Property_map_binder< FIMap, typename Pointer_property_map< typename boost::property_traits< FIMap >::value_type >::type >
  get_face_index_map() const
  {
    if (face_indices.empty())
    {
      face_index_type index = 0;
      face_indices.resize(num_faces(_graph));
      for (std::size_t i=face_patch.find_first(); i < face_patch.npos; i = face_patch.find_next(i))
     {
        face_indices[i] = index++;
     }
    }
    return bind_property_maps(fimap, make_property_map(face_indices) );
  }

  Property_map_binder< VIMap, typename Pointer_property_map< typename boost::property_traits< VIMap >::value_type >::type >
  get_vertex_index_map() const
  {
    if (vertex_indices.empty())
    {
      vertex_index_type index = 0;
      vertex_indices.resize(num_vertices(_graph));
      for (std::size_t i=vertex_patch.find_first(); i < vertex_patch.npos; i = vertex_patch.find_next(i))
     {
        vertex_indices[i] = index++;
     }
    }
    return bind_property_maps(vimap, make_property_map(vertex_indices) );
  }

  Property_map_binder< HIMap, typename Pointer_property_map< typename boost::property_traits< HIMap >::value_type >::type >
  get_halfedge_index_map() const
  {
    if (halfedge_indices.empty())
    {
      halfedge_index_type index = 0;
      halfedge_indices.resize(num_halfedges(_graph));
      for (std::size_t i=halfedge_patch.find_first(); i < halfedge_patch.npos; i = halfedge_patch.find_next(i))
     {
        halfedge_indices[i] = index++;
     }
    }
    return bind_property_maps(himap, make_property_map(halfedge_indices) );
  }

private:
  Graph& _graph;
  FIMap fimap;
  VIMap vimap;
  HIMap himap;
  boost::dynamic_bitset<> face_patch;
  boost::dynamic_bitset<> vertex_patch;
  boost::dynamic_bitset<> halfedge_patch;
  mutable std::vector<face_index_type> face_indices;
  mutable std::vector<vertex_index_type> vertex_indices;
  mutable std::vector<halfedge_index_type> halfedge_indices;
};

} // namespace CGAL

namespace boost
{

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >
{
  typedef CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> G;
  typedef boost::graph_traits<Graph> BGTG;
  typedef typename BGTG::vertex_descriptor vertex_descriptor;
  typedef typename BGTG::halfedge_descriptor halfedge_descriptor;
  typedef typename BGTG::edge_descriptor edge_descriptor;
  typedef typename BGTG::face_descriptor face_descriptor;

  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::vertex_iterator>    vertex_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::halfedge_iterator>  halfedge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::edge_iterator>      edge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::face_iterator>      face_iterator;

  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::out_edge_iterator>  out_edge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::in_edge_iterator>   in_edge_iterator;

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
    return BGTG::null_vertex();
  }

  static halfedge_descriptor null_halfedge()
  {
    return BGTG::null_halfedge();
  }

  static edge_descriptor null_edge()
  {
    return edge_descriptor(BGTG::null_halfedge());
  }

  static face_descriptor null_face()
  {
    return BGTG::null_face();
  }
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct graph_traits< const CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >
    : public graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >
{};


} // namespace boost


namespace CGAL {
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
bool
in_CC(const typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::face_descriptor f,
      const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.is_in_cc(f);
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
bool
in_CC(const typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
      const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return  w.is_in_cc(h);
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
bool
in_CC(const typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
      const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return  w.is_in_cc(halfedge(e, w.graph()));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
bool
in_CC(const typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
      const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.is_in_cc(v);
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::vertices_size_type
num_vertices(const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.number_of_vertices();
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::edges_size_type
num_edges(const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.number_of_halfedges()/2;
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
       const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(in_CC(v, w));
  typename boost::graph_traits<Graph>::degree_size_type v_deg = 0;
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(v, w);
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor hcirc = h;
  do
  {
    if(in_CC(hcirc, w))
      ++v_deg;
    hcirc = opposite(next(hcirc, w.graph()), w.graph());
  }while(hcirc != h);
  return v_deg;
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
           const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(in_CC(v, w));
  return std::distance(out_edges(v, w).first ,out_edges(v, w).second);
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
          const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(in_CC(v, w));
  return std::distance(in_edges(v, w).first ,in_edges(v, w).second);
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
source(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
       const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(e, w));
  return source(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
target(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
       const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(e, w));
  return target(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor, bool>
edge(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor u,
     typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
     const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(u, w) && in_CC(v, w));
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e = edge(u, v, w.graph()).first;
  bool res = in_CC(e, w);
  return std::make_pair(e, res);
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_iterator>
vertices(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph >::vertex_iterator g_vertex_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_vertex_iterator b,e;
  boost::tie(b,e) = vertices(w.graph());
  return make_range(vertex_iterator(predicate, b, e),
                    vertex_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_iterator>
edges(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_iterator edge_iterator;
  typedef typename boost::graph_traits<Graph >::edge_iterator g_edge_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_edge_iterator b,e;
  boost::tie(b,e) = edges(w.graph());
  return make_range(edge_iterator(predicate, b, e),
                    edge_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
          const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::out_edge_iterator out_edge_iterator;
  typedef typename boost::graph_traits<Graph >::out_edge_iterator g_out_edge_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_out_edge_iterator b,e;
  boost::tie(b,e) = out_edges(v, w.graph());
  return make_range(out_edge_iterator(predicate, b, e),
                    out_edge_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::in_edge_iterator>
in_edges(typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::in_edge_iterator in_edge_iterator;
  typedef typename boost::graph_traits<Graph >::in_edge_iterator g_in_edge_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_in_edge_iterator b,e;
  boost::tie(b,e) = in_edges(v, w.graph());
  return make_range(in_edge_iterator(predicate, b, e),
                    in_edge_iterator(predicate, e, e));
}

//
// HalfedgeGraph
//
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor
edge(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(CGAL::in_CC(h, w));
  return edge(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits<  Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(CGAL::in_CC(e, w));
  return halfedge(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(v, w));
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(v, w.graph());
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor hcirc = h;
  do
  {
    if(in_CC(hcirc, w))
      return hcirc;
    hcirc = opposite(next(hcirc, w.graph()), w.graph());
  }while(hcirc != h);
  return boost::graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::null_halfedge();
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor, bool>
halfedge(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor u,
         typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(u, w) && in_CC(v, w));
  typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(u, v, w.graph()).first;
  return std::make_pair(h, in_CC(h, w));
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
opposite(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(h, w) );
     return opposite(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
source(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
       const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(h, w) );
  return source(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
target(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
       const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(h, w) );
  return target(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
next(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(in_CC(h, w));
  if(in_CC(next(h, w.graph()), w))
    return next(h, w.graph());
  //act as a border
  typedef typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h_d;
  BOOST_FOREACH( h_d hcirc,
                CGAL::halfedges_around_target(target(h, w.graph()), w.graph()))
  {
    if(hcirc != h && in_CC(hcirc, w))
    {
      return opposite(hcirc, w.graph());
    }
  }
  return boost::graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::null_halfedge();
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
prev(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  CGAL_assertion(in_CC(h, w));
  if(in_CC(prev(h, w.graph()), w))
    return prev(h, w.graph());

  //act as a border
  typedef typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h_d;
  BOOST_FOREACH(h_d hcirc,
                CGAL::halfedges_around_source(source(h, w.graph()), w.graph()))
  {
    if(hcirc != h && in_CC(hcirc, w))
    {
      return opposite(hcirc, w.graph());
    }
  }
  return boost::graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::null_halfedge();
}

//
// HalfedgeListGraph
//

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator,
typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator>
halfedges(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<Graph >::halfedge_iterator g_halfedge_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  std::pair<g_halfedge_iterator, g_halfedge_iterator> original_halfedges = halfedges(w.graph());

  return make_range(halfedge_iterator(predicate, original_halfedges.first, original_halfedges.second),
                    halfedge_iterator(predicate, original_halfedges.second, original_halfedges.second));
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::halfedges_size_type
num_halfedges(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.number_of_halfedges();
}

// FaceGraph
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::face_descriptor
face(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(CGAL::in_CC(h, w));
  if(in_CC(face(h,w.graph()), w))
    return face(h,w.graph());
  else
    return boost::graph_traits< CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::null_face();
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::face_descriptor f,
         const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(CGAL::in_CC(f, w));
  return halfedge(f,w.graph());
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
Iterator_range<typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::face_iterator>
faces(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Connected_components_graph<Graph, FIMap, VIMap, HIMap> >::face_iterator face_iterator;
  typedef typename boost::graph_traits<Graph >::face_iterator g_face_iterator;

  typename Connected_components_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  std::pair<g_face_iterator, g_face_iterator> original_faces = faces(w.graph());

  return make_range(face_iterator(predicate, original_faces.first, original_faces.second),
                    face_iterator(predicate, original_faces.second, original_faces.second));
}



template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::vertices_size_type
num_faces(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.number_of_faces();
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
bool
in_CC(const Connected_components_graph<Graph, FIMap, VIMap, HIMap> & w, bool verbose = false)
{
  return in_CC(w.graph(),verbose);
}

template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_map<Graph, PropertyTag >::const_type
get(PropertyTag ptag, const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return get(ptag, w.graph());
}

template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_map<Graph, PropertyTag >::type
get(PropertyTag ptag, Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return get(ptag, w.graph());
}

//specializations for indices
template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Connected_components_graph<Graph, FIMap, VIMap, HIMap>, CGAL::face_index_t >::type
get(CGAL::face_index_t, const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_face_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Connected_components_graph<Graph, FIMap, VIMap, HIMap>, boost::vertex_index_t >::type
get(boost::vertex_index_t, const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_vertex_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Connected_components_graph<Graph, FIMap, VIMap, HIMap>, CGAL::halfedge_index_t >::type
get(CGAL::halfedge_index_t, const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_halfedge_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type
get(PropertyTag ptag,
    const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k)
{
  return get(ptag, w.graph(), k);
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
void
put(PropertyTag ptag, const Connected_components_graph<Graph, FIMap, VIMap, HIMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k,
    typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type& v)
{
  put(ptag, w.graph(), k, v);
}

}//end namespace CGAL

namespace boost {
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap,
         typename PropertyTag>
struct property_map<CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap>,PropertyTag> {
  typedef typename boost::property_map<Graph, PropertyTag >::type type;
  typedef typename boost::property_map<Graph, PropertyTag >::const_type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap,
         typename PropertyTag>
struct graph_has_property<CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap>, PropertyTag>
    : graph_has_property<Graph, PropertyTag> {};


//specializations for indices
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct property_map<CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap>, CGAL::face_index_t>{
  typedef typename CGAL::Property_map_binder< FIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< FIMap >::value_type >::type > type;
  typedef type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct property_map<CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap>, boost::vertex_index_t>{
  typedef typename CGAL::Property_map_binder< VIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< VIMap >::value_type >::type > type;
  typedef type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>

struct property_map<CGAL::Connected_components_graph<Graph, FIMap, VIMap, HIMap>, CGAL::halfedge_index_t>{
  typedef typename CGAL::Property_map_binder< HIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< HIMap >::value_type >::type > type;
  typedef type const_type;
};
}// namespace boost
#endif // CGAL_BOOST_GRAPH_CONNECTED_COMPONENTS_GRAPH_H
