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
#ifndef CGAL_BOOST_GRAPH_CONNECTED_COMPONENT_GRAPH_H
#define CGAL_BOOST_GRAPH_CONNECTED_COMPONENT_GRAPH_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <boost/foreach.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/filter_iterator.hpp>

namespace CGAL
{

/*!
\ingroup PkgBGLHelper

The class `Connected_component_graph` wraps a graph  into another graph in such a way that only the specified connected component is seen from the outside.

For example, calling `vertices(graph)` will return an iterator range of all but only the vertices that belong to the connected component whose id in `fccmap` is `pid`.

\attention The functions `num_vertices()`, `num_edges()`, `num_halfedges()`, `num_faces()`, are forwarded from the underlying graph,
which means that `num_vertices(graph)` is different from `std::distance(vertices(graph).first,vertices(graph).second)`

Property maps can be wrapped with `Connected_component_graph_property_map`.
\tparam Graph must be a model of a `FaceListGraph` and `HalfedgeListGraph`.
\tparam FaceComponentMap a model of `WritablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
        `graph_traits<PolygonMesh>::faces_size_type` as value type.

\cgalModels `FaceListGraph`
\cgalModels `HalfedgeListGraph`
*/

template<typename Graph, typename FaceComponentMap>
struct Connected_component_graph
{
    typedef boost::graph_traits<Graph>                  gt;
    typedef typename gt::vertex_descriptor              vertex_descriptor;
    typedef typename gt::halfedge_descriptor            halfedge_descriptor;
    typedef typename gt::edge_descriptor                edge_descriptor;
    typedef typename gt::face_descriptor                face_descriptor;
    typedef Connected_component_graph<Graph, FaceComponentMap>   Self;


    Connected_component_graph(const Graph& graph,
                     FaceComponentMap fccmap,
                     typename boost::property_traits<FaceComponentMap>::value_type pid)
        : _graph(graph), property_map(fccmap), patch_index(pid)
    {
    }


    Connected_component_graph(const Connected_component_graph&f)
        : _graph(f.graph()), property_map(f.propertyMap()), patch_index(f.patchIndex())
    {
    }

    const Graph& graph()const{ return _graph; }

    FaceComponentMap propertyMap()const{ return property_map; }

    typename boost::property_traits<FaceComponentMap>::value_type patchIndex()const{ return patch_index; }

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
            return (is_valid(s, *adapter));
        }
        const Self* adapter;
    };

private:
    const Graph& _graph;
    FaceComponentMap property_map;
    typename boost::property_traits<FaceComponentMap>::value_type patch_index;
};

} // namespace CGAL

namespace boost
{

template<typename Graph, typename FaceComponentMap>
struct graph_traits< CGAL::Connected_component_graph<Graph, FaceComponentMap> >
{
    typedef CGAL::Connected_component_graph<Graph, FaceComponentMap> G;
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
        return vertex_descriptor(BGTG::null_vertex());
    }

    static halfedge_descriptor null_halfedge()
    {
        return halfedge_descriptor(BGTG::null_halfedge());
    }

    static edge_descriptor null_edge()
    {
        return edge_descriptor(BGTG::null_halfedge());
    }

    static face_descriptor null_face()
    {
        return face_descriptor(BGTG::null_face());
    }
};

template<typename Graph, typename FaceComponentMap >
struct graph_traits< const CGAL::Connected_component_graph<Graph, FaceComponentMap> >
        : public graph_traits< CGAL::Connected_component_graph<Graph, FaceComponentMap> >
{};


} // namespace boost


namespace CGAL {
template <class Graph, typename FaceComponentMap >
bool
is_valid(const typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::face_descriptor f,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    return boost::get(w.propertyMap(), f) == w.patchIndex();
}

template <class Graph, typename FaceComponentMap >
bool
is_valid(const typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    return is_valid(face(h, w.graph()), w) ||
            is_valid(face(opposite(h, w.graph()), w.graph()), w);
}

template <class Graph, typename FaceComponentMap >
bool
is_valid(const typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor e,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    return is_valid(halfedge(e, w.graph()), w);
}

template <class Graph, typename FaceComponentMap >
bool
is_valid(const typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h = halfedge(v, w.graph());
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor hcirc = h;
    do
    {
        if(is_valid(face(hcirc, w.graph()), w))
            return true;
        hcirc = opposite(next(hcirc, w.graph()), w.graph());
    }while(hcirc != h);
    return false;
}


template<typename Graph, typename FaceComponentMap>
typename boost::graph_traits<Graph>::vertices_size_type
num_vertices(const Connected_component_graph<Graph, FaceComponentMap>& w)
{
    return num_vertices(w.graph());
}

template<typename Graph, typename FaceComponentMap>
typename boost::graph_traits<Graph>::edges_size_type
num_edges(const Connected_component_graph<Graph, FaceComponentMap>& w)
{
    return num_edges(w.graph());
}

template<typename Graph, typename FaceComponentMap>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
       const Connected_component_graph<Graph, FaceComponentMap>& w)
{
    CGAL_assertion(is_valid(v, w));
    typename boost::graph_traits<Graph>::degree_size_type v_deg = 0;
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h = halfedge(v, w);
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor hcirc = h;
    do
    {
        if(is_valid(hcirc, w))
            ++v_deg;
        hcirc = opposite(next(hcirc, w.graph()), w.graph());
    }while(hcirc != h);
    return v_deg;
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
           const Connected_component_graph<Graph, FaceComponentMap>& w)
{
    CGAL_assertion(is_valid(v, w));
    return std::distance(out_edges(v, w).first ,out_edges(v, w).second);
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
          const Connected_component_graph<Graph, FaceComponentMap>& w)
{
    CGAL_assertion(is_valid(v, w));
    return std::distance(in_edges(v, w).first ,in_edges(v, w).second);
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor
source(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor e,
       const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(e, w));
    return source(e, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor
target(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor e,
       const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(e, w));
    return target(e, w.graph());
}

template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor, bool>
edge(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor u,
     typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
     const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(u, w) && is_valid(v, w));
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor e = edge(u, v, w.graph()).first;
    bool res = is_valid(e, w);
    return std::make_pair(e, res);
}


template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_iterator>
vertices(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<Graph >::vertex_iterator g_vertex_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    std::pair<g_vertex_iterator, g_vertex_iterator> original_vertices = vertices(w.graph());

    return make_range(vertex_iterator(predicate, original_vertices.first, original_vertices.second),
                      vertex_iterator(predicate, original_vertices.second, original_vertices.second));
}

template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_iterator>
edges(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::edge_iterator edge_iterator;
    typedef typename boost::graph_traits<Graph >::edge_iterator g_edge_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    std::pair<g_edge_iterator, g_edge_iterator> original_edges = edges(w.graph());

    return make_range(edge_iterator(predicate, original_edges.first, original_edges.second),
                      edge_iterator(predicate, original_edges.second, original_edges.second));
}

template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::out_edge_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
          const Connected_component_graph<Graph, FaceComponentMap> & w)
{

    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits<Graph >::out_edge_iterator g_out_edge_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    g_out_edge_iterator b,e;
    boost::tie(b,e) = out_edges(v, w.graph());
    return make_range(out_edge_iterator(predicate, b, e),
                      out_edge_iterator(predicate, e, e));
}

template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::in_edge_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::in_edge_iterator>
in_edges(typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{

    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::in_edge_iterator in_edge_iterator;
    typedef typename boost::graph_traits<Graph >::in_edge_iterator g_in_edge_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    g_in_edge_iterator b,e;
    boost::tie(b,e) = in_edges(v, w.graph());
    return make_range(in_edge_iterator(predicate, b, e),
                      in_edge_iterator(predicate, e, e));
}

//
// HalfedgeGraph
//
template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor
edge(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
     const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(CGAL::is_valid(h, w));
    return edge(h, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits<  Connected_component_graph<Graph, FaceComponentMap> >::edge_descriptor e,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(CGAL::is_valid(e, w));
    return halfedge(e, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(v, w));
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h = halfedge(v, w.graph());
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor hcirc = h;
    do
    {
        if(is_valid(hcirc, w))
            return hcirc;
        hcirc = opposite(next(hcirc, w.graph()), w.graph());
    }while(hcirc != h);
    return boost::graph_traits< CGAL::Connected_component_graph<Graph, FaceComponentMap> >::null_halfedge();
}


template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor, bool>
halfedge(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor u,
         typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor v,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(u, w) && is_valid(v, w));
    typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h = halfedge(u, v, w.graph()).first;
    return std::make_pair(h, is_valid(h, w));
}


template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
opposite(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(h, w) );
    return opposite(h, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor
source(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
       const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(h, w) );
    return source(h, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::vertex_descriptor
target(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
       const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(h, w) );
    return target(h, w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
next(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
     const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(is_valid(h, w));
    if(is_valid(face(h, w.graph()), w))
        return next(h, w.graph());

    //act as a border
    typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor hcirc = h;
    do
    {
        if(is_valid(hcirc, w))
        {
            return hcirc;
        }
        hcirc = opposite(next(hcirc,w.graph()),w.graph());
    }while(hcirc != h);
    return boost::graph_traits< CGAL::Connected_component_graph<Graph, FaceComponentMap> >::null_halfedge();
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
prev(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
     const Connected_component_graph<Graph, FaceComponentMap> & w)
{

    CGAL_assertion(is_valid(h, w));
    if(is_valid(face(h, w.graph()), w))
        return prev(h, w.graph());

    //act as a border
    typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor hcirc = h;
    do
    {
        if(is_valid(hcirc, w))
        {
            return hcirc;
        }
        hcirc = opposite(prev(hcirc,w.graph()), w.graph());
    }while(hcirc != h);
    return boost::graph_traits< CGAL::Connected_component_graph<Graph, FaceComponentMap> >::null_halfedge();
}

//
// HalfedgeListGraph
//

template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_iterator>
halfedges(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::halfedge_iterator halfedge_iterator;
    typedef typename boost::graph_traits<Graph >::halfedge_iterator g_halfedge_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    std::pair<g_halfedge_iterator, g_halfedge_iterator> original_halfedges = halfedges(w.graph());

    return make_range(halfedge_iterator(predicate, original_halfedges.first, original_halfedges.second),
                      halfedge_iterator(predicate, original_halfedges.second, original_halfedges.second));
}


template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Graph>::halfedges_size_type
num_halfedges(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    return num_halfedges(w.graph());
}

// FaceGraph
template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::face_descriptor
face(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor h,
     const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(CGAL::is_valid(h, w));
    return face(h,w.graph());
}

template <class Graph, typename FaceComponentMap >
typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Connected_component_graph<Graph, FaceComponentMap> >::face_descriptor f,
         const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    CGAL_assertion(CGAL::is_valid(f, w));
    return halfedge(f,w.graph());
}


template <class Graph, typename FaceComponentMap >
std::pair<typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::face_iterator,
typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::face_iterator>
faces(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    typedef typename boost::graph_traits<Connected_component_graph<Graph, FaceComponentMap> >::face_iterator face_iterator;
    typedef typename boost::graph_traits<Graph >::face_iterator g_face_iterator;

    typename Connected_component_graph<Graph, FaceComponentMap> ::Is_simplex_valid predicate(&w);
    std::pair<g_face_iterator, g_face_iterator> original_faces = faces(w.graph());

    return make_range(face_iterator(predicate, original_faces.first, original_faces.second),
                      face_iterator(predicate, original_faces.second, original_faces.second));
}



template <class Graph, typename FaceComponentMap >
typename boost::graph_traits<Graph>::vertices_size_type
num_faces(const Connected_component_graph<Graph, FaceComponentMap> & w)
{
    return num_faces(w.graph());
}


template <class Graph, typename FaceComponentMap >
bool
is_valid(const Connected_component_graph<Graph, FaceComponentMap> & w, bool verbose = false)
{
    return is_valid(w.graph(),verbose);
}

template <class Graph, typename FaceComponentMap, class PropertyTag>
typename boost::property_map<Graph, PropertyTag >::type
get(PropertyTag ptag, const Connected_component_graph<Graph,FaceComponentMap>& w)
{
  return get(ptag, w.graph());
}


template <class Graph, typename FaceComponentMap, class PropertyTag>
typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type
get(PropertyTag ptag,
    const Connected_component_graph<Graph, FaceComponentMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k)
{
  return get(ptag, w.graph(), k);
}


template <class Graph, typename FaceComponentMap, class PropertyTag>
void
put(PropertyTag ptag, const Connected_component_graph<Graph, FaceComponentMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k,
    typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type& v)
{
  put(ptag, w.graph(), k, v);
}

}//end namespace CGAL

namespace boost {
  template <typename Graph, typename FaceComponentMap, typename PropertyTag>
  struct property_map<CGAL::Connected_component_graph<Graph, FaceComponentMap>,PropertyTag> {
    typedef typename boost::property_map<Graph, PropertyTag >::type type;
    typedef typename boost::property_map<Graph, PropertyTag >::const_type const_type;
  };

  template<typename Graph, typename FaceComponentMap, typename PropertyTag>
  struct graph_has_property<CGAL::Connected_component_graph<Graph, FaceComponentMap>, PropertyTag>
    : graph_has_property<Graph, PropertyTag> {};

}// namespace boost

#endif // CGAL_BOOST_GRAPH_CONNECTED_COMPONENT_GRAPH_H
