// Copyright (c) 2012 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H
#define CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H

#include <CGAL/property_map.h>
#include <boost/type_traits/is_const.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <boost/type_traits/remove_const.hpp>

namespace CGAL{

//property map
template <class Polyhedron,
           class VertexPointPMap >
struct Triangle_from_facet_handle_property_map{  
  typename boost::remove_const<Polyhedron>::type* g;
  VertexPointPMap m_vppm;

  Triangle_from_facet_handle_property_map()
  {}

  Triangle_from_facet_handle_property_map(Polyhedron* p)
    : g( const_cast<typename boost::remove_const<Polyhedron>::type*>(p) ),
      m_vppm( get(vertex_point, *g) )
  {}

  Triangle_from_facet_handle_property_map(Polyhedron* g,
                                          VertexPointPMap vppm )
    : g(const_cast<typename boost::remove_const<Polyhedron>::type*>(g)),
      m_vppm(vppm)
  {}

  typedef typename boost::property_traits< VertexPointPMap >::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;
  //classical typedefs
 
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor key_type;
  typedef Triangle_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  //get function for property map
  inline friend
  Triangle_3
  get(const Triangle_from_facet_handle_property_map<Polyhedron,VertexPointPMap>& m,
      typename Triangle_from_facet_handle_property_map<Polyhedron,VertexPointPMap>::key_type f)
  {
    typedef typename boost::property_traits< VertexPointPMap >::value_type Point_3;
    typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;

    CGAL_precondition(halfedge(f,*m.g) == next(next(next(halfedge(f,*m.g),*m.g),*m.g),*m.g));
    const Point_3& a = get(m.m_vppm, target(halfedge(f,*m.g),*m.g));
    const Point_3& b = get(m.m_vppm, target(next(halfedge(f,*m.g),*m.g),*m.g));
    const Point_3& c = get(m.m_vppm,target(next(next(halfedge(f,*m.g),*m.g),*m.g),*m.g));
 
    return Triangle_3(a,b,c);
  }
};


template < class HalfedgeGraph,
           class VertexPointPMap >
struct Segment_from_edge_descriptor_property_map{

  Segment_from_edge_descriptor_property_map()
  {}

Segment_from_edge_descriptor_property_map(HalfedgeGraph* g)
  : m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
    m_vppm( get(vertex_point, *m_graph) )
  {}

  Segment_from_edge_descriptor_property_map(
    HalfedgeGraph* g,
    VertexPointPMap vppm ) :
  m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
  m_vppm(vppm)
  {}

  //classical typedefs
  typedef typename boost::property_traits< VertexPointPMap >::value_type Point;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor key_type;
  typedef typename Kernel_traits<Point>::Kernel::Segment_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  //data
  typename boost::remove_const<HalfedgeGraph>::type* m_graph;
  VertexPointPMap m_vppm;

  //get function for property map
  inline friend
  value_type
  get(Segment_from_edge_descriptor_property_map<HalfedgeGraph,VertexPointPMap> pmap,
      key_type h)
  {
    typedef typename boost::property_map< HalfedgeGraph, vertex_point_t>::type Point_pmap;
    typedef typename boost::property_traits< Point_pmap >::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel::Segment_3 Segment_3;

    return Segment_3(get(pmap.m_vppm, boost::source(h, *pmap.m_graph) ),
                     get(pmap.m_vppm, boost::target(h, *pmap.m_graph) ) );
  }
};

//property map to access a point from a facet handle
template <class Polyhedron,
          class VertexPointPMap>
struct One_point_from_facet_handle_property_map{

  One_point_from_facet_handle_property_map(Polyhedron* g = NULL)
    : m_graph( const_cast<typename boost::remove_const<Polyhedron>::type*>(g) )
  {}

template <class Polyhedron,
          class VertexPointPMap>
  One_point_from_facet_handle_property_map(Polyhedron* gm, VertexPointPMap vppm )
    : m_graph( const_cast<typename boost::remove_const<Polyhedron>::type*>(g) ),
      m_vppm(vppm)
  {}

  typename boost::remove_const<Polyhedron>::type* m_graph;
  VertexPointPMap m_vppm;

  //classical typedefs
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor key_type;
  typedef typename boost::property_traits< VertexPointPMap >::value_type value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  //get function for property map
  inline friend
  reference
  get(const One_point_from_facet_handle_property_map<Polyhedron,VertexPointPMap>& m,
      key_type f)
  {
    return get(m.m_vppm, target(halfedge(f, *m.m_graph), *m.m_graph));
  }
};

//property map to access a point from an edge
template < class HalfedgeGraph,
           class VertexPointPMap >
struct Source_point_from_edge_descriptor{
  Source_point_from_edge_descriptor(
    HalfedgeGraph* g = NULL ) :
  m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
  m_vppm( get(vertex_point, *m_graph) )
  {}

  Source_point_from_edge_descriptor(
    HalfedgeGraph* g,
    VertexPointPMap vppm ) :
  m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
  m_vppm(vppm)
  {}

  //classical typedefs
  typedef typename boost::property_traits< VertexPointPMap >::value_type value_type;
  typedef typename boost::property_traits< VertexPointPMap >::reference reference;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor key_type;
  typedef boost::readable_property_map_tag category;
  //data
  typename boost::remove_const<HalfedgeGraph>::type* m_graph;
  VertexPointPMap m_vppm;

  //get function for property map
  inline friend
  reference
  get(Source_point_from_edge_descriptor<HalfedgeGraph,VertexPointPMap> pmap,
      key_type h)
  {
    return  get(vertex_point,
                       *pmap.m_graph,
                       boost::source(h, *pmap.m_graph) );
  }
};

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H
