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
template <class Polyhedron>
struct Triangle_from_facet_handle_property_map{
  Triangle_from_facet_handle_property_map(Polyhedron* = NULL){}
  typedef typename Kernel_traits<
    typename Polyhedron::Vertex::Point>::Kernel::Triangle_3 Triangle_3;
  //classical typedefs
  typedef typename boost::mpl::if_<
    typename boost::is_const<Polyhedron>::type,
    typename Polyhedron::Facet_const_handle,
    typename Polyhedron::Facet_handle >::type key_type;
  typedef Triangle_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  //get function for property map
  inline friend
  Triangle_3
  get(Triangle_from_facet_handle_property_map<Polyhedron>,
      typename Triangle_from_facet_handle_property_map<Polyhedron>::key_type f)
  {
    typedef typename Polyhedron::Traits Kernel;
    CGAL_precondition(f->halfedge() == f->halfedge()->next()->next()->next());
    const typename Kernel::Point_3& a = f->halfedge()->vertex()->point();
    const typename Kernel::Point_3& b = f->halfedge()->next()->vertex()->point();
    const typename Kernel::Point_3& c = f->halfedge()->next()->next()->vertex()->point();
    return typename Kernel::Triangle_3(a,b,c);
  }
};


template < class HalfedgeGraph,
           class VertexPointPMap >
struct Segment_from_edge_descriptor_property_map{
  Segment_from_edge_descriptor_property_map(
    HalfedgeGraph* g = NULL ) :
  m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
  m_vppm( boost::get(vertex_point, *m_graph) )
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
    typedef typename Kernel_traits<Point>::Kernel Kernel;

    return typename Kernel::Segment_3(
      get(pmap.m_vppm, boost::source(h, *pmap.m_graph) ),
      get(pmap.m_vppm, boost::target(h, *pmap.m_graph) ) );
  }
};

//property map to access a point from a facet handle
template <class Polyhedron>
struct One_point_from_facet_handle_property_map{
  One_point_from_facet_handle_property_map(Polyhedron* = NULL){}
  //classical typedefs
  typedef typename boost::mpl::if_<
    typename boost::is_const<Polyhedron>::type,
    typename Polyhedron::Facet_const_handle,
    typename Polyhedron::Facet_handle >::type key_type;
  typedef typename Polyhedron::Vertex::Point_3 value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  //get function for property map
  inline friend
  reference
  get(One_point_from_facet_handle_property_map<Polyhedron>,
      key_type f)
  {
    return f->halfedge()->vertex()->point();
  }
};

//property map to access a point from an edge
template < class HalfedgeGraph,
           class VertexPointPMap >
struct Source_point_from_edge_descriptor{
  Source_point_from_edge_descriptor(
    HalfedgeGraph* g = NULL ) :
  m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
  m_vppm( boost::get(vertex_point, *m_graph) )
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
    return  boost::get(vertex_point,
                       *pmap.m_graph,
                       boost::source(h, *pmap.m_graph) );
  }
};

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H
