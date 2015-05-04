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

#ifndef HALFEDGE_AND_FACE_GRAPH_PROPERTY_MAPS_H
#define HALFEDGE_AND_FACE_GRAPH_PROPERTY_MAPS_H

#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/type_traits/is_const.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace CGAL{

//property map
template <class FaceGraph,
           class VertexPointPMap >
struct Triangle_from_face_descriptor_property_map{
  typename boost::remove_const<FaceGraph>::type* m_graph;
  VertexPointPMap m_vppm;

  Triangle_from_face_descriptor_property_map() : m_graph(NULL)
  {}

  Triangle_from_face_descriptor_property_map(FaceGraph* g)
    : m_graph( const_cast<typename boost::remove_const<FaceGraph>::type*>(g) ),
      m_vppm( get(vertex_point, *m_graph) )
  {}

  Triangle_from_face_descriptor_property_map(FaceGraph* g,
                                          VertexPointPMap vppm )
    : m_graph(const_cast<typename boost::remove_const<FaceGraph>::type*>(g)),
      m_vppm(vppm)
  {}

  typedef typename boost::property_traits< VertexPointPMap >::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;
 
  //classical typedefs
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor key_type;
  typedef Triangle_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  //get function for property map
  inline friend
  Triangle_3
  get(const Triangle_from_face_descriptor_property_map<FaceGraph,VertexPointPMap>& pmap,
      typename Triangle_from_face_descriptor_property_map<FaceGraph,VertexPointPMap>::key_type f)
  {
    typedef typename boost::property_traits< VertexPointPMap >::value_type Point_3;
    typedef typename Kernel_traits<Point_3>::Kernel::Triangle_3 Triangle_3;
    typename boost::remove_const<FaceGraph>::type & g = *(pmap.m_graph);

    CGAL_precondition(halfedge(f,g) == next(next(next(halfedge(f,g),g),g),g));
    const Point_3& a = get(pmap.m_vppm, target(halfedge(f,g),g));
    const Point_3& b = get(pmap.m_vppm, target(next(halfedge(f,g),g),g));
    const Point_3& c = get(pmap.m_vppm,target(next(next(halfedge(f,g),g),g),g));
 
    return Triangle_3(a,b,c);
  }
};


template < class HalfedgeGraph,
           class VertexPointPMap >
struct Segment_from_edge_descriptor_property_map{

  Segment_from_edge_descriptor_property_map()  : m_graph(NULL)
  {}

Segment_from_edge_descriptor_property_map(HalfedgeGraph* g)
  : m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
    m_vppm( get(vertex_point, *m_graph) )
  {}

  Segment_from_edge_descriptor_property_map(HalfedgeGraph* g,
                                            VertexPointPMap vppm )
    : m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) ),
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
    typedef typename boost::property_traits< VertexPointPMap >::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel::Segment_3 Segment_3;

    return Segment_3(get(pmap.m_vppm, source(h, *pmap.m_graph) ),
                     get(pmap.m_vppm, target(h, *pmap.m_graph) ) );
  }
};

//property map to access a point from a facet handle
template <class FaceGraph,
          class VertexPointPMap>
struct One_point_from_face_descriptor_property_map{
  One_point_from_face_descriptor_property_map()  : m_graph(NULL)
  {}

  One_point_from_face_descriptor_property_map(FaceGraph* g)
    : m_graph( const_cast<typename boost::remove_const<FaceGraph>::type*>(g) )
    , m_vppm( get(vertex_point, *m_graph) )
  {}

  One_point_from_face_descriptor_property_map(FaceGraph* g, VertexPointPMap vppm )
    : m_graph( const_cast<typename boost::remove_const<FaceGraph>::type*>(g) ),
      m_vppm(vppm)
  {}

  typename boost::remove_const<FaceGraph>::type* m_graph;
  VertexPointPMap m_vppm;

  //classical typedefs
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor key_type;
  typedef typename boost::property_traits< VertexPointPMap >::value_type value_type;
  typedef typename boost::property_traits< VertexPointPMap >::reference reference;
  typedef boost::lvalue_property_map_tag category;

  //get function for property map
  inline friend
  reference
  get(const One_point_from_face_descriptor_property_map<FaceGraph,VertexPointPMap>& m,
      key_type f)
  {
    return get(m.m_vppm, target(halfedge(f, *m.m_graph), *m.m_graph));
  }
};

//property map to access a point from an edge
template < class HalfedgeGraph,
           class VertexPointPMap >
struct Source_point_from_edge_descriptor{
  Source_point_from_edge_descriptor()  : m_graph(NULL)
  {}

  Source_point_from_edge_descriptor(HalfedgeGraph* g)
    : m_graph( const_cast<typename boost::remove_const<HalfedgeGraph>::type*>(g) )
    , m_vppm( get(vertex_point, *m_graph) )
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
                       source(h, *pmap.m_graph) );
  }
};

} //namespace CGAL

#endif //HALFEDGE_AND_FACE_GRAPH_PROPERTY_MAPS_H
