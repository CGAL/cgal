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
// Author(s) : Andreas Fabri

#ifndef CGAL_GRAPH_GEOMETRY_H
#define CGAL_GRAPH_GEOMETRY_H

#include <cmath>
#include <limits>

#include <boost/graph/graph_traits.hpp>
#include <boost/range/distance.hpp>

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Origin.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_concepts.h>
#include <boost/concept/assert.hpp>

namespace CGAL {

/// \ingroup PkgBGLGraphGeometry
/// @{

template<typename HalfedgeGraph, 
         typename PositionMap, 
         typename NormalMap>
void calculate_face_normals(const HalfedgeGraph& g, 
                            PositionMap pm, 
                            NormalMap nm) 
{
  typedef boost::graph_traits<HalfedgeGraph> GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  typedef typename GraphTraits::enclosure_iterator enc_iterator;
  typedef typename boost::property_traits<PositionMap>::value_type position;
  typedef typename boost::property_traits<NormalMap>::value_type normal;

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb)
  {
    edge_descriptor edg = edge(*fb, g);
    edge_descriptor edgb = edg;
    enc_iterator eb, ee;
    boost::tie(eb, ee) = enclosure(*fb, g);

    position p0 = pm[target(edg, g)];
    edg = next(edg, g);
    position p1 = pm[target(edg, g)];
    edg = next(edg, g);
    position p2 = pm[target(edg, g)];
    edg = next(edg, g);
      
    if(edg == edgb) {
      // triangle 
      nm[*fb] = CGAL::unit_normal(p1, p2, p0);
    } else {
      normal n(CGAL::NULL_VECTOR);
      do {
        n = n + CGAL::normal(p1, p2, p0);
        p0 = p1;
        p1 = p2;

        edg = next(edg, g);
        p2 = pm[target(edg, g)];
      } while(edg != edgb);
        
      nm[*fb] = n / CGAL::sqrt(n.squared_length());
    }
  }
}

/// \cond SKIP_IN_MANUAL
namespace internal {
template<typename T>
typename Kernel_traits<T>::Kernel::FT
dot(const T& t1, const T& t2) 
{
  return t1[0]*t2[0]
    + t1[1]*t2[1]
    + t1[2]*t2[2];
}
  
} // internal
/// \endcond

template<typename HalfedgeGraph, 
         typename Position, 
         typename Normal,
         typename Boundary>
void calculate_vertex_normals(const HalfedgeGraph& g, 
                              Position position_map, 
                              Normal normal_map,
                              Boundary boundary_map)
{
  typedef boost::graph_traits<HalfedgeGraph> GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::out_edge_iterator out_edge_iterator;
    
  typedef typename boost::property_traits<Normal>::value_type normal;
  typedef typename boost::property_traits<Position>::value_type position;
  typedef typename Kernel_traits<normal>::Kernel::FT FT;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb)
  {
    normal nn(CGAL::NULL_VECTOR);
      
    position p0 = position_map[*vb];
    out_edge_iterator oeb, oee;
      
    for(boost::tie(oeb, oee) = out_edges(*vb, g); oeb != oee; ++oeb)
    {
      if(!boundary_map[target(*oeb, g)])
      {
        position p1 = position_map[target(*oeb, g)];
        p1 = p1 - (p0 - CGAL::ORIGIN);
        FT length = CGAL::sqrt((p1 - CGAL::ORIGIN).squared_length());
        if(length > (std::numeric_limits<FT>::min)())
          p1 = CGAL::ORIGIN + ((p1 - CGAL::ORIGIN) * 1.0 / length);
          
        position p2 = position_map[source(prev(halfedge(*oeb, g), g), g)];
        p2 = p2 - (p0 - CGAL::ORIGIN);
        length = CGAL::sqrt((p2 - CGAL::ORIGIN).squared_length());
        if(length > (std::numeric_limits<FT>::min)())
          p2 = CGAL::ORIGIN + ((p2 - CGAL::ORIGIN) * 1.0 / length);
          
        FT cosine
          = internal::dot(p1, p2) / CGAL::sqrt(internal::dot(p1, p1) * internal::dot(p2, p2));
        if(cosine < -1.0) cosine = -1.0;
        else if(cosine > 1.0) cosine = 1.0;

        FT angle = std::acos(cosine);

        normal n = CGAL::unit_normal(position(CGAL::ORIGIN), p1, p2);
        n = n * angle;
        nn = nn + n;
      }
    }
      
    FT length = CGAL::sqrt(nn.squared_length());
    if(length > (std::numeric_limits<FT>::min)())
      nn = nn * (1.0 / length);
    normal_map[*vb] = nn;
  }
}


/// Convenience function that calls `triangle()` with the property map
/// obtained by `get(CGAL::vertex_point, g)`.
template <typename HalfedgeGraph>
typename Kernel_traits<
  typename boost::property_traits< 
    typename boost::property_map< HalfedgeGraph, CGAL::vertex_point_t>::const_type 
    >::value_type
    >::Kernel::Triangle_3
triangle(const HalfedgeGraph& g,
         typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h) {
  return triangle(g, h, get(CGAL::vertex_point, g));
}


/// `triangle()` returns a `Triangle_3` constructed
/// from the positions of the `vertex_descriptors` around the face
/// incident to `h` in counter-clockwise direction.
///
/// The Kernel of the returned `Triangle_3` is the same as the
/// Kernel of the `value_type` of the.
///
/// \pre The face incident to h is triangular.
///
/// \tparam PositionMap must be a model of \ref ReadablePropertyMap. 
///        Its value_type must be a model of `Kernel::Point_3`
/// \tparam HalfedgeGraph  must be a model of \ref HalfedgeGraph.
/// \param g The graph 
/// \param h The halfedge
/// \param pm The map used to find the vertices of the triangle.
template <typename HalfedgeGraph,
          typename PositionMap>
typename Kernel_traits<typename boost::property_traits<PositionMap>::value_type>::Kernel::Triangle_3
triangle(const HalfedgeGraph& g,
         typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h,
         const PositionMap& pm)
{
  BOOST_CONCEPT_ASSERT((HalfedgeGraphConcept<HalfedgeGraph>));

  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor vertex_descriptor;
  
  typedef typename Kernel_traits<
    typename boost::property_traits<PositionMap>::value_type
    >::Kernel::Triangle_3 Triangle_3;
  
  CGAL_assertion_code(
    typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor h2 = h;
  )
  vertex_descriptor u = target(h,g);
  h = next(h,g);
  vertex_descriptor v = target(h,g);
  h = next(h,g);
  vertex_descriptor w = target(h,g);
  h = next(h,g);

  // are we really looking at a triangle?
  CGAL_assertion(h == h2);

  return Triangle_3(get(pm, u),
                    get(pm,v),
                    get(pm,w));
}




/// @}

} // CGAL

#endif /* CGAL_GRAPH_GEOMETRY_H */
