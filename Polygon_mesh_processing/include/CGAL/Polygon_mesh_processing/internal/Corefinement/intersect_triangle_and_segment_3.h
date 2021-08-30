// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_PMP_INTERSECT_TRIANGLE_AND_SEGMENT_3_H
#define CGAL_INTERNAL_PMP_INTERSECT_TRIANGLE_AND_SEGMENT_3_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/Intersections_3/Segment_3_Triangle_3.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Intersection_type.h>
#include <CGAL/property_map.h>

namespace CGAL{
namespace Polygon_mesh_processing {
namespace Corefinement{


template<class TriangleMesh,class Point_3>
std::tuple<Intersection_type,
             typename boost::graph_traits<TriangleMesh>::halfedge_descriptor,
             bool,bool>
find_intersection(const Point_3& p, const Point_3& q,  //segment
                  const Point_3& a, const Point_3& b, const Point_3& c, //triangle
                  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd, // halfedge of the triangle face, its target is a
                  const TriangleMesh& tm,
                  bool is_src_coplanar=false,bool is_tgt_coplanar=false)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef std::tuple<Intersection_type,halfedge_descriptor,bool,bool> result_type;

  Orientation ab=orientation(p,q,a,b);
  Orientation bc=orientation(p,q,b,c);
  Orientation ca=orientation(p,q,c,a);

  if ( ab==POSITIVE || bc==POSITIVE || ca==POSITIVE )
    return result_type(EMPTY,GT::null_halfedge(),false,false);

  int nb_coplanar=(ab==COPLANAR?1:0) + (bc==COPLANAR?1:0) + (ca==COPLANAR?1:0);

  if ( nb_coplanar==0 )
    return result_type(ON_FACE,hd,is_src_coplanar,is_tgt_coplanar);

  if (nb_coplanar==1){
    if (ab==COPLANAR)
      // intersection is ab
      return result_type(ON_EDGE,next(hd,tm),is_src_coplanar,is_tgt_coplanar);
    if (bc==COPLANAR)
      // intersection is bc
      return result_type(ON_EDGE,prev(hd,tm),is_src_coplanar,is_tgt_coplanar);
    CGAL_assertion(ca==COPLANAR);
    // intersection is ca
    return result_type(ON_EDGE,hd,is_src_coplanar,is_tgt_coplanar);
  }

  CGAL_assertion(nb_coplanar==2);

  if (ab!=COPLANAR)
    // intersection is c
    return result_type(ON_VERTEX,prev(hd,tm),is_src_coplanar,is_tgt_coplanar);
  if (bc!=COPLANAR)
    // intersection is a
    return result_type(ON_VERTEX,hd,is_src_coplanar,is_tgt_coplanar);
  CGAL_assertion(ca!=COPLANAR);
  // intersection is b
  return result_type(ON_VERTEX,next(hd,tm),is_src_coplanar,is_tgt_coplanar);
}


template<class TriangleMesh, class VertexPointMap1, class VertexPointMap2>
std::tuple<Intersection_type,
             typename boost::graph_traits<TriangleMesh>::halfedge_descriptor,
             bool,bool>
intersection_type(
             typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h_1,
             typename boost::graph_traits<TriangleMesh>::face_descriptor f_2,
             const TriangleMesh& tm1,
             const TriangleMesh& tm2,
             const VertexPointMap1& vpm1,
             const VertexPointMap2& vpm2)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef std::tuple<Intersection_type,halfedge_descriptor,bool,bool> result_type;
  typedef typename boost::property_traits<VertexPointMap1>::reference Point_ref1;
  typedef typename boost::property_traits<VertexPointMap2>::reference Point_ref2;
  typedef typename boost::property_traits<VertexPointMap1>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel Kernel;

  CGAL_static_assertion((std::is_same<Point_3, typename boost::property_traits<VertexPointMap2>::value_type>::value));

  halfedge_descriptor h_2=halfedge(f_2,tm2);

  Point_ref2 a = get(vpm2, target(h_2,tm2) );
  Point_ref2 b = get(vpm2, target(next(h_2,tm2),tm2) );
  Point_ref2 c = get(vpm2, source(h_2,tm2) );
  Point_ref1 p = get(vpm1, source(h_1,tm1) );
  Point_ref1 q = get(vpm1, target(h_1,tm1) );

  const Orientation abcp = orientation(a,b,c,p);
  const Orientation abcq = orientation(a,b,c,q);

  switch ( abcp ) {
  case POSITIVE:
    switch ( abcq ) {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
      return result_type(EMPTY,GT::null_halfedge(),false,false);
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      return find_intersection(p,q,a,b,c,h_2,tm2);
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
    return find_intersection(p,q,a,b,c,h_2,tm2,false,true);

    default: // should not happen.
      CGAL_assertion(false);
      return result_type(EMPTY,GT::null_halfedge(),false,false);
    }
  case NEGATIVE:
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection(q,p,a,b,c,h_2,tm2);
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
      return result_type(EMPTY,GT::null_halfedge(),false,false);
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      return find_intersection(q,p,a,b,c,h_2,tm2,false,true);
    default: // should not happen.
      CGAL_assertion(false);
      return result_type(EMPTY,GT::null_halfedge(),false,false);
    }
  case COPLANAR: // p belongs to the triangle's supporting plane
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection(q,p,a,b,c,h_2,tm2,true);
    case NEGATIVE:
      // q sees the triangle in clockwise order
      return find_intersection(p,q,a,b,c,h_2,tm2,true);
    case COPLANAR:
      // the segment is coplanar with the triangle's supporting plane
      // we test whether the segment intersects the triangle in the common
      // supporting plane
      if ( ::CGAL::Intersections::internal::do_intersect_coplanar(a,b,c,p,q,Kernel()) )
        return result_type(COPLANAR_TRIANGLES,GT::null_halfedge(),true,true);
      return result_type(EMPTY,GT::null_halfedge(),true,true);

    default: // should not happen.
      CGAL_assertion(false);
      return result_type(EMPTY,GT::null_halfedge(),false,false);
    }
  default: // should not happen.
    CGAL_assertion(false);
    return result_type(EMPTY,GT::null_halfedge(),false,false);
  }
}

} } } // CGAL::Polygon_mesh_processing::Corefinement

#endif //CGAL_INTERNAL_PMP_INTERSECT_TRIANGLE_AND_SEGMENT_3_H
