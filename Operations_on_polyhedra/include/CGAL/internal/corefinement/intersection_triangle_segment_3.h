// Copyright (c) 2011 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_H
#define CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_H

#include <CGAL/license/Polygon_mesh_processing.h>


//TODO rename this file when doing proper integration
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>
#include <CGAL/internal/Intersections_3/Triangle_3_Segment_3_intersection.h>
namespace CGAL{
namespace internal_IOP{

enum Intersection_type {FACET,EDGE,VERTEX,EMPTY,COPLNR};


//define shortcut for intersection result types
template <class Polyhedron,class Is_const>
struct Intersection_types{
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle Intersection_info;
  typedef cpp11::tuple<Intersection_type,Intersection_info,bool,bool> Intersection_result;
};


template<class Polyhedron,class Kernel,class Is_const>
typename Intersection_types<Polyhedron,Is_const>::Intersection_result
find_intersection(const typename Kernel::Point_3& p, const typename Kernel::Point_3& q, 
                  const typename Kernel::Point_3& a, const typename Kernel::Point_3& b, const typename Kernel::Point_3& c,
                  typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle /*hh*/,
                  typename Polyhedron_types<Polyhedron,Is_const>::Facet_handle  fh,
                  bool is_vertex_coplanar=false,bool is_vertex_opposite_coplanar=false)
{
  typedef typename Intersection_types<Polyhedron,Is_const>::Intersection_info Intersection_info;
  typedef typename Intersection_types<Polyhedron,Is_const>::Intersection_result Intersection_result;
  
  Orientation ab=orientation(p,q,a,b);
  Orientation bc=orientation(p,q,b,c);
  Orientation ca=orientation(p,q,c,a);
  
  if ( ab==POSITIVE || bc==POSITIVE || ca==POSITIVE )
    return Intersection_result(EMPTY,Intersection_info(),false,false);
  
  int nb_coplanar=(ab==COPLANAR?1:0) + (bc==COPLANAR?1:0) + (ca==COPLANAR?1:0);
  
  if ( nb_coplanar==0 )
    return cpp11::make_tuple(FACET,Intersection_info(fh->halfedge()),is_vertex_coplanar,is_vertex_opposite_coplanar);
  
  if (nb_coplanar==1){
    if (ab==COPLANAR)
      return cpp11::make_tuple(EDGE,Intersection_info(fh->halfedge()->next()),is_vertex_coplanar,is_vertex_opposite_coplanar);
    if (bc==COPLANAR)
      return cpp11::make_tuple(EDGE,Intersection_info(fh->halfedge()->next()->next()),is_vertex_coplanar,is_vertex_opposite_coplanar);
    CGAL_assertion(ca==COPLANAR);
    return cpp11::make_tuple(EDGE,Intersection_info(fh->halfedge()),is_vertex_coplanar,is_vertex_opposite_coplanar);
  }
  
  CGAL_assertion(nb_coplanar==2);
  
  if (ab!=COPLANAR)
    return cpp11::make_tuple(VERTEX,Intersection_info(fh->halfedge()->next()->next()),is_vertex_coplanar,is_vertex_opposite_coplanar);
  if (bc!=COPLANAR)
    return cpp11::make_tuple(VERTEX,Intersection_info(fh->halfedge()),is_vertex_coplanar,is_vertex_opposite_coplanar);
  CGAL_assertion(ca!=COPLANAR);
  return cpp11::make_tuple(VERTEX,Intersection_info(fh->halfedge()->next()),is_vertex_coplanar,is_vertex_opposite_coplanar);
}


template<class Polyhedron, class Kernel, class Is_const, class PolyhedronPointPmap>
typename Intersection_types<Polyhedron,Is_const>::Intersection_result
do_intersect(typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle hh,
             typename Polyhedron_types<Polyhedron,Is_const>::Facet_handle fh,
             PolyhedronPointPmap ppmap)
{
  typedef typename Intersection_types<Polyhedron,Is_const>::Intersection_info   Intersection_info;
  typedef typename Intersection_types<Polyhedron,Is_const>::Intersection_result Intersection_result;
  
  const typename Kernel::Point_3 & a = get(ppmap, fh->halfedge()->vertex() );
  const typename Kernel::Point_3 & b = get(ppmap, fh->halfedge()->next()->vertex() );
  const typename Kernel::Point_3 & c = get(ppmap, fh->halfedge()->next()->next()->vertex() );
  const typename Kernel::Point_3 & p = get(ppmap, hh->vertex() );
  const typename Kernel::Point_3 & q = get(ppmap, hh->opposite()->vertex() );


  const Orientation abcp = orientation(a,b,c,p);
  const Orientation abcq = orientation(a,b,c,q);


  switch ( abcp ) {
  case POSITIVE:
    switch ( abcq ) {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
      return Intersection_result(EMPTY,Intersection_info(),false,false);
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      return find_intersection<Polyhedron,Kernel,Is_const>(p,q,a,b,c,hh,fh);
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
    return find_intersection<Polyhedron,Kernel,Is_const>(p,q,a,b,c,hh,fh,false,true);
      
    default: // should not happen.
      CGAL_assertion(false);
      return Intersection_result(EMPTY,Intersection_info(),false,false);
    }
  case NEGATIVE:
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection<Polyhedron,Kernel,Is_const>(q,p,a,b,c,hh,fh);
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
      return Intersection_result(EMPTY,Intersection_info(),false,false);
    case COPLANAR:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      return find_intersection<Polyhedron,Kernel,Is_const>(q,p,a,b,c,hh,fh,false,true);
    default: // should not happen.
      CGAL_assertion(false);
      return Intersection_result(EMPTY,Intersection_info(),false,false);
    }
  case COPLANAR: // p belongs to the triangle's supporting plane
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection<Polyhedron,Kernel,Is_const>(q,p,a,b,c,hh,fh,true,false);
    case NEGATIVE:
      // q sees the triangle in clockwise order
      return find_intersection<Polyhedron,Kernel,Is_const>(p,q,a,b,c,hh,fh,true,false);    
    case COPLANAR:
      // the segment is coplanar with the triangle's supporting plane
      // we test whether the segment intersects the triangle in the common 
      // supporting plane
      if ( ::CGAL::internal::do_intersect_coplanar(a,b,c,p,q,Kernel()) )
        return Intersection_result(COPLNR,Intersection_info(),true,true);
      return Intersection_result(EMPTY,Intersection_info(),true,true);
      
    default: // should not happen.
      CGAL_assertion(false);
      return Intersection_result(EMPTY,Intersection_info(),false,false);
    }
  default: // should not happen.
    CGAL_assertion(false);
    return Intersection_result(EMPTY,Intersection_info(),false,false);
  }
}

}} //namespace CGAL::internal_IOP

#endif //CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_H
