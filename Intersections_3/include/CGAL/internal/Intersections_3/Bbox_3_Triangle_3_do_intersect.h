// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Triangle_3.h>
#include <CGAL/Bbox_3.h>

// Fast Triangle-Cuboid intersection test, following Tomas Akenine-Moeller description.
// The code looks slightly different from his code because we avoid the translation at
// a minimal cost (and we use C++ ;).

#include <CGAL/number_utils.h>
// ST: is this include really needed ?
#include <CGAL/internal/Intersections_3/Bbox_3_Plane_3_do_intersect.h>

namespace CGAL {

namespace internal {

  template <class K>
  inline
  bool do_bbox_intersect(const typename K::Triangle_3& triangle,
                         const CGAL::Bbox_3& bbox)
  {
    const typename K::Point_3& p = triangle.vertex(0);
    const typename K::Point_3& q = triangle.vertex(1);
    const typename K::Point_3& r = triangle.vertex(2);

    for(int i = 0; i < 3; ++i) {
      if(p[i] <= q[i]) {
        if(q[i] <= r[i]) { // pqr
          if(((bbox.max)(i) < p[i]) || ((bbox.min)(i) > r[i]))
            return false;
        }
        else {
          if(p[i] <= r[i]) { // prq
            if((bbox.max)(i) < p[i] || (bbox.min)(i) > q[i])
              return false;
          }
          else { // rpq
            if((bbox.max)(i) < r[i] || (bbox.min)(i) > q[i])
              return false;
          }
        }
      }
      else {
        if(p[i] <= r[i]) { // qpr
          if((bbox.max)(i) < q[i] || (bbox.min)(i) > r[i])
            return false;
        }
        else {
          if(q[i] <= r[i]) { // qrp
            if((bbox.max)(i) < q[i] || (bbox.min)(i) > p[i])
              return false;
          }
          else { // rqp
            if((bbox.max)(i) < r[i] || (bbox.min)(i) > p[i])
              return false;
          }
        }
      }
    }
    return true;
  }

  // AXE is the axe such that p is orthogonal to it.
  // if you do not know it, or if it does not exist,
  // use get_min_max without the AXE template parameter
  // available in _plane_is_cuboid_do_intersect.h
  template <class K, int AXE>
  inline
    void get_min_max(const typename K::FT& px,
    const typename K::FT& py,
    const typename K::FT& pz,
    const CGAL::Bbox_3& c,
    typename K::Point_3& p_min,
    typename K::Point_3& p_max)
  {
    if(AXE == 0 || px > 0) {
      if(AXE == 1 || py > 0) {
        if(AXE == 2 || pz > 0) {
          p_min = typename K::Point_3(c.xmin(), c.ymin(),c.zmin());
          p_max = typename K::Point_3(c.xmax(), c.ymax(),c.zmax());
        }
        else {
          p_min = typename K::Point_3(c.xmin(), c.ymin(),c.zmax());
          p_max = typename K::Point_3(c.xmax(), c.ymax(),c.zmin());
        }
      }
      else {
        if(AXE == 2 || pz > 0) {
          p_min = typename K::Point_3(c.xmin(), c.ymax(),c.zmin());
          p_max = typename K::Point_3(c.xmax(), c.ymin(),c.zmax());
        }
        else {
          p_min = typename K::Point_3(c.xmin(), c.ymax(),c.zmax());
          p_max = typename K::Point_3(c.xmax(), c.ymin(),c.zmin());
        }
      }
    }
    else {
      if(AXE == 1 || py > 0) {
        if(AXE == 2 || pz > 0) {
          p_min = typename K::Point_3(c.xmax(), c.ymin(),c.zmin());
          p_max = typename K::Point_3(c.xmin(), c.ymax(),c.zmax());
        }
        else {
          p_min = typename K::Point_3(c.xmax(), c.ymin(),c.zmax());
          p_max = typename K::Point_3(c.xmin(), c.ymax(),c.zmin());
        }
      }
      else {
        if(AXE == 2 || pz > 0) {
          p_min = typename K::Point_3(c.xmax(), c.ymax(),c.zmin());
          p_max = typename K::Point_3(c.xmin(), c.ymin(),c.zmax());
        }
        else {
          p_min = typename K::Point_3(c.xmax(), c.ymax(),c.zmax());
          p_max = typename K::Point_3(c.xmin(), c.ymin(),c.zmin());
        }
      }
    }
  }


  template <class K, int AXE, int SIDE>
  inline
  typename K::FT
  do_axis_intersect_aux(const typename K::FT& alpha,
                        const typename K::FT& beta,
                        const typename K::Vector_3* sides)
  {
    switch ( AXE )
    {
      case 0:
        return -sides[SIDE].z()*alpha + sides[SIDE].y()*beta;
      case 1:
        return sides[SIDE].z()*alpha - sides[SIDE].x()*beta;
      case 2:
        return -sides[SIDE].y()*alpha + sides[SIDE].x()*beta;
      default:
        CGAL_error();
        return typename K::FT(0.);
    }
  }

  //given a vector checks whether it is collinear to a base vector
  //of the orthonormal frame. return -1 otherwise
  template <class K>
  inline
  int
  collinear_axis(typename K::Vector_3 side){
    if ( certainly(side[0]==0) ){
      if ( certainly(side[1]==0) )        return 2;
      if ( certainly(side[2]==0) )        return 1;
    }
    else{
        if ( certainly(side[1]==0) && 
            certainly(side[2]==0) )       return 0;
    }
    return -1;
  }

  template <class K, int AXE, int SIDE>
  inline
  Uncertain<bool> do_axis_intersect(const typename K::Triangle_3& triangle,
                                    const typename K::Vector_3* sides,
                                    const CGAL::Bbox_3& bbox)
  {
    const typename K::Point_3* j = & triangle.vertex(SIDE);
    const typename K::Point_3* k = & triangle.vertex((SIDE+2)%3);

    typename K::Point_3 p_min, p_max;
    get_min_max<K, AXE>(AXE==0? 0: AXE==1? sides[SIDE].z(): -sides[SIDE].y(),
                        AXE==0? -sides[SIDE].z(): AXE==1? 0: sides[SIDE].x(),
                        AXE==0? sides[SIDE].y(): AXE==1? -sides[SIDE].x(): 0,
                        bbox, p_min, p_max);

    switch ( AXE )
    {
    case 0: {
      // t_max >= t_min
      Uncertain<bool> b =  do_axis_intersect_aux<K,AXE,SIDE>(k->y()-j->y(), k->z()-j->z(), sides) >= 0;
      if (is_indeterminate(b))
        return b;
      if(b) std::swap(j,k);
      return CGAL_OR( (do_axis_intersect_aux<K,AXE,SIDE>(p_min.y()-j->y(), p_min.z()-j->z(), sides) <= 0),
                      (do_axis_intersect_aux<K,AXE,SIDE>(p_max.y()-k->y(), p_max.z()-k->z(), sides) >= 0) );      
    }
    case 1: {
      // t_max >= t_min
      Uncertain<bool> b =  do_axis_intersect_aux<K,AXE,SIDE>(k->x()-j->x(), k->z()-j->z(), sides) >= 0;
      if (is_indeterminate(b))
        return b;
      if(b) std::swap(j,k);
      return  CGAL_OR( (do_axis_intersect_aux<K,AXE,SIDE>(p_min.x()-j->x(), p_min.z()-j->z(), sides) <= 0),
                       (do_axis_intersect_aux<K,AXE,SIDE>(p_max.x()-k->x(), p_max.z()-k->z(), sides) >= 0) );

    }
    case 2: {
      // t_max >= t_min
      Uncertain<bool> b = do_axis_intersect_aux<K,AXE,SIDE>(k->x()-j->x(), k->y()-j->y(), sides) >= 0;
      if ( is_indeterminate(b))
        return b;
      if(b) std::swap(j,k);
      return  CGAL_OR( (do_axis_intersect_aux<K,AXE,SIDE>(p_min.x()-j->x(), p_min.y()-j->y(), sides) <= 0),
                       (do_axis_intersect_aux<K,AXE,SIDE>(p_max.x()-k->x(), p_max.y()-k->y(), sides) >= 0) );
    }
    default:
      // Should not happen
      CGAL_error();
      return make_uncertain(false);
    }
  }

  template <class K>
  bool do_intersect(const typename K::Triangle_3& a_triangle,
    const CGAL::Bbox_3& a_bbox,
    const K& k)
  {

    if(certainly_not( do_bbox_intersect<K>(a_triangle, a_bbox) ))
      return false;

    if(certainly_not( do_intersect(a_triangle.supporting_plane(), a_bbox, k) ))
      return false;

    typename K::Vector_3 sides[3];
    sides[0] = a_triangle[1] - a_triangle[0];
    sides[1] = a_triangle[2] - a_triangle[1];
    sides[2] = a_triangle[0] - a_triangle[2];    
    int forbidden_axis=-1;
    int forbidden_size=-1;
    //determine whether one vector is collinear with an axis
    int tmp=collinear_axis<K>(sides[0]);
    if ( tmp!= -1){
      forbidden_axis=tmp;
      forbidden_size=0;
    }
    else{
      tmp=collinear_axis<K>(sides[1]);
      if ( tmp!= -1){
        forbidden_axis=tmp;
        forbidden_size=1;
      }
      else{
        tmp=collinear_axis<K>(sides[2]);
        if ( tmp!= -1){
          forbidden_axis=tmp;
          forbidden_size=2;
        }        
      }
    }
    
#if 0
    typename K::Point_3 p(a_bbox.xmin(), a_bbox.ymin(), a_bbox.zmin());
    typename K::Point_3 q(a_bbox.xmax(), a_bbox.ymax(), a_bbox.zmax());
    
    typename K::Point_3 m = CGAL::midpoint(p,q);
    typename K::Vector_3 v = m - CGAL::ORIGIN;
    
    typename K::Triangle_3 triangle(a_triangle[0]-v, a_triangle[1]-v, a_triangle[2]-v);
  
    Bbox_3 bbox( (p-v).bbox() + (q-v).bbox());

#else 
    const typename K::Triangle_3&  triangle = a_triangle;
    const Bbox_3& bbox = a_bbox;
#endif    
    
    Uncertain<bool> ind = make_uncertain(true);
    
    Uncertain<bool> b = make_uncertain(true);
    if (forbidden_axis!=0){
      if (forbidden_size!=0){
        b = do_axis_intersect<K,0,0>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=1){
        b = do_axis_intersect<K,0,1>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=2){
        b = do_axis_intersect<K,0,2>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
    }
    
    if (forbidden_axis!=1){
      if (forbidden_size!=0){
        b = do_axis_intersect<K,1,0>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=1){
        b = do_axis_intersect<K,1,1>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=2){
        b = do_axis_intersect<K,1,2>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
    }
    
    if (forbidden_axis!=2){
      if (forbidden_size!=0){
        b = do_axis_intersect<K,2,0>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=1){
        b = do_axis_intersect<K,2,1>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
      if (forbidden_size!=2){
        b = do_axis_intersect<K,2,2>(triangle, sides, bbox);
        if(is_indeterminate(b)){
          ind = b;
        } else if(! b){
          return false;
        }
      }
    }
    return ind; // throws exception in case it is indeterminate
  }

  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bbox,
                    const typename K::Triangle_3& triangle,
                    const K& k)
  {
    return do_intersect(triangle, bbox, k);
  }

} // namespace internal

template <class K>
bool do_intersect(const CGAL::Triangle_3<K>& triangle,
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(triangle, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox,
		  const CGAL::Triangle_3<K>& triangle)
{
  return typename K::Do_intersect_3()(triangle, bbox);
}

} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H
