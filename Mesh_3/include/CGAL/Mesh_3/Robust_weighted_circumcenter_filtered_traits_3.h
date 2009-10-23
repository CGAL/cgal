// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
#define CGAL_MESH_3_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H


#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Robust_construction.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>


namespace CGAL {


template < typename K_ >
class Robust_filtered_construct_weighted_circumcenter_3
{
public:
  typedef Exact_predicates_exact_constructions_kernel          EK;
  typedef Weighted_converter_3<Cartesian_converter<K_, EK> >   To_exact;
  typedef Weighted_converter_3<Cartesian_converter<EK,K_> >    Back_from_exact;
  
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K_> Rt;
  typedef CGAL::Regular_triangulation_euclidean_traits_3<EK> Exact_Rt;
  
  typedef typename Rt::Weighted_point_3               Weighted_point_3;
  typedef typename Rt::Bare_point                     Bare_point;
  typedef typename Rt::FT                             FT;
  typedef typename Rt::Sphere_3                       Sphere_3;
  
  typedef Bare_point                                  result_type;
  
  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r,
                          const Weighted_point_3 & s ) const
  {
    CGAL_precondition(! Rt().coplanar_3_object()(p,q,r,s) );
    
    typename Rt::Construct_weighted_circumcenter_3 weighted_circumcenter =
      Rt().construct_weighted_circumcenter_3_object();
    typename Rt::Has_on_bounded_side_3 on_bounded_side =
      Rt().has_on_bounded_side_3_object();
    
    // Compute denominator to swith to exact if it is 0
    const FT denom = compute_denom(p,q,r,s);
    if ( ! CGAL_NTS is_zero(denom) )
    {
      result_type point = weighted_circumcenter(p,q,r,s);
      
      // Fast output
      if ( on_bounded_side(Sphere_3(p,q,r,s),point) )
        return point;
    }
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();
    
    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r),
                                                       to_exact(s)));
  }
  
  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r ) const
  {
    CGAL_precondition(! Rt().collinear_3_object()(p,q,r) );
    
    typename Rt::Construct_weighted_circumcenter_3 weighted_circumcenter =
      Rt().construct_weighted_circumcenter_3_object();
    typename Rt::Has_on_bounded_side_3 on_bounded_side =
      Rt().has_on_bounded_side_3_object();
    
    // Compute denominator to swith to exact if it is 0
    const FT denom = compute_denom(p,q,r);
    if ( ! CGAL_NTS is_zero(denom) )
    {
      result_type point = weighted_circumcenter(p,q,r);
      
      // Fast output
      if ( on_bounded_side(Sphere_3(p,q,r),point) )
        return point;
    }
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();
    
    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r)));
  }
  
  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q ) const
  {
    typename Rt::Construct_weighted_circumcenter_3 weighted_circumcenter =
      Rt().construct_weighted_circumcenter_3_object();
    typename Rt::Has_on_bounded_side_3 on_bounded_side =
      Rt().has_on_bounded_side_3_object();
    
    // No division here
    result_type point = weighted_circumcenter(p,q);
    
    // Fast output
    if ( on_bounded_side(Sphere_3(p,q),point) )
      return point;
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();
    
    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q)));
  }

private:

  FT compute_denom(const Weighted_point_3 & p,
                   const Weighted_point_3 & q,
                   const Weighted_point_3 & r,
                   const Weighted_point_3 & s) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z(),
                         s.x(),s.y(),s.z());
  }

  FT compute_denom(const Weighted_point_3 & p,
                   const Weighted_point_3 & q,
                   const Weighted_point_3 & r) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z());
  }

  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &rx, const FT &ry, const FT &rz,
                   const FT &sx, const FT &sy, const FT &sz) const
  {
    const FT qpx = qx-px;
    const FT qpy = qy-py;
    const FT qpz = qz-pz;
    const FT rpx = rx-px;
    const FT rpy = ry-py;
    const FT rpz = rz-pz;
    const FT spx = sx-px;
    const FT spy = sy-py;
    const FT spz = sz-pz;

    return determinant(qpx,qpy,qpz,
                       rpx,rpy,rpz,
                       spx,spy,spz);
  }

  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &rx, const FT &ry, const FT &rz) const
  {
    const FT qpx = qx-px;
    const FT qpy = qy-py;
    const FT qpz = qz-pz;
    const FT rpx = rx-px;
    const FT rpy = ry-py;
    const FT rpz = rz-pz;
    const FT sx = qpy*rpz-qpz*rpy;
    const FT sy = qpz*rpx-qpx*rpz;
    const FT sz = qpx*rpy-qpy*rpx;

    return determinant(qpx,qpy,qpz,
                       rpx,rpy,rpz,
                       sx,sy,sz);
  }

};

/**
 * @class Robust_weighted_circumcenter_filtered_traits_3
 */
template<class K_>
struct Robust_weighted_circumcenter_filtered_traits_3
: public CGAL::Regular_triangulation_euclidean_traits_3<K_>
{
  typedef CGAL::Robust_filtered_construct_weighted_circumcenter_3<K_>
    Construct_weighted_circumcenter_3;
  
  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(); }
  
};  // end class Robust_weighted_circumcenter_filtered_traits_3


}  // end namespace CGAL

#endif // CGAL_MESH_3_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
