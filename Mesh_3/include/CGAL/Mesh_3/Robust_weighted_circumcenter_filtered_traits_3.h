// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
#include <CGAL/constructions/constructions_on_weighted_points_cartesian_3.h>

namespace CGAL {

template < typename K >
class Robust_filtered_compute_squared_radius_3
  : public K::Compute_squared_radius_3
{
public:
  typedef Exact_predicates_exact_constructions_kernel         EK;
  typedef Weighted_converter_3<Cartesian_converter<K, EK> >   To_exact;
  typedef Weighted_converter_3<Cartesian_converter<EK,K> >    Back_from_exact;
  
  typedef typename K::Point_3                         Point_3;
  typedef typename K::FT                              FT;
  typedef FT                                          result_type;

#ifndef  CGAL_CFG_MATCHING_BUG_6 
  using K::Compute_squared_radius_3::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  typedef typename K::Sphere_3                        Sphere_3;
  typedef typename K::Circle_3                        Circle_3;

  result_type
  operator()( const Sphere_3& s) const
  { return K::Compute_squared_radius_3::operator()(s); }

  result_type
  operator()( const Circle_3& c) const
  { return K::Compute_squared_radius_3::operator()(c); }

  FT operator() ( const Point_3 & p,
                  const Point_3 & q,
                  const Point_3 & r) const
  {
    return K::Compute_squared_radius_3::operator()(p,q,r);
  }

  FT operator() ( const Point_3 & p,
                  const Point_3 & q) const
  {
    return K::Compute_squared_radius_3::operator()(p,q);
  }

  FT operator() ( const Point_3 & p) const
  {
    return K::Compute_squared_radius_3::operator()(p);
  }
#endif // CGAL_CFG_MATCHING_BUG_6
  
  FT operator() ( const Point_3 & p,
                  const Point_3 & q,
                  const Point_3 & r,
                  const Point_3 & s ) const
  {
    typename K::Compute_squared_radius_3 sq_radius =
      K().compute_squared_radius_3_object();
    
    // Compute denominator to swith to exact if it is 0
    const FT denom = compute_denom(p,q,r,s);
    if ( ! CGAL_NTS is_zero(denom) )
    {
      return sq_radius(p,q,r,s);
    }
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Compute_squared_radius_3 exact_sq_radius =
      EK().compute_squared_radius_3_object();
    
    return back_from_exact(exact_sq_radius(to_exact(p),
                                           to_exact(q),
                                           to_exact(r),
                                           to_exact(s)));
  }

private:
  FT compute_denom(const Point_3 & p,
                   const Point_3 & q,
                   const Point_3 & r,
                   const Point_3 & s) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z(),
                         s.x(),s.y(),s.z());
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
};

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
    CGAL_precondition(Rt().orientation_3_object()(p,q,r,s) == CGAL::POSITIVE);
    
    // We use Side_of_oriented_sphere_3: it is static filtered and
    // we know that p,q,r,s are positive oriented
    typename Rt::Side_of_oriented_sphere_3 side_of_oriented_sphere =
      Rt().side_of_oriented_sphere_3_object();

    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             s.x(), s.y(), s.z(), s.weight(),
                                             num_x,  num_y, num_z, den);
    
    if ( ! CGAL_NTS is_zero(den) )
    {
      FT inv = FT(1)/(FT(2) * den);
      Bare_point res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);
      
      // Fast output
      if ( side_of_oriented_sphere(p,q,r,s,res) == CGAL::ON_POSITIVE_SIDE )
        return res;
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
        
    typename Rt::Side_of_bounded_sphere_3 side_of_bounded_sphere =
      Rt().side_of_bounded_sphere_3_object();
    
    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             num_x,  num_y, num_z, den);
    
    if ( ! CGAL_NTS is_zero(den) )
    {
      FT inv = FT(1)/(FT(2) * den);
      Bare_point res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);
      
      // Fast output
      if ( side_of_bounded_sphere(p,q,r,res) == CGAL::ON_BOUNDED_SIDE )
        return res;
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
    
    typename Rt::Side_of_bounded_sphere_3 side_of_bounded_sphere =
      Rt().side_of_bounded_sphere_3_object();
    
    // No division here
    result_type point = weighted_circumcenter(p,q);
    
    // Fast output
    if ( side_of_bounded_sphere(p,q,point) == CGAL::ON_BOUNDED_SIDE )
      return point;
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();
    
    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q)));
  }
};
  
  
  
template < typename K_ >
class Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3
{
public:
  typedef Exact_predicates_exact_constructions_kernel          EK;
  typedef Weighted_converter_3<Cartesian_converter<K_, EK> >   To_exact;
  typedef Weighted_converter_3<Cartesian_converter<EK,K_> >    Back_from_exact;
  
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K_> Rt;
  typedef CGAL::Regular_triangulation_euclidean_traits_3<EK> Exact_Rt;
  
  typedef typename Rt::Weighted_point_3               Weighted_point_3;
  typedef typename Rt::FT                             FT;
  
  FT operator() ( const Weighted_point_3& p,
                  const Weighted_point_3& q,
                  const Weighted_point_3& r,
                  const Weighted_point_3& s ) const
  {
    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             s.x(), s.y(), s.z(), s.weight(),
                                             num_x,  num_y, num_z, den);
    if ( ! CGAL_NTS is_zero(den) )
    {
      FT inv = FT(1)/(FT(2) * den);
      
      return (  CGAL_NTS square(num_x) 
              + CGAL_NTS square(num_y)
              + CGAL_NTS square(num_z) ) * CGAL_NTS square(inv) - p.weight();
    }
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
      Exact_Rt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
    
    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q),
                                                to_exact(r),to_exact(s)));
  }
  
  
  FT operator() (const Weighted_point_3& p,
                 const Weighted_point_3& q,
                 const Weighted_point_3& r) const
  {
    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             num_x,  num_y, num_z, den);
    
    if ( ! CGAL_NTS is_zero(den) )
    {
      FT inv = FT(1)/(FT(2) * den);
      
      return (  CGAL_NTS square(num_x) 
              + CGAL_NTS square(num_y)
              + CGAL_NTS square(num_z) ) * CGAL_NTS square(inv) - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
      Exact_Rt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
    
    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q),to_exact(r)));
  }
  
  
  FT operator() (const Weighted_point_3& p,
                 const Weighted_point_3& q) const
  {
    // Compute denominator to swith to exact if it is 0
    FT qpx = q.x() - p.x();
    FT qpy = q.y() - p.y();
    FT qpz = q.z() - p.z();
    FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
    
    if ( ! CGAL_NTS is_zero(qp2) )
    {
      FT inv = FT(1)/(FT(2)*qp2);
      FT alpha = 1/FT(2) + (p.weight()-q.weight())*inv;
      
      return  CGAL_NTS square(alpha)*qp2 - p.weight();      
    }
    
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    Exact_Rt::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
      Exact_Rt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
    
    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q)));
  }
  
  
  FT operator() (const Weighted_point_3& p) const
  {
    return -p.weight();
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
  
  typedef CGAL::Robust_filtered_compute_squared_radius_3<K_> Compute_squared_radius_3;
  
  typedef CGAL::Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3<K_>
    Compute_squared_radius_smallest_orthogonal_sphere_3;

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(); }

  Compute_squared_radius_3
  compute_squared_radius_3_object() const
  { return Compute_squared_radius_3(); }
  
  Compute_squared_radius_smallest_orthogonal_sphere_3
  compute_squared_radius_smallest_orthogonal_sphere_3_object() const
  { return Compute_squared_radius_smallest_orthogonal_sphere_3(); }
  
};  // end class Robust_weighted_circumcenter_filtered_traits_3


}  // end namespace CGAL

#endif // CGAL_MESH_3_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
