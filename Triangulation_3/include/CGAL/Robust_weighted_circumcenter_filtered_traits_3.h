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

#ifndef CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
#define CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/constructions/kernel_ftC3.h>

namespace CGAL
{

template < typename K >
class Robust_filtered_compute_squared_radius_3
    : public K::Compute_squared_radius_3
{
public:
  typedef Exact_predicates_exact_constructions_kernel  EK;
  typedef Cartesian_converter<K, EK>                   To_exact;
  typedef Cartesian_converter<EK,K>                    Back_from_exact;

  typedef typename K::Point_3                          Point_3;
  typedef typename K::FT                               FT;
  typedef FT                                           result_type;

#ifndef  CGAL_CFG_MATCHING_BUG_6
  using K::Compute_squared_radius_3::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Circle_3                         Circle_3;

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
  typedef K_                                                 K;
  typedef Exact_predicates_exact_constructions_kernel        EK;
  typedef Cartesian_converter<K, EK>                         To_exact;
  typedef Cartesian_converter<EK, K>                         Back_from_exact;

  typedef typename K::Weighted_point_3                       Weighted_point_3;
  typedef typename K::Point_3                                Point_3;
  typedef typename K::FT                                     FT;
  typedef typename K::Sphere_3                               Sphere_3;

  typedef Point_3                                           result_type;

  typename K::Construct_point_3 wp2p;
  typename K::Construct_weighted_point_3 p2wp;

  Robust_filtered_construct_weighted_circumcenter_3()
    :
      wp2p(K().construct_point_3_object()),
      p2wp(K().construct_weighted_point_3_object())
  { }

  Point_3  operator() ( const Weighted_point_3 & p,
                        const Weighted_point_3 & q,
                        const Weighted_point_3 & r,
                        const Weighted_point_3 & s,
                        bool force_exact = false) const
  {
    CGAL_precondition(K().orientation_3_object()(
      wp2p(p), wp2p(q), wp2p(r), wp2p(s)) == CGAL::POSITIVE);

    if(! force_exact)
    {
      // We use power_side_of_power_sphere_3: it is static filtered and
      // we know that p,q,r,s are positive oriented
      typename K::Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere =
        K().power_side_of_oriented_power_sphere_3_object();

      // Compute denominator to swith to exact if it is 0
      FT num_x, num_y, num_z, den;
      bool unweighted = (p.weight() == 0) && (q.weight() == 0) &&
                        (r.weight() == 0) && (s.weight() == 0);

      if(unweighted){
        determinants_for_circumcenterC3(p.x(), p.y(), p.z(),
                                        q.x(), q.y(), q.z(),
                                        r.x(), r.y(), r.z(),
                                        s.x(), s.y(), s.z(),
                                        num_x,  num_y, num_z, den);
      } else {
        determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                                 q.x(), q.y(), q.z(), q.weight(),
                                                 r.x(), r.y(), r.z(), r.weight(),
                                                 s.x(), s.y(), s.z(), s.weight(),
                                                 num_x,  num_y, num_z, den);
      }

      if ( ! CGAL_NTS is_zero(den) )
      {
        FT inv = FT(1)/(FT(2) * den);
        Point_3 res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);

        if(unweighted){
          if (side_of_oriented_sphere(wp2p(p), wp2p(q), wp2p(r), wp2p(s), res)
              == CGAL::ON_POSITIVE_SIDE )
            return res;
        } else {
          // Fast output
          if ( power_side_of_oriented_power_sphere(p,q,r,s,p2wp(res)) == CGAL::ON_POSITIVE_SIDE )
            return res;
        }
      }
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r),
                                                       to_exact(s)));
  }

  Point_3 operator() ( const Weighted_point_3 & p,
                       const Weighted_point_3 & q,
                       const Weighted_point_3 & r ) const
  {
    CGAL_precondition(! K().collinear_3_object()(wp2p(p), wp2p(q), wp2p(r)));

    typename K::Side_of_bounded_sphere_3 side_of_bounded_sphere =
        K().side_of_bounded_sphere_3_object();

    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             num_x,  num_y, num_z, den);

    if ( ! CGAL_NTS is_zero(den) )
    {
      FT inv = FT(1)/(FT(2) * den);
      Point_3 res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);

      // Fast output
      if ( side_of_bounded_sphere(wp2p(p),wp2p(q),wp2p(r),res) == CGAL::ON_BOUNDED_SIDE )
        return res;
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r)));
  }

  Point_3 operator() ( const Weighted_point_3 & p,
                       const Weighted_point_3 & q ) const
  {
    typename K::Construct_weighted_circumcenter_3 weighted_circumcenter =
        K().construct_weighted_circumcenter_3_object();
    typename K::Side_of_bounded_sphere_3 side_of_bounded_sphere =
        K().side_of_bounded_sphere_3_object();
    typename K::Construct_point_3 cp =
        K().construct_point_3_object();

    // No division here
    result_type point = weighted_circumcenter(p,q);

    // Fast output
    if ( side_of_bounded_sphere(cp(p), cp(q), point) == CGAL::ON_BOUNDED_SIDE )
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q)));
  }
};

template < typename K_ >
class Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3
{
public:
  typedef K_                                                 K;
  typedef Exact_predicates_exact_constructions_kernel        EK;
  typedef Cartesian_converter<K_, EK>                        To_exact;
  typedef Cartesian_converter<EK,K_>                         Back_from_exact;

  typedef typename K::Weighted_point_3                       Weighted_point_3;
  typedef typename K::FT                                     FT;

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
    EK::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
        EK().compute_squared_radius_smallest_orthogonal_sphere_3_object();

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
    EK::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
        EK().compute_squared_radius_smallest_orthogonal_sphere_3_object();

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
    EK::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
        EK().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q)));
  }

  FT operator() (const Weighted_point_3& p) const
  {
    return -p.weight();
  }
};

template<class K_>
class Robust_weighted_circumcenter_filtered_traits_3
  : public K_
{
public:
  typedef CGAL::Robust_filtered_construct_weighted_circumcenter_3<K_>
    Construct_weighted_circumcenter_3;
  typedef CGAL::Robust_filtered_compute_squared_radius_3<K_>
    Compute_squared_radius_3;
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

  Robust_weighted_circumcenter_filtered_traits_3() { }
  Robust_weighted_circumcenter_filtered_traits_3(const K_& k) : K_(k) { }
};

}  // end namespace CGAL

#endif // CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
