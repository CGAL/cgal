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

#ifndef ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H_
#define ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H_


#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Robust_construction.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>


namespace CGAL {


template < typename K >
class Robust_filtered_construct_weighted_circumcenter_3
{
public:
  typedef typename K::Weighted_point_3               Weighted_point_3;
  typedef typename K::Bare_point                     Bare_point;
  typedef typename K::FT                             FT;
  typedef typename K::Sphere_3                       Sphere_3;
  typedef Bare_point                                 result_type;

  typedef Exact_predicates_exact_constructions_kernel   EK2;
  typedef Regular_triangulation_euclidean_traits_3<EK2> EK;
  typedef Weighted_converter_3<
      Cartesian_converter<typename K::Kernel, EK2> >    To_exact;
  typedef Weighted_converter_3<
      Cartesian_converter<EK2, typename K::Kernel> >    Back_from_exact;


  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r,
                          const Weighted_point_3 & s ) const
  {
    typename K::Construct_weighted_circumcenter_3 weighted_circumcenter =
        K().construct_weighted_circumcenter_3_object();
    typename K::Has_on_bounded_side_3 on_bounded_side =
        K().has_on_bounded_side_3_object();

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
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                   to_exact(q),
                   to_exact(r),
                   to_exact(s)));
  }

  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r ) const
  {
    typename K::Construct_weighted_circumcenter_3 weighted_circumcenter =
      K().construct_weighted_circumcenter_3_object();
    typename K::Has_on_bounded_side_3 on_bounded_side =
      K().has_on_bounded_side_3_object();

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
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                   to_exact(q),
                   to_exact(r)));
  }

  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q ) const
  {
    typename K::Construct_weighted_circumcenter_3 weighted_circumcenter =
      K().construct_weighted_circumcenter_3_object();
    typename K::Has_on_bounded_side_3 on_bounded_side =
      K().has_on_bounded_side_3_object();

    // No division here
    result_type point = weighted_circumcenter(p,q);

    // Fast output
    if ( on_bounded_side(Sphere_3(p,q),point) )
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      EK().construct_weighted_circumcenter_3_object();

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
template<class K>
struct Robust_weighted_circumcenter_filtered_traits_3
: public K
{
  typedef CGAL::Robust_filtered_construct_weighted_circumcenter_3<K>
                                            Construct_weighted_circumcenter_3;

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(); }

};  // end class Robust_weighted_circumcenter_filtered_traits_3


}  // end namespace CGAL

#endif // ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H_
