// Copyright (c) 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ROBUST_PERIODIC_WEIGHTED_CIRCUMCENTER_TRAITS_3_H
#define CGAL_ROBUST_PERIODIC_WEIGHTED_CIRCUMCENTER_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/Exact_kernel_selector.h>

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

namespace CGAL
{

template <typename K_, typename CWC_Base_>
class Robust_periodic_construct_weighted_circumcenter_3
  : public CWC_Base_
{
  typedef K_                                                 Base_traits;
  typedef CWC_Base_                                          Base;

public:
  typedef typename Base_traits::FT                           FT;
  typedef typename Base_traits::Point_3                      Point_3;
  typedef typename Base_traits::Weighted_point_3             Weighted_point_3;
  typedef typename Base_traits::Sphere_3                     Sphere_3;

  typedef typename Base_traits::Offset                       Offset;

  typedef Point_3                                            result_type;

  Robust_periodic_construct_weighted_circumcenter_3(const Base_traits& base_traits)
    : Base(base_traits.construct_weighted_circumcenter_3_object()),
      base_traits(base_traits)
  { }

  Point_3 operator()(const Weighted_point_3& p, const Weighted_point_3& q,
                     const Weighted_point_3& r, const Weighted_point_3& s) const
  {
    return Base::operator()(p,q,r,s);
  }

  Point_3 operator()(const Weighted_point_3& p, const Weighted_point_3& q, const Weighted_point_3& r) const
  {
    return Base::operator()(p,q,r);
  }

  Point_3 operator()(const Weighted_point_3& p, const Weighted_point_3& q) const
  {
    return Base::operator()(p,q);
  }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q,
                     const Weighted_point_3& r,
                     const Weighted_point_3& s,
                     const Offset& o_p, const Offset& o_q,
                     const Offset& o_r, const Offset& o_s) const
  {
    typename Base_traits::Construct_weighted_point_3 p2wp =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Power_side_of_oriented_power_sphere_3 ps =
               base_traits.power_side_of_oriented_power_sphere_3_object();

    Point_3 c = cwc(p, q, r, s, o_p, o_q, o_r, o_s);

    if(ps(p, q, r, s, p2wp(c), o_p, o_q, o_r, o_s, Offset()) != ON_POSITIVE_SIDE)
    {
      // switch to exact
      typedef typename Base_traits::Kernel                         K;
      typedef typename Exact_kernel_selector<K>::Exact_kernel      EK;
      typedef typename Exact_kernel_selector<K>::C2E               C2E;
      typedef typename Exact_kernel_selector<K>::E2C               E2C;

      C2E to_exact;
      E2C back_from_exact;

      typedef Periodic_3_regular_triangulation_traits_3<EK> Exact_traits;
      Exact_traits etraits(to_exact(base_traits.get_domain()));

      c = back_from_exact(
            etraits.construct_weighted_circumcenter_3_object()(
              to_exact(p), to_exact(q), to_exact(r), to_exact(s),
              o_p, o_q, o_r, o_s));

      CGAL_assertion(ps(p, q, r, s, p2wp(c), o_p, o_q, o_r, o_s, Offset()) == ON_POSITIVE_SIDE);
    }

    return c;
  }

  Point_3 operator()(const Weighted_point_3 & p,
                     const Weighted_point_3 & q,
                     const Weighted_point_3 & r,
                     const Offset& o_p, const Offset& o_q,
                     const Offset& o_r) const
  {
    typename Base_traits::Construct_weighted_point_3 p2wp =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Power_side_of_bounded_power_sphere_3 ps =
               base_traits.power_side_of_bounded_power_sphere_3_object();

    Point_3 c = cwc(p, q, r, o_p, o_q, o_r);

    if(ps(p, q, r, p2wp(c), o_p, o_q, o_r, Offset()) != ON_BOUNDED_SIDE)
    {
      // switch to exact
      typedef typename Base_traits::Kernel                         K;
      typedef typename Exact_kernel_selector<K>::Exact_kernel      EK;
      typedef typename Exact_kernel_selector<K>::C2E               C2E;
      typedef typename Exact_kernel_selector<K>::E2C               E2C;

      C2E to_exact;
      E2C back_from_exact;

      typedef Periodic_3_regular_triangulation_traits_3<EK> Exact_traits;
      Exact_traits etraits(to_exact(base_traits.get_domain()));

      c = back_from_exact(
            etraits.construct_weighted_circumcenter_3_object()(
              to_exact(p), to_exact(q), to_exact(r), o_p, o_q, o_r));

      CGAL_assertion(ps(p, q, r, p2wp(c), o_p, o_q, o_r, Offset()) == ON_BOUNDED_SIDE);
    }

    return c;
  }

  Point_3 operator()(const Weighted_point_3 & p,
                     const Weighted_point_3 & q,
                     const Offset& o_p, const Offset& o_q) const
  {
    typename Base_traits::Construct_weighted_point_3 p2wp =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Power_side_of_bounded_power_sphere_3 ps =
               base_traits.power_side_of_bounded_power_sphere_3_object();

    Point_3 c = cwc(p, q, o_p, o_q);

    if(ps(p, q, p2wp(c), o_p, o_q, Offset()) != ON_BOUNDED_SIDE)
    {
      // switch to exact
      typedef typename Base_traits::Kernel                         K;
      typedef typename Exact_kernel_selector<K>::Exact_kernel      EK;
      typedef typename Exact_kernel_selector<K>::C2E               C2E;
      typedef typename Exact_kernel_selector<K>::E2C               E2C;

      C2E to_exact;
      E2C back_from_exact;

      typedef Periodic_3_regular_triangulation_traits_3<EK> Exact_traits;
      Exact_traits etraits(to_exact(base_traits.get_domain()));

      c = back_from_exact(
            etraits.construct_weighted_circumcenter_3_object()(
              to_exact(p), to_exact(q), o_p, o_q));

      CGAL_assertion(ps(p, q, p2wp(c), o_p, o_q, Offset()) == ON_BOUNDED_SIDE);
    }

    return c;
  }

  const Base_traits& base_traits;
};

template <typename K_>
class Robust_periodic_weighted_circumcenter_traits_3
  : public K_
{
  typedef K_                                          Base_traits;

public:
  typedef typename Base_traits::Iso_cuboid_3          Iso_cuboid_3;

  typedef CGAL::Robust_periodic_construct_weighted_circumcenter_3<
            Base_traits, typename Base_traits::Construct_weighted_circumcenter_3>
                                                      Construct_weighted_circumcenter_3;

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(static_cast<const Base_traits&>(*this)); }

  Robust_periodic_weighted_circumcenter_traits_3(const Iso_cuboid_3& domain = Iso_cuboid_3(0,0,0,1,1,1),
                                                 const Base_traits& t = Base_traits())
    : Base_traits(domain, t)
  { }
};

}  // end namespace CGAL

#endif // CGAL_ROBUST_PERIODIC_WEIGHTED_CIRCUMCENTER_TRAITS_3_H
