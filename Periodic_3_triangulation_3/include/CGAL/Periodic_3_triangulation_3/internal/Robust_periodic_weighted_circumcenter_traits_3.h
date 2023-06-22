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
    typename Base_traits::Compute_weight_3 weight =
               base_traits.compute_weight_3_object();
    typename Base_traits::Construct_point_3 point =
               base_traits.construct_point_3_object();
    typename Base_traits::Construct_weighted_point_3 weighted_point =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Orientation_3 orientation =
               base_traits.orientation_3_object();
    typename Base_traits::Power_side_of_oriented_power_sphere_3 ps =
               base_traits.power_side_of_oriented_power_sphere_3_object();

    // Calling Construct_weighted_circumcenter_3(p,q,r,s,o_p,o_q,o_r,o_s) will construct
    // the 3D weighted points (in Functor_with_offset_adaptor) but that is bad because
    // if we're unlucky pqrs is almost flat and the orientation might become wrong
    // (i.e. != POSITIVE) when we construct the Euclidean, non-periodic 3D representations
    // of the points
    const Point_3 euc_p = point(p, o_p);
    const Point_3 euc_q = point(q, o_q);
    const Point_3 euc_r = point(r, o_r);
    const Point_3 euc_s = point(s, o_s);

    bool needs_exact = false;
    Point_3 c;
    if(orientation(euc_p, euc_q, euc_r, euc_s) != POSITIVE) // see comment above
    {
      needs_exact = true;
    }
    else
    {
      c = cwc(weighted_point(euc_p, weight(p)),
              weighted_point(euc_q, weight(q)),
              weighted_point(euc_r, weight(r)),
              weighted_point(euc_s, weight(s)));
    }

    if(needs_exact ||
       ps(p, q, r, s, weighted_point(c), o_p, o_q, o_r, o_s, Offset()) != ON_POSITIVE_SIDE)
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

      CGAL_assertion(ps(p, q, r, s, weighted_point(c), o_p, o_q, o_r, o_s, Offset()) == ON_POSITIVE_SIDE);
    }

    return c;
  }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q,
                     const Weighted_point_3& r,
                     const Offset& o_p,
                     const Offset& o_q,
                     const Offset& o_r) const
  {
    typename Base_traits::Compute_weight_3 weight =
               base_traits.compute_weight_3_object();
    typename Base_traits::Construct_point_3 point =
               base_traits.construct_point_3_object();
    typename Base_traits::Construct_weighted_point_3 weighted_point =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Collinear_3 collinear =
               base_traits.collinear_3_object();
    typename Base_traits::Power_side_of_bounded_power_sphere_3 ps =
               base_traits.power_side_of_bounded_power_sphere_3_object();

    // Calling Construct_weighted_circumcenter_3(p,q,r,o_p,o_q,o_r) will construct
    // the 3D weighted points (in Functor_with_offset_adaptor) but that is bad because
    // if we're unlucky pqr is almost collinear and it becomes so when we construct the Euclidean,
    // non-periodic 3D representations of the points
    const Point_3 euc_p = point(p, o_p);
    const Point_3 euc_q = point(q, o_q);
    const Point_3 euc_r = point(r, o_r);

    bool needs_exact = false;
    Point_3 c;
    if(collinear(euc_p, euc_q, euc_r)) // see comment above
    {
      needs_exact = true;
    }
    else
    {
      c = cwc(weighted_point(euc_p, weight(p)),
              weighted_point(euc_q, weight(q)),
              weighted_point(euc_r, weight(r)));
    }

    if(needs_exact ||
       ps(p, q, r, weighted_point(c), o_p, o_q, o_r, Offset()) != ON_BOUNDED_SIDE)
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

      CGAL_assertion(ps(p, q, r, weighted_point(c), o_p, o_q, o_r, Offset()) == ON_BOUNDED_SIDE);
    }

    return c;
  }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q,
                     const Offset& o_p,
                     const Offset& o_q) const
  {
    typename Base_traits::Equal_3 equal =
               base_traits.equal_3_object();
    typename Base_traits::Compute_weight_3 weight =
               base_traits.compute_weight_3_object();
    typename Base_traits::Construct_point_3 point =
               base_traits.construct_point_3_object();
    typename Base_traits::Construct_weighted_point_3 weighted_point =
               base_traits.construct_weighted_point_3_object();
    typename Base_traits::Construct_weighted_circumcenter_3 cwc =
               base_traits.construct_weighted_circumcenter_3_object();
    typename Base_traits::Power_side_of_bounded_power_sphere_3 ps =
               base_traits.power_side_of_bounded_power_sphere_3_object();

    // Calling Construct_weighted_circumcenter_3(p,q,o_p,o_q) will construct
    // the 3D weighted points (in Functor_with_offset_adaptor) but that is bad because
    // if we're unlucky p and q become equal when we construct the Euclidean,
    // non-periodic 3D representations of the points
    const Point_3 euc_p = point(p, o_p);
    const Point_3 euc_q = point(q, o_q);

    bool needs_exact = false;
    Point_3 c;
    if(equal(euc_p, euc_q)) // see comment above
    {
      needs_exact = true;
    }
    else
    {
      c = cwc(weighted_point(euc_p, weight(p)),
              weighted_point(euc_q, weight(q)));
    }

    if(needs_exact ||
       ps(p, q, weighted_point(c), o_p, o_q, Offset()) != ON_BOUNDED_SIDE)
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

      CGAL_assertion(ps(p, q, weighted_point(c), o_p, o_q, Offset()) == ON_BOUNDED_SIDE);
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
