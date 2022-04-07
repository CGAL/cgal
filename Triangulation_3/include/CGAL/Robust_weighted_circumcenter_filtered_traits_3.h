// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

namespace CGAL {

template <typename K_, typename CSR_Base_>
class Robust_filtered_compute_squared_radius_3
  : public CSR_Base_
{
  typedef CSR_Base_                                    Base;

public:
  typedef K_                                           Kernel;
  typedef Exact_predicates_exact_constructions_kernel  EKernel;
  typedef Cartesian_converter<Kernel, EKernel>         To_exact;
  typedef Cartesian_converter<EKernel, Kernel>         Back_from_exact;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point_3;

  typedef FT                                           result_type;

  Robust_filtered_compute_squared_radius_3(const Kernel& k)
    : Base(k.compute_squared_radius_3_object()),
      traits(k)
  { }

  using Base::operator();

  FT operator()(const Point_3& p,
                const Point_3& q,
                const Point_3& r,
                const Point_3& s) const
  {
    typename Kernel::Compute_squared_radius_3 sq_radius =
      traits.compute_squared_radius_3_object();

    // Compute denominator to swith to exact if it is 0
    const FT denom = compute_denom(p,q,r,s);
    if( ! CGAL_NTS is_zero(denom) )
    {
      return sq_radius(p,q,r,s);
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_3 exact_sq_radius =
        EKernel().compute_squared_radius_3_object();

    return back_from_exact(exact_sq_radius(to_exact(p),
                                           to_exact(q),
                                           to_exact(r),
                                           to_exact(s)));
  }

private:
  FT compute_denom(const Point_3& p,
                   const Point_3& q,
                   const Point_3& r,
                   const Point_3& s) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z(),
                         s.x(),s.y(),s.z());
  }

  FT compute_denom(const FT& px, const FT& py, const FT& pz,
                   const FT& qx, const FT& qy, const FT& qz,
                   const FT& rx, const FT& ry, const FT& rz,
                   const FT& sx, const FT& sy, const FT& sz) const
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

private:
  const Kernel& traits;
};

template <typename K_, typename CWC_Base_>
class Robust_filtered_construct_weighted_circumcenter_3
  : public CWC_Base_
{
  typedef CWC_Base_                                         Base;

public:
  typedef K_                                                Kernel;
  typedef Exact_predicates_exact_constructions_kernel       EKernel;
  typedef Cartesian_converter<Kernel, EKernel>              To_exact;
  typedef Cartesian_converter<EKernel, Kernel>              Back_from_exact;

  typedef typename Kernel::Weighted_point_3                 Weighted_point_3;
  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::FT                               FT;
  typedef typename Kernel::Sphere_3                         Sphere_3;

  typedef Point_3                                           result_type;

  Robust_filtered_construct_weighted_circumcenter_3(const Kernel& k)
    : Base(k.construct_weighted_circumcenter_3_object()),
      traits(k)
  { }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q,
                     const Weighted_point_3& r,
                     const Weighted_point_3& s,
                     bool force_exact = false) const
  {
    typename Kernel::Construct_point_3 cp = traits.construct_point_3_object();
    typename Kernel::Construct_weighted_point_3 cwp = traits.construct_weighted_point_3_object();

    CGAL_precondition(Kernel().orientation_3_object()(
      cp(p), cp(q), cp(r), cp(s)) == CGAL::POSITIVE);

    if(! force_exact)
    {
      // We use power_side_of_power_sphere_3: it is static filtered and
      // we know that p,q,r,s are positive oriented
      typename Kernel::Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere =
        traits.power_side_of_oriented_power_sphere_3_object();

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

      if(! CGAL_NTS is_zero(den))
      {
        FT inv = FT(1)/(FT(2) * den);
        Point_3 res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);

        if(unweighted)
        {
          if(side_of_oriented_sphere(cp(p), cp(q), cp(r), cp(s), res)
              == CGAL::ON_POSITIVE_SIDE )
            return res;
        }
        else
        {
          // Fast output
          if(power_side_of_oriented_power_sphere(p,q,r,s,cwp(res)) == CGAL::ON_POSITIVE_SIDE)
            return res;
        }
      }
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EKernel().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r),
                                                       to_exact(s)));
  }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q,
                     const Weighted_point_3& r) const
  {
    typename Kernel::Construct_point_3 cp = traits.construct_point_3_object();
    CGAL_precondition(! traits.collinear_3_object()(cp(p), cp(q), cp(r)));

    typename Kernel::Side_of_bounded_sphere_3 side_of_bounded_sphere =
      traits.side_of_bounded_sphere_3_object();

    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             num_x,  num_y, num_z, den);

    if(! CGAL_NTS is_zero(den))
    {
      FT inv = FT(1)/(FT(2) * den);
      Point_3 res(p.x() + num_x*inv, p.y() - num_y*inv, p.z() + num_z*inv);

      // Fast output
      if(side_of_bounded_sphere(cp(p),cp(q),cp(r),res) == CGAL::ON_BOUNDED_SIDE)
        return res;
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EKernel().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r)));
  }

  Point_3 operator()(const Weighted_point_3& p,
                     const Weighted_point_3& q) const
  {
    typename Kernel::Construct_point_3 cp =
        traits.construct_point_3_object();
    typename Kernel::Construct_weighted_circumcenter_3 weighted_circumcenter =
        traits.construct_weighted_circumcenter_3_object();
    typename Kernel::Side_of_bounded_sphere_3 side_of_bounded_sphere =
        traits.side_of_bounded_sphere_3_object();

    // No division here
    result_type point = weighted_circumcenter(p,q);

    // Fast output
    if(side_of_bounded_sphere(cp(p), cp(q), point) == CGAL::ON_BOUNDED_SIDE)
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
        EKernel().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p), to_exact(q)));
  }

private:
  const Kernel& traits;
};

template <typename K_, typename CSRSOS_Base_>
class Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3
  : public CSRSOS_Base_
{
  typedef CSRSOS_Base_                                       Base;

public:
  typedef K_                                                 Kernel;
  typedef Exact_predicates_exact_constructions_kernel        EKernel;
  typedef Cartesian_converter<Kernel, EKernel>               To_exact;
  typedef Cartesian_converter<EKernel, Kernel>               Back_from_exact;

  typedef typename Kernel::Weighted_point_3                  Weighted_point_3;
  typedef typename Kernel::FT                                FT;

  typedef FT                                                 result_type;

  Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3(const Kernel& k)
    : Base(k.compute_squared_radius_smallest_orthogonal_sphere_3_object())
  { }

  FT operator()(const Weighted_point_3& p,
                const Weighted_point_3& q,
                const Weighted_point_3& r,
                const Weighted_point_3& s) const
  {
    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             s.x(), s.y(), s.z(), s.weight(),
                                             num_x,  num_y, num_z, den);
    if(! CGAL_NTS is_zero(den))
    {
      FT inv = FT(1)/(FT(2) * den);

      return (CGAL_NTS square(num_x) +
              CGAL_NTS square(num_y) +
              CGAL_NTS square(num_z)) * CGAL_NTS square(inv) - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
      EKernel().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q),
                                                to_exact(r),to_exact(s)));
  }

  FT operator()(const Weighted_point_3& p,
                const Weighted_point_3& q,
                const Weighted_point_3& r) const
  {
    // Compute denominator to swith to exact if it is 0
    FT num_x, num_y, num_z, den;
    determinants_for_weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                             q.x(), q.y(), q.z(), q.weight(),
                                             r.x(), r.y(), r.z(), r.weight(),
                                             num_x,  num_y, num_z, den);

    if(! CGAL_NTS is_zero(den))
    {
      FT inv = FT(1)/(FT(2) * den);

      return (CGAL_NTS square(num_x) +
              CGAL_NTS square(num_y) +
              CGAL_NTS square(num_z) ) * CGAL_NTS square(inv) - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
      EKernel().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q),to_exact(r)));
  }

  FT operator()(const Weighted_point_3& p,
                const Weighted_point_3& q) const
  {
    // Compute denominator to swith to exact if it is 0
    FT qpx = q.x() - p.x();
    FT qpy = q.y() - p.y();
    FT qpz = q.z() - p.z();
    FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);

    if(! CGAL_NTS is_zero(qp2))
    {
      FT inv = FT(1)/(FT(2)*qp2);
      FT alpha = 1/FT(2) + (p.weight()-q.weight())*inv;

      return  CGAL_NTS square(alpha)*qp2 - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_smallest_orthogonal_sphere_3 exact_compute_radius =
        EKernel().compute_squared_radius_smallest_orthogonal_sphere_3_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q)));
  }

  FT operator()(const Weighted_point_3& p) const
  {
    return -p.weight();
  }
};

template<class K_>
class Robust_weighted_circumcenter_filtered_traits_3
  : public K_
{
  typedef Robust_weighted_circumcenter_filtered_traits_3<K_>          Self;

public:
  typedef K_                                                          Kernel;

  typedef CGAL::Robust_filtered_construct_weighted_circumcenter_3<
            Kernel, typename Kernel::Construct_weighted_circumcenter_3>
                                                            Construct_weighted_circumcenter_3;

  typedef CGAL::Robust_filtered_compute_squared_radius_3<
            Kernel, typename Kernel::Compute_squared_radius_3>  Compute_squared_radius_3;

  typedef CGAL::Robust_filtered_compute_squared_radius_smallest_orthogonal_sphere_3<
            Kernel, typename Kernel::Compute_squared_radius_smallest_orthogonal_sphere_3>
                                                            Compute_squared_radius_smallest_orthogonal_sphere_3;

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(static_cast<const Kernel&>(*this)); }

  Compute_squared_radius_3
  compute_squared_radius_3_object() const
  { return Compute_squared_radius_3(static_cast<const Kernel&>(*this)); }

  Compute_squared_radius_smallest_orthogonal_sphere_3
  compute_squared_radius_smallest_orthogonal_sphere_3_object() const
  { return Compute_squared_radius_smallest_orthogonal_sphere_3(static_cast<const Kernel&>(*this)); }

  Robust_weighted_circumcenter_filtered_traits_3(const Kernel& k = Kernel()) : Kernel(k) { }
};

}  // end namespace CGAL

#endif // CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_3_H
