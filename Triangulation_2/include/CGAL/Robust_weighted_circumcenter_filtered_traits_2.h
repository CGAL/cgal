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
#ifndef CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_2_H
#define CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/constructions/kernel_ftC2.h>

namespace CGAL {

template <typename K_, typename CC_Base_>
class Robust_filtered_construct_circumcenter_2
  : public CC_Base_
{
  typedef CC_Base_                                          Base;

public:
  typedef K_                                                Kernel;
  typedef Exact_predicates_exact_constructions_kernel       EKernel;
  typedef Cartesian_converter<Kernel, EKernel>              To_exact;
  typedef Cartesian_converter<EKernel, Kernel>              Back_from_exact;

  typedef typename Kernel::FT                               FT;
  typedef typename Kernel::Point_2                          Point_2;

  typedef Point_2                                           result_type;

  Robust_filtered_construct_circumcenter_2(const Kernel& k)
    : Base(k.construct_circumcenter_2_object()),
      traits(k)
  { }

  Point_2 operator()(const Point_2& p,
                     const Point_2& q,
                     const Point_2& r,
                     bool force_exact = false) const
  {
    CGAL_precondition(traits.orientation_2_object()(p, q, r) == CGAL::POSITIVE);

    if(! force_exact)
    {
      // Compute denominator to switch to exact if it is 0
      FT num_x, num_y, den;
      determinants_for_circumcenterC2(p.x(), p.y(),
                                      q.x(), q.y(),
                                      r.x(), r.y(),
                                      num_x, num_y, den);

      if(! CGAL_NTS is_zero(den))
      {
        FT inv = FT(1)/(FT(2) * den);
        Point_2 res(p.x() + num_x*inv, p.y() - num_y*inv);

        if(side_of_oriented_circle(p, q, r, res) == CGAL::ON_POSITIVE_SIDE)
          return res;
      }
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_circumcenter_2 exact_circumcenter = EKernel().construct_circumcenter_2_object();

    return back_from_exact(exact_circumcenter(to_exact(p), to_exact(q), to_exact(r)));
  }

  Point_2 operator()(const Point_2& p,
                     const Point_2& q) const
  {
    typename Kernel::Construct_circumcenter_2 circumcenter = traits.construct_circumcenter_2_object();
    typename Kernel::Side_of_bounded_circle_2 side_of_bounded_circle = traits.side_of_bounded_circle_2_object();

    // No division here
    Point_2 point = circumcenter(p,q);

    // Fast output
    if(side_of_bounded_circle(p, q, point) == CGAL::ON_BOUNDED_SIDE)
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_circumcenter_2 exact_circumcenter = EKernel().construct_circumcenter_2_object();

    return back_from_exact(exact_circumcenter(to_exact(p), to_exact(q)));
  }

private:
  const Kernel& traits;
};

template <typename K_, typename CSR_Base_>
class Robust_filtered_compute_squared_radius_2
  : public CSR_Base_
{
  typedef CSR_Base_                                    Base;

public:
  typedef K_                                           Kernel;
  typedef Exact_predicates_exact_constructions_kernel  EKernel;
  typedef Cartesian_converter<Kernel, EKernel>         To_exact;
  typedef Cartesian_converter<EKernel, Kernel>         Back_from_exact;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_2                     Point_2;

  typedef FT                                           result_type;

  Robust_filtered_compute_squared_radius_2(const Kernel& k)
    : Base(k.compute_squared_radius_2_object()),
      traits(k)
  { }

  using Base::operator();

  FT operator()(const Point_2& p,
                const Point_2& q,
                const Point_2& r) const
  {
    typename Kernel::Compute_squared_radius_2 sq_radius =
      traits.compute_squared_radius_2_object();

    // Compute denominator to switch to exact if it is 0
    const FT denom = compute_denom(p,q,r);
    if( ! CGAL_NTS is_zero(denom) )
    {
      return sq_radius(p,q,r);
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_2 exact_sq_radius = EKernel().compute_squared_radius_2_object();

    return back_from_exact(exact_sq_radius(to_exact(p),
                                           to_exact(q),
                                           to_exact(r)));
  }

private:
  FT compute_denom(const Point_2& p,
                   const Point_2& q,
                   const Point_2& r) const
  {
    return compute_denom(p.x(),p.y(),
                         q.x(),q.y(),
                         r.x(),r.y());
  }

  FT compute_denom(const FT& px, const FT& py,
                   const FT& qx, const FT& qy,
                   const FT& rx, const FT& ry) const
  {
    const FT qpx = qx-px;
    const FT qpy = qy-py;

    const FT rpx = rx-px;
    const FT rpy = ry-py;

    return determinant(qpx,qpy,
                       rpx,rpy);
  }

private:
  const Kernel& traits;
};

template <typename K_, typename CWC_Base_>
class Robust_filtered_construct_weighted_circumcenter_2
  : public CWC_Base_
{
  typedef CWC_Base_                                         Base;

public:
  typedef K_                                                Kernel;
  typedef Exact_predicates_exact_constructions_kernel       EKernel;
  typedef Cartesian_converter<Kernel, EKernel>              To_exact;
  typedef Cartesian_converter<EKernel, Kernel>              Back_from_exact;

  typedef typename Kernel::Weighted_point_2                 Weighted_point_2;
  typedef typename Kernel::Point_2                          Point_2;
  typedef typename Kernel::FT                               FT;
  typedef typename Kernel::Circle_2                         Circle_2;

  typedef Point_2                                           result_type;

  Robust_filtered_construct_weighted_circumcenter_2(const Kernel& k)
    : Base(k.construct_weighted_circumcenter_2_object()),
      traits(k)
  { }

  Point_2 operator()(const Weighted_point_2& p,
                     const Weighted_point_2& q,
                     const Weighted_point_2& r,
                     bool force_exact = false) const
  {
    typename Kernel::Construct_point_2 cp = traits.construct_point_2_object();
    typename Kernel::Construct_weighted_point_2 cwp = traits.construct_weighted_point_2_object();

    CGAL_precondition(traits.orientation_2_object()(cp(p), cp(q), cp(r)) == CGAL::POSITIVE);

    if(! force_exact)
    {
      // Compute denominator to switch to exact if it is 0
      FT num_x, num_y, den;
      bool unweighted = (p.weight() == 0) && (q.weight() == 0) && (r.weight() == 0);

      if(unweighted){
        determinants_for_circumcenterC2(p.x(), p.y(),
                                        q.x(), q.y(),
                                        r.x(), r.y(),
                                        num_x,  num_y, den);
      } else {
        determinants_for_weighted_circumcenterC2(p.x(), p.y(), p.weight(),
                                                 q.x(), q.y(), q.weight(),
                                                 r.x(), r.y(), r.weight(),
                                                 num_x, num_y, den);
      }

      if(! CGAL_NTS is_zero(den))
      {
        FT inv = FT(1)/(FT(2) * den);
        Point_2 res(p.x() + num_x*inv, p.y() - num_y*inv);

        if(unweighted)
        {
          typename Kernel::Side_of_oriented_circle_2 side_of_oriented_circle =
            traits.side_of_oriented_circle_2_object();

          if(side_of_oriented_circle(cp(p), cp(q), cp(r), res) == CGAL::ON_POSITIVE_SIDE)
            return res;
        }
        else
        {
          // We use power_side_of_oriented_circle_2: it is static filtered and
          // we know that p,q, and r are positively oriented
          typename Kernel::Power_side_of_oriented_power_circle_2 power_side_of_oriented_power_circle =
            traits.power_side_of_oriented_power_circle_2_object();

          // Fast output
          if(power_side_of_oriented_power_circle(p,q,r,cwp(res)) == CGAL::ON_POSITIVE_SIDE)
            return res;
        }
      }
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_weighted_circumcenter_2 exact_weighted_circumcenter =
        EKernel().construct_weighted_circumcenter_2_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r)));
  }

  Point_2 operator()(const Weighted_point_2& p,
                     const Weighted_point_2& q) const
  {
    typename Kernel::Construct_point_2 cp =
        traits.construct_point_2_object();
    typename Kernel::Construct_weighted_circumcenter_2 weighted_circumcenter =
        traits.construct_weighted_circumcenter_2_object();
    typename Kernel::Side_of_bounded_circle_2 side_of_bounded_circle =
        traits.side_of_bounded_circle_2_object();

    // No division here
    Point_2 point = weighted_circumcenter(p,q);

    // Fast output
    if(side_of_bounded_circle(cp(p), cp(q), point) == CGAL::ON_BOUNDED_SIDE)
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Construct_weighted_circumcenter_2 exact_weighted_circumcenter =
        EKernel().construct_weighted_circumcenter_2_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p), to_exact(q)));
  }

private:
  const Kernel& traits;
};

template <typename K_, typename CSRSOS_Base_>
class Robust_filtered_compute_squared_radius_smallest_orthogonal_circle_2
  : public CSRSOS_Base_
{
  typedef CSRSOS_Base_                                       Base;

public:
  typedef K_                                                 Kernel;
  typedef Exact_predicates_exact_constructions_kernel        EKernel;
  typedef Cartesian_converter<Kernel, EKernel>               To_exact;
  typedef Cartesian_converter<EKernel, Kernel>               Back_from_exact;

  typedef typename Kernel::Weighted_point_2                  Weighted_point_2;
  typedef typename Kernel::FT                                FT;

  typedef FT                                                 result_type;

  Robust_filtered_compute_squared_radius_smallest_orthogonal_circle_2(const Kernel& k)
    : Base(k.compute_squared_radius_smallest_orthogonal_circle_2_object())
  { }

  FT operator()(const Weighted_point_2& p,
                const Weighted_point_2& q,
                const Weighted_point_2& r) const
  {
    // Compute denominator to switch to exact if it is 0
    FT num_x, num_y, den;
    determinants_for_weighted_circumcenterC2(p.x(), p.y(), p.weight(),
                                             q.x(), q.y(), q.weight(),
                                             r.x(), r.y(), r.weight(),
                                             num_x,  num_y, den);
    if(! CGAL_NTS is_zero(den))
    {
      FT inv = FT(1)/(FT(2) * den);

      return (CGAL_NTS square(num_x) +
              CGAL_NTS square(num_y)) * CGAL_NTS square(inv) - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_smallest_orthogonal_circle_2 exact_compute_radius =
      EKernel().compute_squared_radius_smallest_orthogonal_circle_2_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q),to_exact(r)));
  }

  FT operator()(const Weighted_point_2& p,
                const Weighted_point_2& q) const
  {
    // Compute denominator to switch to exact if it is 0
    FT qpx = q.x() - p.x();
    FT qpy = q.y() - p.y();
    FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy);

    if(! CGAL_NTS is_zero(qp2))
    {
      FT inv = FT(1)/(FT(2)*qp2);
      FT alpha = 1/FT(2) + (p.weight()-q.weight())*inv;

      return  CGAL_NTS square(alpha)*qp2 - p.weight();
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EKernel::Compute_squared_radius_smallest_orthogonal_circle_2 exact_compute_radius =
        EKernel().compute_squared_radius_smallest_orthogonal_circle_2_object();

    return back_from_exact(exact_compute_radius(to_exact(p),to_exact(q)));
  }

  FT operator()(const Weighted_point_2& p) const
  {
    return -p.weight();
  }
};

template<class K_>
class Robust_circumcenter_filtered_traits_2
  : public K_
{
  typedef Robust_circumcenter_filtered_traits_2<K_>                   Self;

public:
  typedef K_                                                          Kernel;

  typedef CGAL::Robust_filtered_construct_circumcenter_2<
            Kernel, typename Kernel::Construct_circumcenter_2>
                                                                      Construct_circumcenter_2;

  typedef CGAL::Robust_filtered_compute_squared_radius_2<
            Kernel, typename Kernel::Compute_squared_radius_2>        Compute_squared_radius_2;

  Construct_circumcenter_2
  construct_circumcenter_2_object() const
  { return Construct_circumcenter_2(static_cast<const Kernel&>(*this)); }

  Compute_squared_radius_2
  compute_squared_radius_2_object() const
  { return Compute_squared_radius_2(static_cast<const Kernel&>(*this)); }

  Robust_circumcenter_filtered_traits_2(const Kernel& k = Kernel()) : Kernel(k) { }
};


template < class BaseGt >
struct Triangulation_structural_filtering_traits<Robust_circumcenter_filtered_traits_2<BaseGt> > {
    typedef typename Triangulation_structural_filtering_traits<BaseGt>::Use_structural_filtering_tag  Use_structural_filtering_tag;
};


template<class K_>
class Robust_weighted_circumcenter_filtered_traits_2
  : public Robust_circumcenter_filtered_traits_2<K_>
{
  typedef Robust_circumcenter_filtered_traits_2<K_>                   Base;
  typedef Robust_weighted_circumcenter_filtered_traits_2<Base>        Self;

public:
  typedef K_                                                          Kernel;

  typedef CGAL::Robust_filtered_construct_weighted_circumcenter_2<
            Kernel, typename Kernel::Construct_weighted_circumcenter_2>
                                                                      Construct_weighted_circumcenter_2;

  typedef CGAL::Robust_filtered_compute_squared_radius_smallest_orthogonal_circle_2<
            Kernel, typename Kernel::Compute_squared_radius_smallest_orthogonal_circle_2>
                                                                      Compute_squared_radius_smallest_orthogonal_circle_2;

  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
  { return Construct_weighted_circumcenter_2(static_cast<const Kernel&>(*this)); }

  Compute_squared_radius_smallest_orthogonal_circle_2
  compute_squared_radius_smallest_orthogonal_circle_2_object() const
  { return Compute_squared_radius_smallest_orthogonal_circle_2(static_cast<const Kernel&>(*this)); }

  Robust_weighted_circumcenter_filtered_traits_2(const Kernel& k = Kernel()) : Base(k) { }
};

} // namespace CGAL

#endif // CGAL_ROBUST_WEIGHTED_CIRCUMCENTER_FILTERED_TRAITS_2_H
