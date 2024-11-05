// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#pragma once
#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Frechet_distance/internal/geometry_basics.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

// let a := v.x()^2 + v.y()^2,
// let b := line_start.x() * v.x() + line_start.y() * v.y(),
// let c := line_start.x()^2 + line_start.y()^2 - radius^2
// <=> lambda^2 * a + lambda * 2 b + c = 0
// <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
// <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
//TODO: use optional!
template <class Traits, class C, class Point, class AFT>
bool
fill_lambda(const Point& circle_center,
            const Point& line_start,
            const Point& line_end,
            const AFT& radius,
            std::pair<Lambda<C>, Lambda<C>>& I,
            C const& curve1,
            typename C::PointID const& center_id,
            C const& curve2,
            typename C::PointID const& seg_start_id)
{
  using FT = typename Traits::FT;
  FT a(0), b(0), c(0);
  auto ccci = typename Traits::Construct_cartesian_const_iterator_d();
  auto it_cc = ccci(circle_center), it_s = ccci(line_start), it_e = ccci(line_end);

  for (auto i = 0; i < C::dimension; ++i, ++it_cc, ++it_s, ++it_e)
  {
    FT start_end_diff = *it_e - *it_s;
    a += CGAL::square(start_end_diff);
    FT start_center_diff = *it_s - *it_cc;
    b -= start_center_diff * start_end_diff;
    c += square(start_center_diff);
  }
  c -= CGAL::square(radius);

  FT minus_b_div_a = b / a;
  FT d = CGAL::square(minus_b_div_a) - c / a;
  if (! is_negative(d))
  {
    auto start = make_root_of_2(minus_b_div_a, -1, d);
    auto end = make_root_of_2(minus_b_div_a, 1, d);
    // Question: can it be negative >1 even in exact?
    if (is_negative(start)) start = decltype(start)(FT(0));
    if (end > FT(1)) end =decltype(end)( FT(1));
    if (start <= FT(1) && end >= FT(0))
    {
      I = std::make_pair(
          (start == FT(0)) ? Lambda<C>(0)
                           : Lambda<C>(start, curve1, center_id,
                                       curve2, seg_start_id, radius, true),
          (end == FT(1)) ? Lambda<C>(1)
                         : Lambda<C>(end, curve1, center_id,
                                     curve2, seg_start_id, radius, false));
      return true;
    }
  }
  return false;
}

namespace HLPred
{

// filtered version
template <typename T>
Interval<Curve<T, true>>
intersection_interval(Curve<T, true> const& curve1,
                      typename Curve<T, true>::PointID const& center_id,
                      Curve<T, true> const& curve2,
                      typename Curve<T, true>::PointID seg_start_id,
                      typename Curve<T, true>::IFT const& radius)
{
  using C = Curve<T, true>;
  Interval<C> I;

  try {
    std::pair<Lambda<C>, Lambda<C>> II;
    // if not empty
    if (fill_lambda<typename T::first_type>(
        curve1[center_id], curve2[seg_start_id], curve2[seg_start_id + 1],
        radius.inf(), II, curve1, center_id, curve2,  seg_start_id))
    {
      I = Interval<C>(II.first, II.second);
    }
  } catch (const CGAL::Uncertain_conversion_exception& e) {
    CGAL_assertion(radius.is_point());
    std::pair<Lambda<C>, Lambda<C>> II;
      // if not empty
    if (fill_lambda<typename T::second_type>(
          curve1.rpoint(center_id), curve2.rpoint(seg_start_id), curve2.rpoint(seg_start_id + 1),
          radius.inf(), II, curve1, center_id, curve2,  seg_start_id))
    {
      I = Interval<C>(II.first, II.second);
    }
  }

  return I;
}

// non-filtered version
template <typename T>
Interval<Curve<T, false>>
intersection_interval(Curve<T, false> const& curve1,
                      typename Curve<T, false>::PointID const& center_id,
                      Curve<T, false> const& curve2,
                      typename Curve<T, false>::PointID seg_start_id,
                      typename Curve<T, false>::IFT const& radius)
{
  using C = Curve<T, false>;

  Interval<C> I;

    std::pair<Lambda<C>, Lambda<C>> II;
    // if not empty
    if (fill_lambda<T>(
        curve1[center_id], curve2[seg_start_id], curve2[seg_start_id + 1],
        radius.inf(), II, curve1, center_id, curve2,  seg_start_id))
    {
      I = Interval<C>(II.first, II.second);
    }

  return I;
}

}  // namespace HLPred

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL
