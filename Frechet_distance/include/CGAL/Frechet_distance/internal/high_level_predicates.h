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

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {



template <typename C>
bool approximate_reals(C const& curve1,
                       typename C::PointID const& center_id,
                       C const& curve2,
                       typename C::PointID const& seg_start_id,
                       typename C::distance_t const& radius,
                       std::pair<Lambda<C>, Lambda<C>>& I)
{
    using Point = typename C::Point;
    const Point& circle_center = curve1[center_id];
    const Point& line_start = curve2[seg_start_id];
    const Point& line_end  = curve2[seg_start_id + 1];

    // let a := v.x()^2 + v.y()^2,
    // let b := line_start.x() * v.x() + line_start.y() * v.y(),
    // let c := line_start.x()^2 + line_start.y()^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)

    using Approx = typename C::distance_t;
    Approx a(0), b(0), c(0);
    for (auto i = 0; i < C::dimension; ++i) {
        Approx start_end_diff = line_end[i] - line_start[i];
        a += CGAL::square(start_end_diff);
        Approx start_center_diff = line_start[i] - circle_center[i];
        b -= start_center_diff * start_end_diff;
        c += CGAL::square(start_center_diff);
    }
    c -= CGAL::square(radius);

    Approx minus_b_div_a = b / a;
    Approx d = CGAL::square(minus_b_div_a) - c / a;
    if (! is_negative(d)) {
        Approx sqrtd = sqrt(d);
        Approx start = minus_b_div_a - sqrtd;
        Approx end = minus_b_div_a + sqrtd;
        if (is_negative(start)) start = Approx(0);
        if (end > Approx(1)) end = Approx(1);
        if (start <= Approx(1) && end >= Approx(0)) {
            I = std::make_pair(
                (start == Approx(0)) ? Lambda<C>(0)
                                     : Lambda<C>(start, curve1, center_id,
                                              curve2, seg_start_id, radius, true),
                (end == Approx(1)) ? Lambda<C>(1)
                                   : Lambda<C>(end, curve1, center_id,
                                           curve2, seg_start_id, radius, false));
            return true;
        }
    }
    return false;
}

template <typename C>
bool exact_reals(C const& curve1,
                 typename C::PointID const& center_id,
                 C const& curve2,
                 typename C::PointID const& seg_start_id,
                 typename C::distance_t const& radius,
                 std::pair<Lambda<C>, Lambda<C>>& I)
{
  using Rational_point = typename C::Rational_point;

  const Rational_point& circle_center = curve1.rpoint(center_id);
  const Rational_point& line_start = curve2.rpoint(seg_start_id);
  const Rational_point& line_end  = curve2.rpoint(seg_start_id + 1);

    using Exact = typename Lambda<C>::Exact;
    using Rational = typename Lambda<C>::Rational;

    Rational a(0), b(0), c(0);
    for (auto i = 0; i < C::dimension; ++i) {
      Rational start_end_diff = line_end[i] - line_start[i];
        a += CGAL::square(start_end_diff);
        Rational start_center_diff = line_start[i] - circle_center[i];
        b -= start_center_diff * start_end_diff;
        c += CGAL::square(start_center_diff);
    }
    assert(radius.is_point());
    c -= CGAL::square(Rational(radius.inf()));

    Rational minus_b_div_a = b / a;
    Rational d = CGAL::square(minus_b_div_a) - c / a;
    if (! is_negative(d)) {
        Exact start(minus_b_div_a, -1, d);
        Exact end(minus_b_div_a, 1, d);
        if (is_negative(start)) start = Exact(0);
        if (end > Exact(1)) end = Exact(1);
        if (start <= Exact(1) && end >= Exact(0)) {
            I = std::make_pair((start == Exact(0)) ? Lambda<C>(0) : Lambda<C>(start),
                               (end == Exact(1) ? Lambda<C>(1) : Lambda<C>(end)));
            return true;
        }
    }
    return false;
}



namespace HLPred
{

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * computes
*/
template <typename C>
Interval<C> intersection_interval(C const& curve1,
                               typename C::PointID const& center_id,
                               C const& curve2,
                               typename C::PointID seg_start_id,
                               typename C::distance_t const& radius)
{
    Interval<C> I;

    try {
      std::pair<Lambda<C>, Lambda<C>> II;
        // if not empty
        if (approximate_reals(curve1, center_id, curve2,  seg_start_id, radius, II)) {
            I = Interval<C>(II.first, II.second);
        }
    } catch (const CGAL::Uncertain_conversion_exception& e) {
      std::pair<Lambda<C>, Lambda<C>> II;
        // if not empty
        if (exact_reals(curve1, center_id, curve2,  seg_start_id, radius, II)) {
            I = Interval<C>(II.first, II.second);
        }
    }

    return I;
}

}  // namespace HLPred

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL
