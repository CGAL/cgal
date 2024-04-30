// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#pragma once
#include <CGAL/license/Polyline_distance.h>

#include <CGAL/Polyline_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Polyline_distance {
namespace internal {



template <typename K>
bool approximate_reals(Curve<K> const& curve1,
                       PointID const& center_id,
                       Curve<K> const& curve2,
                       PointID const& seg_start_id,
                       distance_t const& radius,
                       std::pair<Lambda<K>, Lambda<K>>& I)
{
    const Point& circle_center = curve1[center_id];
    const Point& line_start = curve2[seg_start_id];
    const Point& line_end  = curve2[seg_start_id + 1];

    // let a := v.x()^2 + v.y()^2,
    // let b := line_start.x() * v.x() + line_start.y() * v.y(),
    // let c := line_start.x()^2 + line_start.y()^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)

    typedef CGAL::Interval_nt<false> Approx;
    Approx a(0), b(0), c(0);
    for (auto i = 0; i < 2; ++i) {
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
                (start == Approx(0)) ? Lambda<K>(0)
                                     : Lambda<K>(start, curve1, center_id,
                                              curve2, seg_start_id, radius, true),
                (end == Approx(1)) ? Lambda<K>(1)
                                   : Lambda<K>(end, curve1, center_id,
                                           curve2, seg_start_id, radius, false));
            return true;
        }
    }
    return false;
}

template <typename K>
bool exact_reals(Curve<K> const& curve1,
                 PointID const& center_id,
                 Curve<K> const& curve2,
                 PointID const& seg_start_id,
                 distance_t const& radius,
                 std::pair<Lambda<K>, Lambda<K>>& I)
{
  typedef typename Curve<K>::Rational_point Rational_point;
  const Rational_point& circle_center = curve1.rpoint(center_id);
  const Rational_point& line_start = curve2.rpoint(seg_start_id);
  const Rational_point& line_end  = curve2.rpoint(seg_start_id + 1);

    typedef typename Lambda<K>::Exact Exact;
    typedef typename Lambda<K>::Rational Rational;

    Rational a(0), b(0), c(0);
    for (auto i = 0; i < 2; ++i) {
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
            I = std::make_pair((start == Exact(0)) ? Lambda<K>(0) : Lambda<K>(start),
                               (end == Exact(1) ? Lambda<K>(1) : Lambda<K>(end)));
            return true;
        }
    }
    return false;
}



namespace HLPred
{

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * computes
*/
template <typename K>
Interval<K> intersection_interval(Curve<K> const& curve1,
                               PointID const& center_id,
                               Curve<K> const& curve2,
                               PointID seg_start_id,
                               distance_t const& radius)
{
    Interval<K> I;

    try {
      std::pair<RealType<K>, RealType<K>> II;
        // if not empty
        if (approximate_reals(curve1, center_id, curve2,  seg_start_id, radius, II)) {
            I = Interval<K>(II.first, II.second);
        }
    } catch (const CGAL::Uncertain_conversion_exception& e) {
      std::pair<RealType<K>, RealType<K>> II;
        // if not empty
        if (exact_reals(curve1, center_id, curve2,  seg_start_id, radius, II)) {
            I = Interval<K>(II.first, II.second);
        }
    }

    return I;
}

}  // namespace HLPred

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL
