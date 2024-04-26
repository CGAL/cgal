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

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <optional>

#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Polyline_distance/internal/id.h>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

//
// distance_t
//
  using distance_t = CGAL::Interval_nt<false>;

//
// Point
//
using Point = CGAL::Simple_cartesian<distance_t>::Point_2;

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a value in the interval `[0,1]`....
*/
struct Lambda {
    typedef CGAL::Exact_rational Rational;
    typedef CGAL::Interval_nt<false> Approx;
    typedef CGAL::Sqrt_extension<Rational, Rational, CGAL::Tag_true,
                                 CGAL::Tag_false>
        Exact;

    mutable Approx approx;
    mutable std::optional<Exact> exact;
    Point circle_center;
    Point line_start, line_end;
    distance_t radius;
    mutable bool is_zero, is_one, is_exact, is_start;

    Lambda() {}

    Lambda(int zero_one)
        : approx(zero_one),
          is_zero(zero_one == 0),
          is_one(zero_one == 1),
          is_exact(true)
    {
    }

    Lambda(Approx approx, const Point& circle_center, distance_t radius,
           const Point& line_start, const Point& line_end, bool is_start)
        : approx(approx),
          circle_center(circle_center),
          line_start(line_start),
          line_end(line_end),
          radius(radius),
          is_zero(false),
          is_one(false),
          is_exact(false),
          is_start(is_start)
    {
    }

    Lambda(const Exact& exact)
        : approx(CGAL::to_interval(exact)),
          exact(exact),
          is_zero(CGAL::is_zero(exact)),
          is_one(CGAL::is_one(exact)),
          is_exact(true)
    {}


    bool update_exact() const
    {
        if (is_exact) {
            return true;
        }

        Rational a, b, c;
        for (auto i = 0; i < 2; ++i) {
            assert(line_start[i].is_point());
            assert(line_end[i].is_point());
            assert(circle_center[i].is_point());
            Rational start_end_diff = Rational(line_end[i].inf()) - Rational(line_start[i].inf());
            a += CGAL::square(start_end_diff);
            Rational start_center_diff = Rational(line_start[i].inf()) - Rational(circle_center[i].inf());
            b -= start_center_diff * start_end_diff;
            c += CGAL::square(start_center_diff);
        }
        c -= CGAL::square(Rational(to_double(radius)));

        Rational minus_b_div_a = b / a;
        Rational d = CGAL::square(minus_b_div_a) - c / a;
        if (! is_negative(d)) {
            if (is_start) {
                Exact start(minus_b_div_a, -1, d);
                if (is_negative(start)) start = Exact(0);
                if (start <= Exact(1)) {
                    exact = std::make_optional(start);
                    approx = CGAL::to_interval(*exact);
                    is_exact = true;
                    return true;
                }
            } else {
                Exact end(minus_b_div_a, 1, d);
                if (end > Exact(1)) end = Exact(1);
                if (! is_negative(end)) {
                    exact = std::make_optional(end);
                    approx = CGAL::to_interval(*exact);
                    is_exact = true;
                    return true;
                }
            }
        }
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, const Lambda&)
    {
        return os;
    }

    bool operator<=(const Lambda& other) const
    {
        return (*this < other) || (*this == other);
    }

    bool operator>=(const Lambda& other) const
    {
        return (other < *this) || (*this == other);
    }

    bool operator>(const Lambda& other) const
    {
      return (other < *this);
    }

    bool operator==(const Lambda& other) const
    {
        return (!(*this < other)) && (!(other < *this));
    }

    bool operator!=(const Lambda& other) const
    {
      return !(*this == other);
    }

    bool operator<(const Lambda& other) const
    {
        if ((is_zero && other.is_zero) || (is_one && other.is_one))
            return false;
        if ((is_zero && (!other.is_zero)) || (!is_one && other.is_one))
            return true;
        CGAL::Uncertain<bool> res = approx < other.approx;
        if (CGAL::is_certain(res)) {
            return CGAL::make_certain(res);
        }
        update_exact();
        other.update_exact();
        return exact < other.exact;
    }
};

} // namespace Polyline_distance
} // namespace internal


inline
bool is_one(const Polyline_distance::internal::Lambda& lambda)
{
    if (lambda.is_one) return true;
    if (lambda.is_zero) return false;
    Uncertain<bool> ub = is_one(lambda.approx);
    if (is_certain(ub)) {
        return make_certain(ub);
    }
    lambda.update_exact();
    return is_one(*lambda.exact);
}

inline
bool is_zero(const Polyline_distance::internal::Lambda& lambda)
{
    if (lambda.is_zero) return true;
    if (lambda.is_one) return false;
    Uncertain<bool> ub = is_zero(lambda.approx);
    if (is_certain(ub)) {
        return make_certain(ub);
    }
    lambda.update_exact();
    return is_zero(*lambda.exact);
}

inline CGAL::Interval_nt<false> to_interval(const Polyline_distance::internal::Lambda& lambda)
{
  return lambda.approx;
}



namespace Polyline_distance
{
namespace internal {


using RealType = Lambda;

inline
bool approximate_reals(const Point& circle_center, distance_t radius,
                       const Point& line_start, const Point& line_end,
                       std::pair<Lambda, Lambda>& I)
{

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
                (start == Approx(0)) ? Lambda(0)
                                     : Lambda(start, circle_center, radius,
                                              line_start, line_end, true),
                (end == Approx(1)) ? Lambda(1)
                                   : Lambda(end, circle_center, radius,
                                            line_start, line_end, false));
            return true;
        }
    }
    return false;
}

inline
bool exact_reals(const Point& circle_center, distance_t radius,
                 const Point& line_start, const Point& line_end,
                 std::pair<Lambda, Lambda>& I)
{
    typedef typename Lambda::Exact Exact;
    typedef typename Lambda::Rational Rational;

    Rational a(0), b(0), c(0);
    for (auto i = 0; i < 2; ++i) {
        assert(line_start[i].is_point());
        assert(line_end[i].is_point());
        assert(circle_center[i].is_point());
        Rational start_end_diff = Rational(line_end[i].inf()) - Rational(line_start[i].inf());
        a += CGAL::square(start_end_diff);
        Rational start_center_diff = Rational(line_start[i].inf()) - Rational(circle_center[i].inf());
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
            I = std::make_pair((start == Exact(0)) ? Lambda(0) : Lambda(start),
                               (end == Exact(1) ? Lambda(1) : Lambda(end)));
            return true;
        }
    }
    return false;
}



using PointID = ID<Point>;


struct ContinuousPoint {
    PointID point;
    distance_t fraction;

    bool operator<(ContinuousPoint const& other) const
    {
        return point < other.point ||
               (point == other.point && fraction < other.fraction);
    }
    bool operator<(PointID point_id) const { return point < point_id; }
    bool operator>(PointID point_id) const { return point >= point_id; }

    bool valid() const { return point.valid(); }
    void invalidate() { point.invalidate(); }
};



// Orientation and Direction etc.

enum class Direction : uint8_t { Up = 0, Down = 1, Left = 2, Right = 3 };
const std::array<Direction, 4> Directions = {
    {Direction::Up, Direction::Down, Direction::Left, Direction::Right}};

enum class Orientation { Horizontal = 0, Vertical = 1 };
Orientation operator!(Orientation orientation);

// short for: backward-forward direction
enum class BFDirection { Backward = 0, Forward = 1 };
BFDirection operator!(BFDirection direction);

// some helper functions for Orientation and Direction
Orientation getOrientation(Direction direction);
Direction getForwardDirection(Orientation orientation);
Direction getBackwardDirection(Orientation orientation);
bool isForward(Direction direction);
bool isBackward(Direction direction);
std::array<Direction, 2> getDirections(Orientation orientation);
BFDirection toBFDirection(Direction direction);

//
// Interval
//

// TODO: does CGAL have any replacement for this or do we want our custom type
// here?
/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
struct Interval {
    RealType begin;
    RealType end;

  Interval()
    : begin(1.), end(0.)
  {}

  Interval(RealType const& begin, RealType const& end)
    : begin(begin), end(end)
  {}

  bool is_empty() const
  {
    return begin > end;
  }

  bool intersects(Interval const& other) const
  {
    if (is_empty() || other.is_empty()) {
      return false;
    }

    return (other.begin >= begin && other.begin <= end) ||
      (other.end >= begin && other.end <= end) ||
      (other.begin <= begin && other.end >= end);
  }
};


std::ostream& operator<<(std::ostream& out, const Interval& interval);

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * Data Types for FrechetLight:
 *
 * CPoint is integral part + rational [0,1] number given as sqrt_extension
 */
class CPoint
{
private:
    PointID point;
    RealType fraction;

    void normalize()
    {
        assert(fraction >= 0. && fraction <= 1.);
        if (CGAL::is_one(fraction)) {
            fraction = 0.;
            ++point;
        }
    }

public:
    CPoint(PointID point, RealType const& fraction)
        : point(point), fraction(fraction)
    {
        normalize();
    }

  CPoint() : point((std::numeric_limits<PointID::IDType>::max)()), fraction(0.)
    {
    }

    PointID getPoint() const { return point; }

    RealType const& getFraction() const { return fraction; }

    double getFractionLB() const
    {
        // returns a lower bound to the fraction
        return CGAL::to_interval(fraction).inf();
    }

    void setPoint(PointID point) { this->point = point; }

    void setFraction(RealType const& frac)
    {
        fraction = frac;
        normalize();
    }

    bool operator<(CPoint const& other) const
    {
        return point < other.point ||
               (point == other.point && fraction < other.fraction);
    }

    bool operator<=(CPoint const& other) const
    {
        return point < other.point ||
               (point == other.point && fraction <= other.fraction);
    }

    bool operator>(CPoint const& other) const
    {
        return point > other.point ||
               (point == other.point && fraction > other.fraction);
    }

    bool operator>=(CPoint const& other) const
    {
        return point > other.point ||
               (point == other.point && fraction >= other.fraction);
    }

    bool operator==(CPoint const& other) const
    {
        return point == other.point && fraction == other.fraction;
    }

    bool operator!=(CPoint const& other) const
    {
        return point != other.point || fraction != other.fraction;
    }

    bool operator<(PointID other) const { return point < other; }

    bool operator>(PointID other) const
    {
        return point > other || (point == other && fraction > 0.);
    }

    bool operator<=(PointID other) const
    {
        return point < other || (point == other && fraction == 0.);
    }

    bool operator>=(PointID other) const { return point >= other; }

    bool operator==(PointID other) const
    {
        return point == other && fraction == 0.;
    }

    bool operator!=(size_t other) const { return !(point == other); }

    CPoint ceil() const
    {
      return (! CGAL::is_zero(fraction)) ? CPoint(point + 1, 0.) : CPoint(point, 0.);
    }

    CPoint floor() const { return CPoint(point, 0.); }

    std::string to_string() const
    {
        // return std::to_string( (double) point + fraction);
        std::stringstream stream;
        // AF  avoid +    stream << std::fixed << std::setprecision(10) <<
        // (double) point + fraction;
        return stream.str();
    }

    friend std::ostream& operator<<(std::ostream& out, const CPoint& p);
};

struct CInterval;
using CIntervals = std::vector<CInterval>;
using CIntervalsID = ID<CIntervals>;
using CIntervalID = std::size_t;


using CPoints = std::vector<CPoint>;

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

using CurveID = std::size_t;
using CurveIDs = std::vector<CurveID>;

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
struct CInterval {
    CPoint begin;
    CPoint end;

    const CInterval* reach_parent = nullptr;
  CPoint fixed = CPoint((std::numeric_limits<PointID::IDType>::max)(), 0.);
    CurveID fixed_curve = -1;

    CPosition getLowerRightPos() const
    {
        if (fixed_curve == 0) {
            CPosition ret = {{fixed, begin}};
            return ret;
        } else {
            CPosition ret = {{end, fixed}};
            return ret;
        }
    }
    CPosition getUpperLeftPos() const
    {
        if (fixed_curve == 0) {
            CPosition ret = {{fixed, end}};
            return ret;
        } else {
            CPosition ret = {{begin, fixed}};
            return ret;
        }
    }

    CInterval(CPoint const& begin, CPoint const& end, CPoint const& fixed,
              CurveID fixed_curve)
        : begin(begin), end(end), fixed(fixed), fixed_curve(fixed_curve)
    {
    }

    CInterval()
      : begin((std::numeric_limits<PointID::IDType>::max)(), 0.),
          end(std::numeric_limits<PointID::IDType>::lowest(), 0.)
    {
    }

    CInterval(CInterval const& other) = default;

    CInterval(CPoint const& begin, CPoint const& end) : begin(begin), end(end)
    {
    }

    CInterval(PointID point1, RealType fraction1, PointID point2,
              RealType const& fraction2)
        : begin(point1, fraction1), end(point2, fraction2)
    {
    }

    CInterval(PointID begin, PointID end) : begin(begin, 0.), end(end, 0.) {}

    bool operator<(CInterval const& other) const
    {
        return begin < other.begin || (begin == other.begin && end < other.end);
    }

    bool is_empty() const { return end < begin; }
    void make_empty()
    {
      begin = {(std::numeric_limits<PointID::IDType>::max)(), 0};
        end = {std::numeric_limits<PointID::IDType>::lowest(), 0};
    }
    void clamp(CPoint const& min, CPoint const& max)
    {
      begin = (std::max)(min, begin);
      end = (std::min)(max, end);
    }
};

std::ostream& operator<<(std::ostream& out, const CInterval& interval);

// TODO: CGAL certainly has a replacement for this; do we want this though as
// it's (probably) only used for visualization? Ellipse
struct Ellipse {
    Point center;
    distance_t width;
    distance_t height;
    double alpha;

    void invalidate()
    {
        width = -1.;
        height = -1.;
    }
    bool is_valid() { return width >= 0 && height >= 0; }
};



std::ostream& operator<<(std::ostream& out, const Point& p)
{
    out << std::setprecision(15) << "(" << p.x() << ", " << p.y() << ")";

    return out;
}

std::ostream& operator<<(std::ostream& out, const CPoint& p)
{
    out << std::setprecision(15) << "(" << (size_t)p.point << " + "
        << p.fraction << ")";

    return out;
}


//
// Interval
//

std::ostream& operator<<(std::ostream& out, const Interval& interval)
{
    // out << std::setprecision (15)
    //      << "(" << interval.begin << ", " << interval.end << ")";

    return out;
}

std::ostream& operator<<(std::ostream& out, const CInterval& interval)
{
    out << std::setprecision(15) << "(" << interval.begin << ", "
        << interval.end << ")";

    return out;
}

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL
