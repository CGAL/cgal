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

#include <CGAL/Frechet_distance/internal/id.h>
#include <CGAL/Frechet_distance/internal/curve.h>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {



/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a value in the interval `[0,1]`....
*/
template<typename C>
struct Lambda {
    using Curve = C;
    using distance_t = typename C::distance_t;
    using Approx = distance_t;
    using PointID = typename C::PointID;
    using Rational = typename Curve::Rational;
    using Exact = CGAL::Sqrt_extension<Rational, Rational, CGAL::Tag_true, CGAL::Tag_false>;

    mutable Approx approx;
    mutable std::optional<Exact> exact;
    const Curve * curve1;
    PointID circle_center;
    const Curve * curve2;
    PointID line_start;
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

  Lambda(const Approx& approx, const Curve& curve1, const PointID& circle_center,
         const Curve& curve2, const PointID& line_start, const distance_t& radius, bool is_start)
        : approx(approx),
          curve1(&curve1),
          circle_center(circle_center),
          curve2(&curve2),
          line_start(line_start),
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

        const typename Curve::Rational_point&  ls = curve2->rpoint(line_start);
        const typename Curve::Rational_point& le = curve2->rpoint(line_start+1);
        const typename Curve::Rational_point& cc = curve1->rpoint(circle_center);
        Rational a, b, c;
        for (auto i = 0; i < Curve::dimension; ++i) {
            Rational start_end_diff = le[i] - ls[i];
            a += CGAL::square(start_end_diff);
            Rational start_center_diff = ls[i] - cc[i];
            b -= start_center_diff * start_end_diff;
            c += CGAL::square(start_center_diff);
        }
        CGAL_assertion(radius.is_point());
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

    friend std::ostream& operator<<(std::ostream& os, const Lambda<C>&)
    {
        return os;
    }

    bool operator<=(const Lambda<C>& other) const
    {
        return (*this < other) || (*this == other);
    }

    bool operator>=(const Lambda<C>& other) const
    {
        return (other < *this) || (*this == other);
    }

    bool operator>(const Lambda<C>& other) const
    {
      return (other < *this);
    }

    bool operator==(const Lambda<C>& other) const
    {
        return (!(*this < other)) && (!(other < *this));
    }

    bool operator!=(const Lambda<C>& other) const
    {
      return !(*this == other);
    }

    bool operator<(const Lambda<C>& other) const
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

} } // namespace Frechet_distance_::internal

template <typename C>
bool is_one(const Frechet_distance_::internal::Lambda<C>& lambda)
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

template <typename C>
bool is_zero(const Frechet_distance_::internal::Lambda<C>& lambda)
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

/*
template <typename K>
CGAL::Interval_nt<false> to_interval(const Frechet_distance_::internal::Lambda<K>& lambda)
{
  return lambda.approx;
}
*/


namespace Frechet_distance_
{
namespace internal {

template <typename C>
struct ContinuousPoint {
    using distance_t = typename C::distance_t;
    using PointID = typename C::PointID;
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
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a
*/

template <typename C>
struct Interval {
    Lambda<C> begin;
    Lambda<C> end;

  Interval()
    : begin(1.), end(0.)
  {}

  Interval(Lambda<C> const& begin, Lambda<C> const& end)
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

template <typename C>
std::ostream& operator<<(std::ostream& out, const Interval<C>& interval);

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * Data Types for FrechetLight:
 *
 * CPoint is integral part + rational [0,1] number given as sqrt_extension
 */
template <typename C>
class CPoint
{
  private:
  using PointID = typename C::PointID;
    PointID point;
    Lambda<C> fraction;

    void normalize()
    {
        assert(fraction >= Lambda<C>(0) && fraction <= Lambda<C>(1));
        if (CGAL::is_one(fraction)) {
            fraction = 0.;
            ++point;
        }
    }

public:
    CPoint(PointID point, Lambda<C> const& fraction)
        : point(point), fraction(fraction)
    {
        normalize();
    }

  CPoint() : point((std::numeric_limits<typename PointID::IDType>::max)()), fraction(0.)
    {
    }

    PointID getPoint() const { return point; }

    Lambda<C> const& getFraction() const { return fraction; }

    double getFractionLB() const
    {
        // returns a lower bound to the fraction
        return CGAL::to_interval(fraction).inf();
    }

    void setPoint(PointID point) { this->point = point; }

    void setFraction(Lambda<C> const& frac)
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
        return point > other || (point == other && fraction > Lambda<C>(0));
    }

    bool operator<=(PointID other) const
    {
        return point < other || (point == other && is_zero(fraction));
    }

    bool operator>=(PointID other) const { return point >= other; }

    bool operator==(PointID other) const
    {
        return point == other &&  is_zero(fraction);
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

    friend std::ostream& operator<<(std::ostream& out, const CPoint<C>& p)
    {
      out << std::setprecision(15) << "(" << (size_t)p.point << " + "
          << p.fraction << ")";

      return out;
    }
};



template <typename C>
using CPosition = std::array<CPoint<C>, 2>;

template <typename C>
using CPositions = std::vector<CPosition<C>>;

using CurveID = std::size_t;

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a
*/
template <typename C>
struct CInterval {
    using PointID =  typename C::PointID;

    CPoint<C> begin;
    CPoint<C> end;

    const CInterval* reach_parent = nullptr;
  CPoint<C> fixed = CPoint<C>((std::numeric_limits<typename PointID::IDType>::max)(), 0.);
    CurveID fixed_curve = -1;

    CPosition<C> getLowerRightPos() const
    {
        if (fixed_curve == 0) {
            CPosition<C> ret = {{fixed, begin}};
            return ret;
        } else {
            CPosition<C> ret = {{end, fixed}};
            return ret;
        }
    }
    CPosition<C> getUpperLeftPos() const
    {
        if (fixed_curve == 0) {
            CPosition<C> ret = {{fixed, end}};
            return ret;
        } else {
            CPosition<C> ret = {{begin, fixed}};
            return ret;
        }
    }

    CInterval(CPoint<C> const& begin, CPoint<C> const& end, CPoint<C> const& fixed,
              CurveID fixed_curve)
        : begin(begin), end(end), fixed(fixed), fixed_curve(fixed_curve)
    {
    }

    CInterval()
      : begin((std::numeric_limits<typename PointID::IDType>::max)(), 0.),
          end(std::numeric_limits<typename PointID::IDType>::lowest(), 0.)
    {
    }

    CInterval(CInterval const& other) = default;

    CInterval(CPoint<C> const& begin, CPoint<C> const& end) : begin(begin), end(end)
    {
    }

    CInterval(PointID point1, Lambda<C> fraction1, PointID point2,
              Lambda<C> const& fraction2)
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
      begin = {(std::numeric_limits<typename PointID::IDType>::max)(), 0};
        end = {std::numeric_limits<typename PointID::IDType>::lowest(), 0};
    }

    void clamp(CPoint<C> const& min, CPoint<C> const& max)
    {
      begin = (std::max)(min, begin);
      end = (std::min)(max, end);
    }
};

template <typename C>
std::ostream& operator<<(std::ostream& out, const CInterval<C>& interval);


template <typename C>
using CIntervals = std::vector<CInterval<C>>;

template <typename C>
using CIntervalsID = ID<CIntervals<C>>;

#if 0
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
#endif

//
// Interval
//
template <typename C>
std::ostream& operator<<(std::ostream& out, const Interval<C>& interval)
{
    out << std::setprecision (15)
        << "(" << interval.begin << ", " << interval.end << ")";

    return out;
}
template <typename C>
std::ostream& operator<<(std::ostream& out, const CInterval<C>& interval)
{
    out << std::setprecision(15) << "(" << interval.begin << ", "
        << interval.end << ")";

    return out;
}

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL
