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

#ifndef CGAL_FRECHET_DISTANCE_INTERNAL_GEOMETRY_BASICS_H
#define CGAL_FRECHET_DISTANCE_INTERNAL_GEOMETRY_BASICS_H

#include <CGAL/license/Frechet_distance.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <optional>

#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Frechet_distance/internal/Id.h>
#include <CGAL/Frechet_distance/internal/Curve.h>

namespace CGAL {
namespace Frechet_distance {
namespace internal {



/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a value in the interval `[0,1]`....
*/
template<typename C>
struct Lambda;

// filtered version
template<typename FilteredTraits>
struct Lambda<Curve<FilteredTraits,true>>
{
    using Curve = ::CGAL::Frechet_distance::internal::Curve<FilteredTraits,true>;
    using C = Curve;
    using distance_t = typename C::distance_t;
    using Approx = distance_t;
    using PointID = typename C::PointID;
    using Rational = typename Curve::Rational;
    using Exact = typename Root_of_traits<Rational>::Root_of_2;

    mutable Approx approx;
    mutable std::optional<Exact> exact;
    const Curve * curve1 = nullptr;
    PointID circle_center = -1;
    const Curve * curve2 = nullptr;
    PointID line_start=-1;
    distance_t radius=0;
    mutable bool is_zero, is_one, is_exact, is_start;

    Lambda(int zero_one)
        : approx(zero_one),
          is_zero(zero_one == 0),
          is_one(zero_one == 1),
          is_exact(true),
          is_start(false)
    {
    }

    Lambda() : Lambda(0) {}

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

    Lambda(const Exact& exact, const Curve&, const PointID&,
                               const Curve&, const PointID&,
                               const distance_t&, bool)
        : approx(CGAL::to_interval(exact)),
          exact(exact),
          is_zero(CGAL::is_zero(exact)),
          is_one(CGAL::is_one(exact)),
          is_exact(true)
    {}

    //     fill_lambda returns a pair and we are only interested in a bound
    bool update_exact() const
    {
        if (is_exact){
            if (!exact.has_value()){
                exact = (is_one) ? std::make_optional(Exact(1)) : std::make_optional(Exact(0));
            }
            return true;
        }

        const typename Curve::Rational_point& ls = curve2->rpoint(line_start);
        const typename Curve::Rational_point& le = curve2->rpoint(line_start+1);
        const typename Curve::Rational_point& cc = curve1->rpoint(circle_center);

        auto ccci = curve1->rational_traits().construct_cartesian_const_iterator_d_object();
        Rational a, b, c;

        auto it_le=ccci(le), it_ls=ccci(ls), it_cc=ccci(cc);

        for (auto i = 0; i < Curve::dimension; ++i) {
            Rational start_end_diff = *it_le - *it_ls;
            a += CGAL::square(start_end_diff);
            Rational start_center_diff = *it_ls - *it_cc;
            b -= start_center_diff * start_end_diff;
            c += CGAL::square(start_center_diff);
            ++it_le; ++it_ls; ++it_cc;
        }
        CGAL_assertion(radius.is_point());
        c -= CGAL::square(Rational(to_double(radius)));

        CGAL_assertion(a!=0);
        Rational minus_b_div_a = b / a;
        Rational d = CGAL::square(minus_b_div_a) - c / a;
        if (! is_negative(d)) {
            if (is_start) {
                Exact start = make_root_of_2(minus_b_div_a, -1, d);
                if (is_negative(start)) start = Exact(0);
                if (start <= Exact(1)) {
                    exact = std::make_optional(start);
                    approx = CGAL::to_interval(*exact);
                    is_exact = true;
                    return true;
                }
            } else {
                Exact end = make_root_of_2(minus_b_div_a, 1, d);
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
        // AF this may be wrong if approx is [0.9,1.1]  and other.is_one:
        // if ((is_zero && (!other.is_zero)) || (!is_one && other.is_one))
        if(is_zero && other.is_one)
            return true;
        CGAL::Uncertain<bool> res = approx < other.approx;
        if (CGAL::is_certain(res)) {
            return CGAL::make_certain(res);
        }
        update_exact();
        other.update_exact();
        bool eres = exact.value() < other.exact.value();
        return eres;
    }
};

// non-filtered version
template<typename T>
struct Lambda<Curve<T,false>>
{
    using Curve = ::CGAL::Frechet_distance::internal::Curve<T,false>;
    using C = Curve;
    using FT = typename C::FT;
    using RO2 = typename Root_of_traits<FT>::Root_of_2;
    using PointID = typename C::PointID;

    RO2 value;
    FT approx; // TODO: isn't it an issue to use that with an exact FT?  Ask @SL
               // TODO: we could use the interval in the RO2?            Ask @SL

    bool is_zero, is_one;

    Lambda() {}

    Lambda(int zero_one)
      : value(zero_one)
      , approx(zero_one)
      , is_zero(zero_one == 0)
      , is_one(zero_one == 1)
    {}

    Lambda(const RO2& v, const Curve& /* curve1 */, const PointID& /* circle_center */,
         const Curve& /* curve2 */, const PointID& /* line_start */, const FT& /* radius */, bool /* is_start */)
      : value(v)
      , is_zero(false)
      , is_one(false)
    {
      if constexpr (!std::is_same_v<RO2, FT>)
      {
        std::pair<double, double> iv = to_interval(v);
        approx = FT((iv.first+iv.second)/2.);
      }
      else
        approx=v;
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
        return value < other.value;
    }
};

} } // namespace Frechet_distance::internal

template <typename FilteredTraits>
bool is_one(const Frechet_distance::internal::Lambda<Frechet_distance::internal::Curve<FilteredTraits,true>>& lambda)
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

template <typename FilteredTraits>
bool is_zero(const Frechet_distance::internal::Lambda<Frechet_distance::internal::Curve<FilteredTraits,true>>& lambda)
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

template <typename T>
bool is_one(const Frechet_distance::internal::Lambda<Frechet_distance::internal::Curve<T,false>>& lambda)
{
    if (lambda.is_one) return true;
    if (lambda.is_zero) return false;
    return is_one(lambda.value);
}

template <typename T>
bool is_zero(const Frechet_distance::internal::Lambda<Frechet_distance::internal::Curve<T,false>>& lambda)
{
    if (lambda.is_zero) return true;
    if (lambda.is_one) return false;
    return is_zero(lambda.value);
}

/*
template <typename K>
CGAL::Interval_nt<false> to_interval(const Frechet_distance::internal::Lambda<K>& lambda)
{
  return lambda.approx;
}
*/


namespace Frechet_distance
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

// TODO: CGAL has no Interval<T> class template yet, but it would look like this probably

template <typename C>
struct Interval {
    Lambda<C> begin;
    Lambda<C> end;

  Interval()
    : begin(1), end(0)
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
            fraction = Lambda<C>(0);
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
  CPoint<C> fixed = CPoint<C>((std::numeric_limits<typename PointID::IDType>::max)(), 0);
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
      : begin((std::numeric_limits<typename PointID::IDType>::max)(), 0),
          end(std::numeric_limits<typename PointID::IDType>::lowest(), 0)
    {
    }

    CInterval(CPoint<C> const& begin, CPoint<C> const& end) : begin(begin), end(end)
    {
    }

    CInterval(PointID point1, Lambda<C> fraction1, PointID point2,
              Lambda<C> const& fraction2)
        : begin(point1, fraction1), end(point2, fraction2)
    {
    }

    CInterval(PointID begin, PointID end) : begin(begin, 0), end(end, 0) {}

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
} // namespace Frechet_distance
} // namespace CGAL

#endif // CGAL_FRECHET_DISTANCE_INTERNAL_GEOMETRY_BASICS_H
