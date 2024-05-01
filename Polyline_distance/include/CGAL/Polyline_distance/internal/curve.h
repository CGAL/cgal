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
#include <CGAL/Polyline_distance/internal/id.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>

#include <vector>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

// TODO: Clean up Curve

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a trajectory. Additionally to the points given in the
 * input file, we also store the length of any prefix of the trajectory.
 */
template <typename T>
class Curve
{
    using K = typename T::BaseTraits;
    static auto get_type()
    {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            return typename K::Exact_kernel{};
        } else {
            return CGAL::Simple_cartesian<CGAL::Exact_rational>{};
        }
    }

    using Rational_kernel = decltype(get_type());

public:
    using distance_t = typename T::distance_t;
    using Point = typename T::iPoint;
    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<typename T::Point>;

    using Rational = typename Rational_kernel::FT;
    using Rational_point = typename Rational_kernel::Point_2;

    Curve() = default;

    template <class PointRange>
    Curve(const PointRange& point_range) : prefix_length(point_range.size())
    {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            input.reserve(point_range.size());
            for (auto const& p : point_range) {
                input.push_back(p);
            }
        }

        points.reserve(point_range.size());
        for (auto const& p : point_range) {
            if constexpr (K::Has_filtered_predicates_tag::value) {
                points.push_back(typename K::C2F()(p));
            } else if constexpr (std::is_floating_point<
                                     typename K::FT>::type::value) {
                points.push_back(Point(p.x(), p.y()));
            } else {
                points.push_back(Point_2(CGAL::to_interval(p.x()),
                                         CGAL::to_interval(p.y())));
            }
        }
        if (points.empty()) {
            return;
        }

        auto const& front = points.front();
        prefix_length[0] = 0;

        Bbox_2 bb;
        for (PointID i = 1; i < points.size(); ++i) {
            auto segment_distance = distance(points[i - 1], points[i]);
            prefix_length[i] = prefix_length[i - 1] + segment_distance;
            bb += points[i].bbox();
        }
        extreme_points = bb;
    }

    //  Curve(const Points& points);

    void reserve(std::size_t n) { points.reserve(n); }

    std::size_t size() const { return points.size(); }

    bool empty() const { return points.empty(); }

    Point const& operator[](PointID const& i) const { return points[i]; }

    Point const& point(PointID const& i) const { return points[i]; }

    Rational_point rpoint(PointID const& i) const
    {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            return typename K::C2E()(input[i]);
        }
        if constexpr (std::is_floating_point<typename K::FT>::type::value) {
            return Rational_point(points[i].x().inf(), points[i].y().inf());
        }
        assert(false);
        return Rational_point();
    }

    bool operator==(Curve const& other) const
    {
        return std::equal(points.cbegin(), points.cend(), other.points.cbegin(),
                          other.points.cend());
    }

    bool operator!=(Curve const& other) const { return !(*this == other); }

    distance_t distance(Point const& p, Point const& q) const
    {
        return CGAL::sqrt(CGAL::squared_distance(p, q));
    }

    template <class P>
    Point interpolate_at(P const& pt) const
    {
        assert(pt.getFraction() >= Lambda<T>(0) &&
               pt.getFraction() <= Lambda<T>(1));
        assert((
            pt.getPoint() < points.size() - 1 ||
            (pt.getPoint() == points.size() - 1 && is_zero(pt.getFraction()))));

        if (is_zero(pt.getFraction())) {
            return points[pt.getPoint()];
        }
        auto fraction = pt.getFraction().approx;
        return points[pt.getPoint()] +
               (fraction * (points[pt.getPoint() + 1] - points[pt.getPoint()]));
    }

    Point interpolate_at(PointID const& pt) const { return points[pt]; }

    distance_t curve_length(PointID const& i, PointID const& j) const
    {
        return prefix_length[j] - prefix_length[i];
    }

    Point const& front() const { return points.front(); }

    Point const& back() const { return points.back(); }

    void push_back(Point const& point);

    typename Points::const_iterator begin() { return points.begin(); }

    typename Points::const_iterator end() { return points.end(); }

    typename Points::const_iterator begin() const { return points.cbegin(); }

    typename Points::const_iterator end() const { return points.cend(); }

    Bbox_2 const& bbox() const { return extreme_points; }

    distance_t getUpperBoundDistance(Curve const& other) const;

private:
    Points points;
    InputPoints input;
    std::vector<distance_t> prefix_length;
    Bbox_2 extreme_points;
};

template <typename K>
std::ostream& operator<<(std::ostream& out, const Curve<K>& curve);

template <typename K>
void Curve<K>::push_back(Point const& point)
{
    if (prefix_length.size()) {
        auto segment_distance = distance(points.back(), point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    } else {
        prefix_length.push_back(0);
    }
    if (points.empty()){
      extreme_points = p.bbox();
    } else {
      extreme_points += p.bbox();
    }
    points.push_back(point);
}


template<typename K>
typename Curve<K>::distance_t Curve<K>::getUpperBoundDistance(Curve const& other) const
{
    Bbox_2 bb = this->bbox() + other.bbox();
    Point min_point = {bb.xmin(), bb.ymin()};
    Point max_point = {bb.xmax(), bb.ymax()};

    return distance(min_point, max_point);
}

template<typename K>
std::ostream& operator<<(std::ostream& out, const Curve<K>& curve)
{
    out << "[";
    for (auto const& point : curve) {
        out << point << ", ";
    }
    out << "]";

    return out;
}

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL
