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
#include <CGAL/Interval_nt.h>

#include <vector>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

// TODO: Clean up Curve


/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a trajectory. Additionally to the points given in the input file,
 * we also store the length of any prefix of the trajectory.
 */
class Curve
{
public:
  using distance_t = CGAL::Interval_nt<false>;
  using Point = CGAL::Simple_cartesian<distance_t>::Point_2;
  using PointID = ID<Point>;
  using Points = std::vector<Point>;

  Curve() = default;

  Curve(const Points& points);

  void reserve(std::size_t n)
  {
    points.reserve(n);
  }

    std::size_t size() const { return points.size(); }

  bool empty() const { return points.empty(); }

  Point const& operator[](PointID const& i) const { return points[i]; }

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
        assert(pt.getFraction() >= 0. && pt.getFraction() <= 1.);
        assert(
            (pt.getPoint() < points.size() - 1 ||
             (pt.getPoint() == points.size() - 1 && pt.getFraction() == 0.)));

        if(pt.getFraction() == 0){
          return points[pt.getPoint()];
        }
        auto fraction = to_interval(pt.getFraction());
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

    Points::const_iterator begin() { return points.begin(); }

    Points::const_iterator end() { return points.end(); }

    Points::const_iterator begin() const { return points.cbegin(); }

    Points::const_iterator end() const { return points.cend(); }

    struct ExtremePoints {
        distance_t min_x, min_y, max_x, max_y;
    };

    ExtremePoints const& getExtremePoints() const;

    distance_t getUpperBoundDistance(Curve const& other) const;

private:
    Points points;
    std::vector<distance_t> prefix_length;
    ExtremePoints extreme_points;

};

std::ostream& operator<<(std::ostream& out, const Curve& curve);

Curve::Curve(const Points& points)
    : points(points), prefix_length(points.size())
{
    if (points.empty()) {
        return;
    }

    auto const& front = points.front();
    extreme_points = {front.x(), front.y(), front.x(), front.y()};
    prefix_length[0] = 0;

    for (PointID i = 1; i < points.size(); ++i) {
        auto segment_distance = distance(points[i - 1], points[i]);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points.min_x = (std::min)(extreme_points.min_x, points[i].x());
        extreme_points.min_y = (std::min)(extreme_points.min_y, points[i].y());
        extreme_points.max_x = (std::max)(extreme_points.max_x, points[i].x());
        extreme_points.max_y = (std::max)(extreme_points.max_y, points[i].y());
    }
}

void Curve::push_back(Point const& point)
{
    if (prefix_length.size()) {
        auto segment_distance = distance(points.back(), point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    } else {
        prefix_length.push_back(0);
    }
    if (points.empty()) {
        extreme_points = {point.x(), point.y(), point.x(), point.y()};
    } else {
        extreme_points.min_x = (std::min)(extreme_points.min_x, point.x());
        extreme_points.min_y = (std::min)(extreme_points.min_y, point.y());
        extreme_points.max_x = (std::max)(extreme_points.max_x, point.x());
        extreme_points.max_y = (std::max)(extreme_points.max_y, point.y());
    }
    points.push_back(point);
}

auto Curve::getExtremePoints() const -> ExtremePoints const&
{
    assert(!points.empty());
    return extreme_points;
}

Curve::distance_t Curve::getUpperBoundDistance(Curve const& other) const
{
    auto const& extreme1 = this->getExtremePoints();
    auto const& extreme2 = other.getExtremePoints();
    Point min_point{(std::min)(extreme1.min_x, extreme2.min_x),
                    (std::min)(extreme1.min_y, extreme2.min_y)};
    Point max_point = {(std::max)(extreme1.max_x, extreme2.max_x),
                       (std::max)(extreme1.max_y, extreme2.max_y)};

    return distance(min_point, max_point);
}

std::ostream& operator<<(std::ostream& out, const Curve& curve)
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
