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
#include <CGAL/Frechet_distance/internal/id.h>

#include <CGAL/Interval_nt.h>
#include <CGAL/Kernel/Type_mapper.h>
#include <CGAL/Bbox.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <vector>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

double length_of_diagonal(const Bbox_2& bb)
{
    using I = Interval_nt<false>;
    I d = square(I(bb.xmax()) - I(bb.xmin()));
    d +=  square(I(bb.ymax()) - I(bb.ymin()));
    return sqrt(d).sup();
}

double length_of_diagonal(const Bbox_3& bb)
{
    using I = Interval_nt<false>;
    I d = square(I(bb.xmax()) - I(bb.xmin()));
    d +=  square(I(bb.ymax()) - I(bb.ymin()));
    d +=  square(I(bb.zmax()) - I(bb.zmin()));
    return sqrt(d).sup();
}

template<unsigned int N, typename T>
double length_of_diagonal(const Bbox<Dimension_tag<N>,T>& bb)
{
    using I = Interval_nt<false>;
    I d = square(I((bb.max)(0)) - I((bb.min)(0)));
    for(int i = 1; i < bb.dimension(); ++i){
      d += square(I((bb.max)(i)) - I((bb.min)(i)));
    }
    return sqrt(d).sup();
}

// TODO: Clean up Curve



/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a trajectory. Additionally to the points given in the
 * input file, we also store the length of any prefix of the trajectory.
 */
template <typename T>
class Curve
{
public:
    using Traits = T;
    using Self = Curve<T>;
    using K    = typename Traits::Kernel;

    static constexpr bool is_filtered = Traits::is_filtered;
    static constexpr bool is_floating_point = Traits::is_floating_point;
public:
    static constexpr int dimension =  T::dimension;

    using distance_t = typename Traits::distance_t;

    using iKernel = typename Traits::Approximate_kernel;

    using Point = typename Traits::Approximate_point;

    using Bbox = typename Traits::Bbox;

    using Construct_bbox = typename Traits::Construct_bbox;

    using Squared_distance = typename Traits::Squared_distance;

    using Difference_of_points = typename Traits::Difference_of_points;

    using Scaled_vector = typename Traits::Scaled_vector;

    using Translated_point = typename Traits::Translated_point;

    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<typename T::Point>;

    using Rational_kernel = typename Traits::Exact_kernel;
    using Rational = typename Rational_kernel::FT;
    using Rational_point = typename Traits::Exact_point;

    using I2R = typename Traits::A2E;

    using K2I = typename Traits::K2A;

    Curve() = default;

/*
      K           Has_filtered_predicates     is_floating_point(K::FT)      store input
    Epick                 Y                             Y                   N as we can  do to_double of interval
    Epeck                 Y                             N                   Y as we later need exact
    SC<double>            N                             Y                   N as we can  do to_double of interval
    SC<Rational>          N                             N                   Y as we later need exact
    SC<Interval_nt>       N                             N                   Y  to simplify life
    SC<Expr>              N                             N                   Y as we later need exact
    == Epeck_with_sqrt
*/
    template <class PointRange>
    Curve(const PointRange& point_range)
    : prefix_length(point_range.size())
    {
        if constexpr ( ! is_floating_point) {
            input.reserve(point_range.size());
            for (auto const& p : point_range) {
                input.push_back(p);
            }
        }

        points.reserve(point_range.size());
        if constexpr (is_filtered){
             for (auto const& p : point_range) {
                 points.push_back(typename K::C2F()(p));
            }
        }else{
            K2I convert;
             for (auto const& p : point_range) {
                 points.push_back(convert(p));
            }
        }

        if (points.empty()) {
            return;
        }

        prefix_length[0] = 0;

        Construct_bbox cbb;
        Bbox bb;
        for (PointID i = 1; i < points.size(); ++i) {
            auto segment_distance = distance(point(i - 1), point(i));
            prefix_length[i] = prefix_length[i - 1] + segment_distance;
            bb += cbb(point(i));
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
        if constexpr (is_floating_point) {
            I2R convert;
            return convert(point(i));
        }else{
          if constexpr (is_filtered) {
              return typename K::C2E()(input[i]);
          }else{
              return input[i];
          }
        }
        CGAL_assertion(false);
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
        return CGAL::sqrt(Squared_distance()(p, q));
    }


    template <class P>
    Point interpolate_at(P const& pt) const
    {
//        assert(pt.getFraction() >= Lambda<Self>(0) &&
//               pt.getFraction() <= Lambda<Self>(1));
        assert((
            pt.getPoint() < points.size() - 1 ||
            (pt.getPoint() == points.size() - 1 &&
//            is_zero(pt.getFraction()))));
            pt.getFraction()==0)));

//        if (is_zero(pt.getFraction())) {
        if (pt.getFraction()==0) {
            return point(pt.getPoint());
        }
        Difference_of_points difference;
        Scaled_vector scale;
        Translated_point translate;
        auto fraction = pt.getFraction().approx;
        return translate(point(pt.getPoint()) ,
                         scale(difference(point(pt.getPoint() ), point(pt.getPoint() + 1)), fraction));
    }

    Point interpolate_at(PointID const& pt) const { return point(pt); }

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

    Bbox const& bbox() const { return extreme_points; }

    distance_t getUpperBoundDistance(Curve const& other) const;

private:
    Points points;
    InputPoints input;
    std::vector<distance_t> prefix_length;
    Bbox extreme_points;
};

template <typename K>
std::ostream& operator<<(std::ostream& out, const Curve<K>& curve);

template <typename K>
void Curve<K>::push_back(Point const& point)
{
    if (prefix_length.empty()) {
        prefix_length.push_back(0);
    } else {
        auto segment_distance = distance(points.back(), point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }
    if (points.empty()){
      extreme_points = point.bbox();
    } else {
      extreme_points += point.bbox();
    }
    points.push_back(point);
}

template<typename K>
typename Curve<K>::distance_t Curve<K>::getUpperBoundDistance(Curve const& other) const
{
    Bbox bb = this->bbox() + other.bbox();
    return length_of_diagonal(bb);
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

} } } // namespace CGAL::Frechet_distance_::internal
