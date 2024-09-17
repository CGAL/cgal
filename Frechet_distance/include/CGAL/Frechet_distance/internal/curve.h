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
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Kernel/Type_mapper.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Bbox.h>
#include <CGAL/STL_Extension/internal/Has_nested_type_Has_filtered_predicates_tag.h>
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
    using K = typename T::Kernel;
    using Self  = Curve<T>;

    static auto get_is_filtered()
    {
      if constexpr (::CGAL::internal::Has_nested_type_Has_filtered_predicates_tag<K>::value)
      {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            return std::true_type();
        }
        else
          return std::false_type();
      }
      else
        return std::false_type();
    }

    static constexpr bool is_filtered = decltype(get_is_filtered())::value;

    static auto get_type()
    {
      if constexpr (is_filtered)
      {
        return typename K::Exact_kernel{};
      }
      else
      {
        if constexpr (std::is_floating_point_v<typename K::FT>)
            return CGAL::Simple_cartesian<CGAL::Exact_rational>{};
        else
            return K{};
      }
    }

    using Rational_kernel = decltype(get_type());

public:
    static constexpr int dimension =  T::dimension;

    using distance_t = CGAL::Interval_nt<false>;
    using iKernel = std::conditional_t<(dimension == 2) || (dimension == 3),
                                       CGAL::Simple_cartesian<distance_t>,
                                       CGAL::Kernel_d_interface<Cartesian_base_d<distance_t,Dimension_tag<dimension>>>>;

    using Point = std::conditional_t<(dimension == 2),
                                     typename iKernel::Point_2 ,
                                     std::conditional_t<(dimension == 3),
                                     typename iKernel::Point_3,
                                     typename iKernel::Point_d>>;

    using Bbox = std::conditional_t<(dimension == 2),
                                     Bbox_2,
                                     std::conditional_t<(dimension == 3),
                                     Bbox_3,
                                     Bbox<Dimension_tag<dimension>,double>>>;

    using  Construct_bbox = std::conditional_t<(dimension == 2) || (dimension == 3),
                                                typename iKernel::Construct_bbox,
                                                typename iKernel::Construct_bbox_d>;

    using  Squared_distance = std::conditional_t<(dimension == 2) || (dimension == 3),
                                                typename iKernel::Compute_squared_distance,
                                                typename iKernel::Squared_distance_d>;
    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<typename T::Point>;

    using Rational = typename Rational_kernel::FT;
    using Rational_point = std::conditional_t<(dimension == 2),
                                     typename Rational_kernel::Point_2 ,
                                     std::conditional_t<(dimension == 3),
                                     typename Rational_kernel::Point_3,
                                     typename Rational_kernel::Point_d>>;
    using D2D = NT_converter<distance_t,double>;
    using I2R = std::conditional_t<(dimension == 2) || (dimension == 3),
                                    Cartesian_converter<iKernel, Rational_kernel, D2D>,
                                    KernelD_converter<iKernel, Rational_kernel>>;


    struct K2I
    {
      Point operator()(const typename T::Point& p) const
      {
        if constexpr (dimension == 2)
          return Point(p[0], p[1]);
        else
          return Point(p[0], p[1], p[2]);
      }
    };

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
    Curve(const PointRange& point_range) : prefix_length(point_range.size())
    {
        if constexpr ( ! std::is_floating_point<typename K::FT>::type::value) {
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
            auto segment_distance = distance(points[i - 1], points[i]);
            prefix_length[i] = prefix_length[i - 1] + segment_distance;
            bb += cbb(points[i]);
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
        if constexpr (std::is_floating_point<typename K::FT>::type::value) {
            I2R convert;
            return convert(points[i]);
        }
        if constexpr (is_filtered) {
            return typename K::C2E()(input[i]);
        }
        if constexpr (! std::is_floating_point<typename K::FT>::type::value &&
                      ! is_filtered) {
            return input[i];
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
    if (prefix_length.size()) {
        auto segment_distance = distance(points.back(), point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    } else {
        prefix_length.push_back(0);
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
