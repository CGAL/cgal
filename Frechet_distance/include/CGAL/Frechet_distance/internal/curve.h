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
#include <CGAL/Cartesian_converter.h>
#include <vector>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

//TODO: move that in Kernel_23/Kernel_d?
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

// TODO: Clean up Curve + factorize filtered vs non-filtered version with a base class

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a trajectory. Additionally to the points given in the
 * input file, we also store the length of any prefix of the trajectory.
 */
template <typename Traits, bool is_filtered=false>
class Curve;

//filtered version
template <typename Approximate_traits, typename Rational_traits>
class Curve<std::pair<Approximate_traits,Rational_traits>, true>
{
public:
  using distance_t = typename Approximate_traits::FT;
  using Point = typename Approximate_traits::Point;
  using Squared_distance = typename Approximate_traits::Squared_distance;

//TODO: handle that
//    static constexpr bool is_floating_point = Traits::is_floating_point;
//    using Construct_cartesian_const_iterator_d = typename Traits::Construct_cartesian_const_iterator_d;

  //TODO: we expect Dimension_tag for now.
  static constexpr int dimension =  Approximate_traits::Dimension::value;
  using Bbox = std::conditional_t<dimension==2,
                                  Bbox_2, std::conditional_t<dimension==3,
                                                             Bbox_3, ::CGAL::Bbox<typename Approximate_traits::Dimension, double>>>;
//////

    using PointID = ID<Point>;
    using Points = std::vector<Point>;
//    using InputPoints = std::vector<typename T::Point>;

    using Rational = typename Rational_traits::FT;
    using Rational_point = typename Rational_traits::Point;

    // TODO: this assumes that input interval are all tight --> need a PM!
    using I2R = Cartesian_converter< typename Kernel_traits<Point>::Kernel,
                                     typename Kernel_traits<Rational_point>::Kernel, NT_converter<distance_t,double>>;

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
        using IPoint = typename std::iterator_traits<typename PointRange::const_iterator>::value_type;
        using K2I = Cartesian_converter< typename Kernel_traits<IPoint>::Kernel,
                                         typename Kernel_traits<Point>::Kernel>;

        //TODO: not working for now with EPECK or SC<Rational>
        //~ if constexpr ( ! is_floating_point) {
            //~ input.reserve(point_range.size());
            //~ for (auto const& p : point_range) {
                //~ input.push_back(p);
            //~ }
        //~ }

        points.reserve(point_range.size());
        K2I convert;
        for (auto const& p : point_range)
          points.push_back(convert(p));

        if (points.empty()) return;

        prefix_length[0] = 0;

        std::array<double, dimension> min_c, max_c;
        for (int d=0; d<dimension; ++d)
        {
          min_c[d] = point(0)[d].inf();
          max_c[d] = point(0)[d].sup();
        }

        for (PointID i = 1; i < points.size(); ++i) {
            auto segment_distance = distance(point(i - 1), point(i));
            prefix_length[i] = prefix_length[i - 1] + segment_distance;
            for (int d=0; d<dimension; ++d)
            {
              min_c[d] = (std::min)(min_c[d], point(i)[d].inf());
              max_c[d] = (std::max)(max_c[d], point(i)[d].sup());
            }
        }
        if constexpr (dimension==2)
          extreme_points=Bbox_2(min_c[0], min_c[1], max_c[0], max_c[1]);
        else if constexpr (dimension==3)
          extreme_points=Bbox_3(min_c[0], min_c[1], min_c[2], max_c[0], max_c[1], max_c[2]);
        else
        {
          for (int d=0; d<dimension; ++d)
          {
            extreme_points.min(d)=min_c[d];
            extreme_points.max(d)=max_c[d];
          }
        }
    }

    //  Curve(const Points& points);

    void reserve(std::size_t n) { points.reserve(n); }

    std::size_t size() const { return points.size(); }

    bool empty() const { return points.empty(); }

    Point const& operator[](PointID const& i) const { return points[i]; }

    Point const& point(PointID const& i) const { return points[i]; }

    Rational_point rpoint(PointID const& i) const
    {
        //~ if constexpr (is_floating_point) {
            I2R convert;
            return convert(point(i));
        //~ }else{
          //~ if constexpr (is_filtered) {
              //~ return typename K::C2E()(input[i]);
          //~ }else{
              //~ return input[i];
          //~ }
        //~ }
        //~ CGAL_assertion(false);
        //~ return Rational_point();
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

        std::array<distance_t, dimension> b_coords;
        auto fraction = pt.getFraction().approx;
        distance_t one_m_f = distance_t(1) - fraction;
        for (int d=0; d<dimension; ++d)
          //TODO: use Construct_cartesian_const_iterator_d
          //TODO: use contruct point and add it in the traits concept
          b_coords[d] = one_m_f * point(pt.getPoint())[d] +  point(pt.getPoint() + 1)[d] * fraction;

        if constexpr (dimension==2)
          return Point(b_coords[0], b_coords[1]);
        else if constexpr (dimension==3)
          return Point(b_coords[0], b_coords[1], b_coords[2]);
        else
          return Point(b_coords.begin(), b_coords.end());
    }

    Point interpolate_at(PointID const& pt) const { return point(pt); }

    distance_t curve_length(PointID const& i, PointID const& j) const
    {
        return prefix_length[j] - prefix_length[i];
    }

    Point const& front() const { return points.front(); }

    Point const& back() const { return points.back(); }

    void push_back(Point const& point)
    {
      if (prefix_length.empty()) {
          prefix_length.push_back(0);
      } else {
          auto segment_distance = distance(points.back(), point);
          prefix_length.push_back(prefix_length.back() + segment_distance);
      }
  //TODO update that code
      if (points.empty()){
        extreme_points = point.bbox();
      } else {
        extreme_points += point.bbox();
      }
      points.push_back(point);
    }

    typename Points::const_iterator begin() { return points.begin(); }

    typename Points::const_iterator end() { return points.end(); }

    typename Points::const_iterator begin() const { return points.cbegin(); }

    typename Points::const_iterator end() const { return points.cend(); }

    Bbox const& bbox() const { return extreme_points; }

    distance_t getUpperBoundDistance(Curve const& other) const
    {
      Bbox bb = this->bbox() + other.bbox();
      return length_of_diagonal(bb).sup();
    }

private:
    Points points;
//    InputPoints input;
    std::vector<distance_t> prefix_length;
    Bbox extreme_points;
};

template <typename T>
class Curve<T, false>
{
public:
    using Traits = T;
//    using K    = typename Traits::Kernel;

//   using Construct_cartesian_const_iterator_d = typename Traits::Construct_cartesian_const_iterator_d;
  //TODO: we expect Dimension_tag for now.
  static constexpr int dimension =  Traits::Dimension::value;
  using Bbox = std::conditional_t<dimension==2,
                                  Bbox_2, std::conditional_t<dimension==3,
                                                             Bbox_3, ::CGAL::Bbox<typename Traits::Dimension, double>>>;


//////
    using FT = typename Traits::FT;
    using distance_t = FT;

    using Point = typename Traits::Point;

    using Squared_distance = typename Traits::Squared_distance;

    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<typename T::Point>;


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
    , points(point_range.begin(), point_range.end())
    {
        if (points.empty()) {
            return;
        }

        prefix_length[0] = 0;

        std::array<double, dimension> min_c, max_c;
        for (int d=0; d<dimension; ++d)
        {
          //TODO: use already know intervals
          auto inter = to_interval(point(0)[d]);
          min_c[d] = inter.first;
          max_c[d] = inter.second;
        }

        for (PointID i = 1; i < points.size(); ++i) {
            auto segment_distance = distance(point(i - 1), point(i));
            prefix_length[i] = prefix_length[i - 1] + segment_distance;
            for (int d=0; d<dimension; ++d)
            {
              //TOOD: use already know intervals
              auto inter = to_interval(point(i)[d]);
              min_c[d] = (std::min)(min_c[d], inter.first);
              max_c[d] = (std::max)(max_c[d], inter.second);
            }
        }
        if constexpr (dimension==2)
          extreme_points=Bbox_2(min_c[0], min_c[1], max_c[0], max_c[1]);
        else if constexpr (dimension==3)
          extreme_points=Bbox_3(min_c[0], min_c[1], min_c[2], max_c[0], max_c[1], max_c[2]);
        else
        {
          for (int d=0; d<dimension; ++d)
          {
            extreme_points.min(d)=min_c[d];
            extreme_points.max(d)=max_c[d];
          }
        }
    }

    //  Curve(const Points& points);

    void reserve(std::size_t n) { points.reserve(n); }

    std::size_t size() const { return points.size(); }

    bool empty() const { return points.empty(); }

    Point const& operator[](PointID const& i) const { return points[i]; }

    Point const& point(PointID const& i) const { return points[i]; }

    bool operator==(Curve const& other) const
    {
        return std::equal(points.cbegin(), points.cend(), other.points.cbegin(),
                          other.points.cend());
    }

    bool operator!=(Curve const& other) const { return !(*this == other); }

    distance_t distance(Point const& p, Point const& q) const
    {
      // TODO: is that an issue for robustness?
        return CGAL::approximate_sqrt(Squared_distance()(p, q));
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

        std::array<distance_t, dimension> b_coords;
        auto fraction = pt.getFraction().approx;
        distance_t one_m_f = distance_t(1) - fraction;
        for (int d=0; d<dimension; ++d)
          //TODO: use Construct_cartesian_const_iterator_d
          //TODO: use contruct point and add it in the traits concept
          b_coords[d] = one_m_f * point(pt.getPoint())[d] +  point(pt.getPoint() + 1)[d] * fraction;

        if constexpr (dimension==2)
          return Point(b_coords[0], b_coords[1]);
        else if constexpr (dimension==3)
          return Point(b_coords[0], b_coords[1], b_coords[2]);
        else
          return Point(b_coords.begin(), b_coords.end());
    }

    Point interpolate_at(PointID const& pt) const { return point(pt); }

    distance_t curve_length(PointID const& i, PointID const& j) const
    {
        return prefix_length[j] - prefix_length[i];
    }

    Point const& front() const { return points.front(); }

    Point const& back() const { return points.back(); }

    void push_back(Point const& point)
    {
      if (prefix_length.empty()) {
          prefix_length.push_back(0);
      } else {
          auto segment_distance = distance(points.back(), point);
          prefix_length.push_back(prefix_length.back() + segment_distance);
      }
  //TODO update that code
      if (points.empty()){
        extreme_points = point.bbox();
      } else {
        extreme_points += point.bbox();
      }
      points.push_back(point);
    }

    typename Points::const_iterator begin() { return points.begin(); }

    typename Points::const_iterator end() { return points.end(); }

    typename Points::const_iterator begin() const { return points.cbegin(); }

    typename Points::const_iterator end() const { return points.cend(); }

    Bbox const& bbox() const { return extreme_points; }

    FT getUpperBoundDistance(Curve const& other) const
    {
      Bbox bb = this->bbox() + other.bbox();
      return length_of_diagonal(bb);
    }


private:
    Points points;
    std::vector<distance_t> prefix_length;
    Bbox extreme_points;
};

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
