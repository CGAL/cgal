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
  : public Curve<Approximate_traits, false>
{
public:
    using Base = Curve<Approximate_traits, false>;
    using distance_t = typename Approximate_traits::FT;
    using Point = typename Approximate_traits::Point_d;
    using PointID = typename Base::PointID;

    using Rational = typename Rational_traits::FT;
    using Rational_point = typename Rational_traits::Point_d;

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

        this->points.reserve(point_range.size());
        K2I convert;
        for (auto const& p : point_range)
          this->points.push_back(convert(p));

        this->init();
    }

    Rational_point rpoint(PointID const& i) const
    {
        //~ if constexpr (is_floating_point) {
            I2R convert;
            return convert(this->point(i));
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
};

template <typename T>
class Curve<T, false>
{
public:
    using Traits = T;

  //TODO: we expect Dimension_tag for now.
  static constexpr int dimension =  Traits::Dimension::value;
  using Bbox = std::conditional_t<dimension==2,
                                  Bbox_2, std::conditional_t<dimension==3,
                                                             Bbox_3, ::CGAL::Bbox<typename Traits::Dimension, double>>>;


//////
    using FT = typename Traits::FT;
    using IFT = std::conditional_t<std::is_floating_point_v<FT>, FT, CGAL::Interval_nt<false>>;
    using distance_t = FT;

    using Point = typename Traits::Point_d;

    using Construct_bbox = typename Traits::Construct_bbox_d;
    using Construct_cartesian_const_iterator_d = typename Traits::Construct_cartesian_const_iterator_d;

    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<Point>;


    Curve() = default;

    void init()
    {
      if (points.empty()) return;

      prefix_length.resize(points.size());
      prefix_length[0] = 0;

      Construct_bbox bbox = traits_.construct_bbox_d_object();
      extreme_points = bbox(point(0));

      for (PointID i = 1; i < points.size(); ++i)
      {
        IFT segment_distance = distance(point(i - 1), point(i), traits_);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points += bbox(point(i));
      }
    }

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
    Curve(const PointRange& point_range, const Traits& traits = Traits())
     : points(point_range.begin(), point_range.end())
     , traits_(traits)
    {
      init();
    }

    //  Curve(const Points& points);

    void reserve(std::size_t n) { points.reserve(n); }

    std::size_t size() const { return points.size(); }

    bool empty() const { return points.empty(); }

    Point const& operator[](PointID const& i) const { CGAL_assertion(i<points.size()); return points[i]; }

    Point const& point(PointID const& i) const { return points[i]; }

    bool operator==(Curve const& other) const
    {
        return std::equal(points.cbegin(), points.cend(), other.points.cbegin(),
                          other.points.cend());
    }

    bool operator!=(Curve const& other) const { return !(*this == other); }

    static IFT squared_distance(Point const& p, Point const& q, const Traits& traits)
    {
      IFT res(0);
      Construct_cartesian_const_iterator_d ccci = traits.construct_cartesian_const_iterator_d_object();

      auto itp = ccci(p), itq = ccci(q);

      for (int d=0;d<Traits::Dimension::value; ++d)
      {
        res+=square(to_ift(*itp)-to_ift(*itq));
        ++itp; ++itq;
      }
      return res;
    }
    static IFT distance(Point const& p, Point const& q, const Traits& traits)
    {
      return sqrt(squared_distance(p,q, traits));
    }
    static IFT to_ift(const FT& n)
    {
      if constexpr(std::is_floating_point_v<FT>)
      {
        return n;
      }
      else
      {
        return IFT(to_interval(n));
      }
    }
    static FT inf(const IFT& n)
    {
      if constexpr(std::is_floating_point_v<FT>)
      {
        return n;
      }
      else
      {
        return n.inf();
      }
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
        Construct_cartesian_const_iterator_d ccci = traits_.construct_cartesian_const_iterator_d_object();
        auto itp = ccci(point(pt.getPoint()));
        auto itq = ccci(point(pt.getPoint()+1));

        for (int d=0; d<dimension; ++d)
          b_coords[d] = one_m_f * (*itp++) +  (*itq++) * fraction;

        if constexpr (dimension==2)
          return Point(b_coords[0], b_coords[1]);
        else if constexpr (dimension==3)
          return Point(b_coords[0], b_coords[1], b_coords[2]);
        else
          return Point(b_coords.begin(), b_coords.end());
    }

    Point interpolate_at(PointID const& pt) const { return point(pt); }

    IFT curve_length(PointID const& i, PointID const& j) const
    {
      CGAL_assertion(i>=0 && i<prefix_length.size());
      CGAL_assertion(j>=0 && j<prefix_length.size());

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
      Construct_bbox bbox;

      if (points.empty()){
        extreme_points = bbox(point);
      } else {
        extreme_points += bbox(point);
      }
      points.push_back(point);
    }

    typename Points::const_iterator begin() { return points.begin(); }

    typename Points::const_iterator end() { return points.end(); }

    typename Points::const_iterator begin() const { return points.cbegin(); }

    typename Points::const_iterator end() const { return points.cend(); }

    Bbox const& bbox() const { return extreme_points; }

    double getUpperBoundDistance(Curve const& other) const
    {
      Bbox bb = this->bbox() + other.bbox();
      return length_of_diagonal(bb);
    }

    const Traits& traits() const
    {
      return traits_;
    }

protected:

    Points points;
    std::vector<IFT> prefix_length;
    Bbox extreme_points;
    Traits traits_;
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
