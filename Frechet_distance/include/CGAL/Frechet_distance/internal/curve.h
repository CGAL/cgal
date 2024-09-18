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

template <typename T>
struct Get_rational_kernel {

    using K = typename T::Kernel;

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

    using type = decltype(get_type());

};

template <typename T, int dim>
struct Curve_base
{
    static constexpr bool is_filtered = Get_rational_kernel<T>::is_filtered;
    using distance_t = Interval_nt<false>;
    using Kernel = typename T::Kernel;
    using iKernel = Kernel_d_interface<Cartesian_base_d<distance_t,Dimension_tag<dim>>>;
    using Rational_kernel = typename Get_rational_kernel<T>::type;
    using Point = typename iKernel::Point_d;
    using Rational_point = typename Rational_kernel::Point_d;
    using Bbox = CGAL::Bbox<Dimension_tag<dim>,double>;
    using Construct_bbox = typename iKernel::Construct_bbox_d;
    using Squared_distance = typename iKernel::Squared_distance_d;
    using Difference_of_points = typename iKernel::Construct_vector_d;
    using Scaled_vector = typename iKernel::Scaled_vector_d;
    using Translated_point = typename iKernel::Translated_point_d;

    using D2D = NT_converter<distance_t,double>;
    using I2R = KernelD_converter<iKernel, Rational_kernel, Default, D2D>;

    using FT2I = NT_converter<typename Kernel::FT,distance_t>;
    using K2I = KernelD_converter<Kernel, iKernel, Default, FT2I>;
};

template <typename T>
struct Curve_base<T,2>
{
    static constexpr bool is_filtered = Get_rational_kernel<T>::is_filtered;

    using distance_t = CGAL::Interval_nt<false>;
    using Kernel = typename T::Kernel;
    using iKernel = CGAL::Simple_cartesian<distance_t>;
    using Rational_kernel = typename Get_rational_kernel<T>::type;
    using Point = typename iKernel::Point_2;
    using Rational_point = typename Rational_kernel::Point_2;
    using Bbox = Bbox_2;
    using Construct_bbox = typename iKernel::Construct_bbox_2;
    using Squared_distance = typename iKernel::Compute_squared_distance_2;
    using Difference_of_points = typename iKernel::Construct_vector_2;
    using Scaled_vector = typename iKernel::Construct_scaled_vector_2;
    using Translated_point = typename iKernel::Construct_translated_point_2;

    using D2D = NT_converter<distance_t,double>;
    using I2R = Cartesian_converter<iKernel, Rational_kernel, D2D>;

    using FT2I = NT_converter<typename Kernel::FT,distance_t>;
    using K2I = Cartesian_converter<Kernel, iKernel, FT2I>;
};

template <typename T>
struct Curve_base<T,3>
{
    static constexpr bool is_filtered = Get_rational_kernel<T>::is_filtered;
    using distance_t = CGAL::Interval_nt<false>;
    using Kernel = typename T::Kernel;
    using iKernel = CGAL::Simple_cartesian<distance_t>;
    using Rational_kernel = typename Get_rational_kernel<T>::type;
    using Point = typename iKernel::Point_3;
    using Rational_point = typename Rational_kernel::Point_3;
    using Bbox = Bbox_3;
    using Construct_bbox = typename iKernel::Construct_bbox_3;
    using Squared_distance = typename iKernel::Compute_squared_distance_3;
    using Difference_of_points = typename iKernel::Construct_vector_3;
    using Scaled_vector = typename iKernel::Construct_scaled_vector_3;
    using Translated_point = typename iKernel::Construct_translated_point_3;

    using D2D = NT_converter<distance_t,double>;
    using I2R = Cartesian_converter<iKernel, Rational_kernel, D2D>;

    using FT2I = NT_converter<typename Kernel::FT,distance_t>;
    using K2I = Cartesian_converter<Kernel, iKernel, FT2I>;
};


/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a trajectory. Additionally to the points given in the
 * input file, we also store the length of any prefix of the trajectory.
 */
template <typename T>
class Curve : public Curve_base<T, T::dimension>
{
public:
    using K = typename T::Kernel;
    using Base = Curve_base<T, T::dimension>;
    using Self  = Curve<T>;


    static constexpr bool is_filtered = Base::is_filtered;
public:
    static constexpr int dimension =  T::dimension;

    using distance_t = CGAL::Interval_nt<false>;

    using iKernel = typename Base::iKernel;

    using Point = typename Base::Point;

    using Bbox = typename Base::Bbox;

    using Construct_bbox = typename Base::Construct_bbox;

    using Squared_distance = typename Base::Squared_distance;

    using Difference_of_points = typename Base::Difference_of_points;

    using Scaled_vector = typename Base::Scaled_vector;

    using Translated_point = typename Base::Translated_point;

    using PointID = ID<Point>;
    using Points = std::vector<Point>;
    using InputPoints = std::vector<typename T::Point>;

    using Rational_kernel = typename Base::Rational_kernel;
    using Rational = typename Rational_kernel::FT;
    using Rational_point = typename Base::Rational_point;

    using I2R = typename Base::I2R;

    using K2I = typename Base::K2I;

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
        Difference_of_points difference;
        Scaled_vector scale;
        Translated_point translate;
        auto fraction = pt.getFraction().approx;
        return translate(points[pt.getPoint()] ,
                         scale(difference(points[pt.getPoint() + 1], points[pt.getPoint()]), fraction));
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
