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

#ifndef CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
#define CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/basic.h>
#include <CGAL/number_utils.h>

#include <array>
#include <iterator>
#include <tuple>
#include <vector>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {


template <class Traits>
class FrechetKdTree
{
    using PT = Traits; // Polyline_traits_2<Traits, double>;
    using FT = typename PT::FT;
    using Point = typename PT::Point_d;
    using Polyline = std::vector<Point>;
    using Polylines = std::vector<Polyline>;
    using PolylineID = std::size_t;
    using PolylineIDs = std::vector<PolylineID>;

    using D = Dimension_tag<4*Traits::Dimension::value>;
    // FIXME: is fixing Cartesian_d too non-general here?
    using Traits_d = Cartesian_d<typename Traits::FT>;
    using Tree_traits_base = Search_traits_d<Traits_d, D>;
    using Point_d = typename Tree_traits_base::Point_d;
    using Point_and_id = std::tuple<Point_d, std::size_t>;
    using Tree_traits = CGAL::Search_traits_adapter<
        Point_and_id, CGAL::Nth_of_tuple_property_map<0, Point_and_id>,
        Tree_traits_base>;
    using Tree = Kd_tree<Tree_traits>;

public:
    FrechetKdTree() = default;

    void insert(Polylines const& curves);  // TODO: also add other methods of
                                           // inserting points
    void build();
    PolylineIDs search(Polyline const& curve, FT distance);

private:
    Kd_tree<Tree_traits> kd_tree;

    static Point_d to_kd_tree_point(const Polyline& curve);

    class QueryItem
    {
        using D = FrechetKdTree::D;
        using Point_d = FrechetKdTree::Point_d;
        using Point_and_id = FrechetKdTree::Point_and_id;

        // const Polyline& curve;
        const Point_d p;
        const FT distance;
        const FT distance_sqr;

    public:
        QueryItem(Polyline const& curve, FT distance)
            : p(to_kd_tree_point(curve)),
              distance(distance),
              distance_sqr(distance * distance)
        {
        }

        bool contains(Point_and_id const& point) const
        {
            auto const& q = std::get<0>(point);

            for (size_t i = 0; i < 4; i += 2) {
                // AF deal with dimension
                auto a = Point(p[i], p[i + 1]);
                auto b = Point(q[i], q[i + 1]);
                // AF: In case Point stays the input point type we have
                // to robustify with interval arithmetic
                // here: certainly
                if (compare_squared_distance(a, b, distance_sqr) == LARGER) {
                    return false;
                }
            }
            for (size_t i = 4; i < 8; ++i) {
                // AF: certainly
                if (CGAL::abs(p[i] - q[i]) > distance) {
                    return false;
                }
            }

            return true;
        }
        bool inner_range_intersects(const Kd_tree_rectangle<FT, D>& rect) const
        {
            for (int d = 0; d < D::value; ++d) {
                if (rect.min_coord(d) > p[d] + distance &&
                    rect.max_coord(d) + distance < p[d]) {
                        // AF: certainly
                    return false;
                }
            }
            return true;
        }
        bool outer_range_contains(const Kd_tree_rectangle<FT, D>& rect) const
        {
            // check if rect[0:2] is contained in circle of start point AND
            // check if rect[2:4] is contained in circle of end point
            for (size_t i = 0; i < 4; i += 2) {
                // TODO: this is a manual test if a rectangle is contained in a
                // circle. Does CGAL offer anything handy for that?
                // AF: deal with dimension
                auto point = Point(p[i], p[i + 1]);
                for (auto x : {rect.min_coord(i), rect.max_coord(i)}) {
                    for (auto y :
                         {rect.min_coord(i + 1), rect.max_coord(i + 1)}) {
                        if (compare_squared_distance(Point(x, y), point,
                                                     distance_sqr) == LARGER) {
                                                        // AF: certainly
                            return false;
                        }
                    }
                }
            }

            // rect[4:8] is contained in intersection rect (see notebook)
            // TODO: this is a manual test if a rectangle is contained in
            // another rectangle. Does CGAL offer anything handy for that?
            for (std::size_t i = 4; i < 8; ++i) {
                if (p[i] - distance > rect.min_coord(i) ||
                    p[i] + distance < rect.max_coord(i)) {
                        // AF: certainly
                    return false;
                }
            }

            return true;
        }
    };
};

template <class Traits>
auto FrechetKdTree<Traits>::to_kd_tree_point(const Polyline& curve) -> Point_d
{
    CGAL_precondition(!curve.empty());

    // TODO: remove this when curves are preprocessed first
    FT x_min, y_min, x_max, y_max;
    x_min = y_min = (std::numeric_limits<FT>::max)();
    x_max = y_max = (std::numeric_limits<FT>::min)();
    for (auto const& point : curve) {
      x_min = (std::min)(x_min, point.x());
      y_min = (std::min)(y_min, point.y());
      x_max = (std::max)(x_max, point.x());
      y_max = (std::max)(y_max, point.y());
    }

    std::array<FT, D::value> values = {curve.front().x(),
                                       curve.front().y(),
                                       curve.back().x(),
                                       curve.back().y(),
                                       x_min,
                                       y_min,
                                       x_max,
                                       y_max};
    return Point_d(D::value, values.begin(), values.end());
}

template <class Traits>
void FrechetKdTree<Traits>::insert(Polylines const& curves)
{
    for (std::size_t id = 0; id < curves.size(); ++id) {
        auto kd_tree_point = to_kd_tree_point(curves[id]);
        kd_tree.insert(Point_and_id{kd_tree_point, id});
    }
}

template <class Traits>
void FrechetKdTree<Traits>::build()
{
    kd_tree.build();
}

template <class Traits>
auto FrechetKdTree<Traits>::search(Polyline const& curve, FT distance)
    -> PolylineIDs
{
    // TODO: This is the best way I found to not copy the 8-dimensional point,
    // but only the id. Is there an even better way?
    class back_insert_it
    {
        PolylineIDs* curve_ids;

    public:
        back_insert_it(PolylineIDs& curve_ids) : curve_ids(&curve_ids) {}
        back_insert_it(back_insert_it& it) = default;
        back_insert_it(back_insert_it&& it) = default;

        back_insert_it& operator*() { return *this; }
        back_insert_it& operator++() { return *this; }
        back_insert_it& operator++(int) { return *this; }
        back_insert_it& operator=(back_insert_it const& it) = default;
        back_insert_it& operator=(const Point_and_id& point_and_id)
        {
            curve_ids->push_back(std::get<1>(point_and_id));
            return *this;
        }
    };

    PolylineIDs result;
    kd_tree.search(back_insert_it(result), QueryItem(curve, distance));

    return result;
}

} // namespace internal
} // namespace Frechet_distance_
}  // end of namespace CGAL

#endif  // CGAL_INTERNAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
