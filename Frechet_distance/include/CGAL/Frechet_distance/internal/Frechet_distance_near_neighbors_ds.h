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
#include <CGAL/Dimension.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/Bbox.h>
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
    using PT = Traits; // model of FrechetDistanceTraits
    using FT = typename PT::FT;
    using Point = typename PT::Point_d;
    using Construct_cartesian_const_iterator = typename PT::Construct_cartesian_const_iterator_d;
    using Compare_squared_distance = typename PT::Compare_squared_distance_d;
    using Construct_bbox = typename PT::Construct_bbox_d;
    using Polyline = std::vector<Point>;
    using Polylines = std::vector<Polyline>;
    using PolylineID = std::size_t;
    using PolylineIDs = std::vector<PolylineID>;

    static constexpr int dim = Traits::Dimension::value;

    using D = Dimension_tag<4*dim>;


    struct Point_d {
        using BB = Bbox<Dimension_tag<dim>,FT>;
        Point ends[2];
        BB bbox;
        using PP = Concatenate_iterator<typename Point::Cartesian_const_iterator, typename Point::Cartesian_const_iterator>;
        using Cartesian_const_iterator_d = Concatenate_iterator<PP, typename BB::Cartesian_const_iterator>;

        Cartesian_const_iterator_d cartesian_begin() const
        {
            Construct_cartesian_const_iterator ccc;
            PP ppb(ccc(ends[0],0), ccc(ends[1]), ccc(ends[0]));
            PP ppe(ccc(ends[0],0), ccc(ends[1]), ccc(ends[1],0));
          return Cartesian_const_iterator_d(ppe, bbox.cartesian_begin(), ppb);
        }

        Cartesian_const_iterator_d cartesian_end() const
        {
            Construct_cartesian_const_iterator ccc;
            PP ppe(ccc(ends[0],0), ccc(ends[1]), ccc(ends[1],0));
            return Cartesian_const_iterator_d(ppe, bbox.cartesian_begin(), bbox.cartesian_end(), 0);
        }

        struct Construct_cartesian_const_iterator_d {
          Cartesian_const_iterator_d operator()(const Point_d& p) const
          {
            return p.cartesian_begin();
          }
          Cartesian_const_iterator_d operator()(const Point_d& p, int) const
          {
            return p.cartesian_end();
          }
        };
    };

    using Tree_traits_base = Search_traits<FT, Point_d, typename Point_d::Cartesian_const_iterator_d, typename Point_d::Construct_cartesian_const_iterator_d, D>;
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

            for (size_t i = 0; i < 2; ++i) {
                const Point& a = p.ends[i];
                const Point& b = q.ends[i];
                // AF: In case Point stays the input point type we have
                // to robustify with interval arithmetic
                // here: certainly
                // AN: yes, here we need certainly!
                if (Compare_squared_distance()(a, b, distance_sqr) == LARGER) {
                    return false;
                }
            }
            typename Point_d::BB::Cartesian_const_iterator pb = p.bbox.cartesian_begin();
            typename Point_d::BB::Cartesian_const_iterator qb = q.bbox.cartesian_begin();
            for (size_t i = 0; i < 2*dim; ++i, ++pb, ++qb) {
                // AF: certainly
                // AN: yes, here we need certainly!
                if (CGAL::abs(*pb - *qb) > distance) {
                    return false;
                }
            }

            return true;
        }
        bool inner_range_intersects(const Kd_tree_rectangle<FT, D>& rect) const
        {
            typename Point_d::Cartesian_const_iterator_d pb = p.cartesian_begin();
            for (int d = 0; d < D::value; ++d, ++pb) {
                if (rect.min_coord(d) > *pb + distance &&
                    rect.max_coord(d) + distance < *pb) {
                    // AF: certainly
                    // AN: yes, here is really certainly!
                    return false;
                }
            }
            return true;
        }

        // AN: here we cannot really do anything wrong; if we return false, we continue searching;
        // if we return true, we return all points inside which is also fine as we check them later.
        bool outer_range_contains(const Kd_tree_rectangle<FT, D>& rect) const
        {
            // check if rect[0 : dim/2] is contained in circle of start point AND
            // check if rect[dim/2 : dim] is contained in circle of end point
            for (size_t i = 0; i < 4; i += 2) {
                // TODO: this is a manual test if a rectangle is contained in a
                // circle. Does CGAL offer anything handy for that?
                // AF: deal with dimension
                const Point& point = p.ends[i/2];
                for (auto x : {rect.min_coord(i), rect.max_coord(i)}) {
                    for (auto y : {rect.min_coord(i + 1), rect.max_coord(i + 1)}) {
                        Point corner(x,y,0);
                        if (Compare_squared_distance()(corner, point,
                                                       distance_sqr) == LARGER) {
                            return false;
                        }
                    }
                }
            }

            // rect[4:8] is contained in intersection rect (see notebook)
            // TODO: this is a manual test if a rectangle is contained in
            // another rectangle. Does CGAL offer anything handy for that?
            typename Point_d::BB::Cartesian_const_iterator pb = p.bbox.cartesian_begin();
            for (std::size_t i = dim/2; i < dim; ++i, ++pb) {
                if (*pb - distance > rect.min_coord(i) ||
                    *pb + distance < rect.max_coord(i)) {
                    return false;
                }
            }

            // AN: This has to be certainly!!
            return true;
        }
    };
};

template <class Traits>
auto FrechetKdTree<Traits>::to_kd_tree_point(const Polyline& curve) -> Point_d
{
    CGAL_precondition(!curve.empty());

    Construct_bbox bbox = Traits().construct_bbox_d_object(); //  AF Do we have to pass a Traits object?
    Point_d res;

    res.ends[0] = curve.front();
    res.ends[1] = curve.back();
    for (auto const& point : curve) {
        Bbox<Dimension_tag<dim>,FT> bb = bbox(point);
        res.bbox +=  bb;
    }
    return res;
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
