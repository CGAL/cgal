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
#include <CGAL/Dimension.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/Bbox_d.h>
#include <array>
#include <iterator>
#include <tuple>
#include <vector>

namespace CGAL {
namespace Frechet_distance {
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
        using BB = Bbox_d<Dimension_tag<dim>>;
        Point ends[2];
        BB endsbbox;
        BB bbox;
        // using PP = Concatenate_iterator<typename Point::Cartesian_const_iterator, typename Point::Cartesian_const_iterator>;
        using Bbcci = typename BB::Cartesian_const_iterator;
        using Cartesian_const_iterator_d = Concatenate_iterator<Bbcci, Bbcci>;

        Cartesian_const_iterator_d cartesian_begin() const
        {
            Bbcci ppb = endsbbox.cartesian_begin();
            Bbcci ppe = endsbbox.cartesian_end();
          return Cartesian_const_iterator_d(ppe, bbox.cartesian_begin(), ppb);
        }

        Cartesian_const_iterator_d cartesian_end() const
        {
            Bbcci ppe = endsbbox.cartesian_end();
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

        // AF: should that go in the traits class?
        Point construct(const std::array<FT,dim>& coords) const
        {
          if constexpr (dim==2){ return Point(coords[0], coords[1]);}
          else if constexpr (dim==3){ return Point(coords[0], coords[1], coords[2]);}
          else { return Point(coords);}
        }

        // AN: here we cannot really do anything wrong; if we return false, we continue searching;
        // if we return true, we return all points inside which is also fine as we check them later.
        bool outer_range_contains(const Kd_tree_rectangle<FT, D>& rect) const
        {
            typename Point_d::Cartesian_const_iterator_d pb = p.cartesian_begin();
            std::array<FT,dim> cornercoords;
            {
              const Point& point = p.ends[0];
              for(size_t i = 0; i < dim; ++i, ++pb){
                  cornercoords[i] = ((rect.min_coord(i) + rect.max_coord(i))/2 > *pb) ?
                                      rect.max_coord(i) : rect.min_coord(i);
              }
              Point corner = construct(cornercoords);
              if (Compare_squared_distance()(corner, point,
                                             distance_sqr) == LARGER) {
                return false;
              }
            }
            {
              const Point& point = p.ends[1];
              for(size_t i = 0; i < dim; ++i, ++pb){
                  cornercoords[i] = ((rect.min_coord(i+dim) + rect.max_coord(i+dim))/2 > *pb) ?
                                      rect.max_coord(i+dim) : rect.min_coord(i+dim);
              }
              Point corner = construct(cornercoords);
              if (Compare_squared_distance()(corner, point,
                                             distance_sqr) == LARGER) {
                return false;
              }
            }

            // rect[2*dim : 4*dim] is contained in intersection rect (see notebook)
            // TODO: this is a manual test if a rectangle is contained in
            // another rectangle. Does CGAL offer anything handy for that?
            for (std::size_t i = 2*dim; i < 4*dim; ++i, ++pb) {
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
    res.endsbbox = bbox(res.ends[0]) + bbox(res.ends[1]);
    for (auto const& point : curve) {
        Bbox_d<Dimension_tag<dim>> bb = bbox(point);
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
} // namespace Frechet_distance
}  // end of namespace CGAL


