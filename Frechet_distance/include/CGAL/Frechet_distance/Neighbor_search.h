// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), CNRS (France), GeometryFactory (France)
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

#ifndef CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
#define CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H

#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/basic.h>
#include <CGAL/Frechet_distance/internal/Neighbor_search.h>

#include <iterator>
#include <vector>

namespace CGAL {
namespace Frechet_distance {

// TODO: hide away in Frechet_distance_::internal (different naming but nvm)
template <typename PointRange>
using PointRangeKernel = typename CGAL::Kernel_traits<
    typename std::iterator_traits<typename PointRange::iterator>::value_type>::
    Kernel;  // TODO: replace by CORE type for initial testing

template <class PointRange, class Traits = PointRangeKernel<PointRange>>
class Neighbor_search
{
    using PT = Traits; //  Polyline_traits_2<Traits, double>;
    using FT = typename PT::FT;
    using Point = typename PT::Point_d;
    using Polyline = std::vector<Point>;
    using Polylines = std::vector<Polyline>;
    using PolylineID = std::size_t;
    using PolylineIDs = std::vector<PolylineID>;

public:
    Neighbor_search() = default;

    void insert(const Polyline& curve);
    void insert(const Polylines& curves);

    PolylineIDs get_close_curves(const Polyline& curve, double distance);

private:
    Polylines curves;
    Frechet_distance::internal::FrechetKdTree<Traits> kd_tree;
};

// TODO: store preprocessed curves after CGALization
template <class PointRange, class Traits>
void Neighbor_search<PointRange, Traits>::insert(
    const Polylines& curves)
{
    this->curves = curves;  // FIXME: copies all the curves...

    kd_tree.insert(curves);
    kd_tree.build();
}

template <class PointRange, class Traits>
auto Neighbor_search<PointRange, Traits>::get_close_curves(
    const Polyline& curve, double distance) -> PolylineIDs
{
    auto result = kd_tree.search(curve, distance);

    auto predicate = [&](PolylineID id) {
      CGAL_assertion(id < curves.size());
        return  is_Frechet_distance_larger(
            curve, curves[id], distance, parameters::geom_traits(Traits()));
    };
    std::cout << result.size() << " curves returned from tree" << std::endl;
    auto new_end = std::remove_if(result.begin(), result.end(), predicate);
    result.erase(new_end, result.end());

    return result;
}

}
}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
