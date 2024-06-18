// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), CNRS (France), GeometryFactory (France)
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
#include <CGAL/Frechet_distance/internal/Frechet_distance_near_neighbors_ds.h>

#include <iterator>
#include <vector>

namespace CGAL
{

// TODO: hide away in Frechet_distance::internal (different naming but nvm)
template <typename PointRange>
using PointRangeKernel = typename CGAL::Kernel_traits<
    typename std::iterator_traits<typename PointRange::iterator>::value_type>::
    Kernel;  // TODO: replace by CORE type for initial testing

template <class PointRange, class Traits = PointRangeKernel<PointRange>>
class FrechetDistanceNearNeighborsDS
{
    using PT = Traits; //  Polyline_traits_2<Traits, double>;
    using FT = typename PT::FT;
    using Point = typename PT::Point;
    using Polyline = typename PT::Polyline;
    using Polylines = typename PT::Polylines;
    using PolylineID = typename PT::PolylineID;
    using PolylineIDs = typename PT::PolylineIDs;

public:
    FrechetDistanceNearNeighborsDS() = default;

    void insert(const Polyline& curve);
    void insert(const Polylines& curve);

    PolylineIDs get_close_curves(const Polyline& curve, double distance);

private:
    Polylines curves;
    Frechet_distance::internal::FrechetKdTree<Traits> kd_tree;
};

// TODO: store preprocessed curves after CGALization
template <class PointRange, class Traits>
void FrechetDistanceNearNeighborsDS<PointRange, Traits>::insert(
    const Polylines& curves)
{
    this->curves = curves;  // FIXME: copies all the curves...

    kd_tree.insert(curves);
    kd_tree.build();
}

template <class PointRange, class Traits>
auto FrechetDistanceNearNeighborsDS<PointRange, Traits>::get_close_curves(
    const Polyline& curve, double distance) -> PolylineIDs
{
    auto result = kd_tree.search(curve, distance);

    auto predicate = [&](PolylineID id) {
        return ! continuous_Frechet_distance_less_than<Traits, PointRange>(
            curve, curves[id], distance);
    };
    auto new_end = std::remove_if(result.begin(), result.end(), predicate);
    result.erase(new_end, result.end());

    return result;
}

}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_NEAR_NEIGHBORS_DS_H
