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

#pragma once
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance/internal/curve.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>
#include <CGAL/Frechet_distance/internal/high_level_predicates.h>

#include <optional>
#include <vector>

namespace CGAL {
namespace Frechet_distance {
namespace internal {

template <typename C>
class FrechetClassical
{
    using Curve = C;
    using distance_t = typename Curve::distance_t;
    using Lambda = CGAL::Frechet_distance::internal::Lambda<C>;

public:
    FrechetClassical() = default;

    bool lessThan(distance_t const& distance, Curve const& curve1, Curve const& curve2);
    std::pair<double,double> calcDistance(Curve const& curve1, Curve const& curve2, double epsilon);
};


template <typename C>
bool FrechetClassical<C>::lessThan(distance_t const& distance, Curve const& curve1,
                               Curve const& curve2)
{
    using OptLambda = std::optional<Lambda>;

    assert(curve1.size() >= 2);
    assert(curve2.size() >= 2);
    distance_t dist_sqr = distance * distance;

    if (Curve::squared_distance(curve1[0], curve2[0], curve1.traits()) > dist_sqr ||
        Curve::squared_distance(curve1.back(), curve2.back(), curve1.traits()) > dist_sqr) {
        return false;
    }

    std::vector<std::vector<OptLambda>> reachable1(curve1.size()-1, std::vector<OptLambda>(curve2.size(), std::nullopt));
    std::vector<std::vector<OptLambda>> reachable2(curve1.size(), std::vector<OptLambda>(curve2.size()-1, std::nullopt));

    for (size_t i = 0; i < curve1.size() - 1; ++i) {
        reachable1[i][0] = 0.;
        if (Curve::squared_distance(curve2[0], curve1[i+1], curve1.traits()) > dist_sqr) { break; }
    }
    for (size_t j = 0; j < curve2.size() - 1; ++j) {
        reachable2[0][j] = 0.;
        if (Curve::squared_distance(curve1[0], curve2[j+1], curve1.traits()) > dist_sqr) { break; }
    }

    for (size_t i = 0; i < curve1.size(); ++i) {
        for (size_t j = 0; j < curve2.size(); ++j) {
            if (i < curve1.size() - 1 && j > 0) {
                auto free_int = HLPred::intersection_interval(curve2, j, curve1, i, distance);
                if (!free_int.is_empty()) {
                    if (reachable2[i][j-1]) {
                        reachable1[i][j] = free_int.begin;
                    }
                    else if (reachable1[i][j-1] && reachable1[i][j-1] <= free_int.end) {
                        reachable1[i][j] = CGAL::max(free_int.begin, reachable1[i][j-1].value());
                    }
                }
            }
            if (j < curve2.size() - 1 && i > 0) {
                auto free_int = HLPred::intersection_interval(curve1, i, curve2, j, distance);
                if (!free_int.is_empty()) {
                    if (reachable1[i-1][j]) {
                        reachable2[i][j] = free_int.begin;
                    }
                    else if (reachable2[i-1][j] && reachable2[i-1][j] <= free_int.end) {
                        reachable2[i][j] = CGAL::max(free_int.begin, reachable2[i-1][j].value());
                    }
                }
            }
        }
    }

    assert(reachable1.back().back().has_value() == reachable2.back().back().has_value());

    return reachable1.back().back().has_value();
}

template <typename C>
std::pair<double,double> FrechetClassical<C>::calcDistance(Curve const& curve1, Curve const& curve2, double epsilon)
{
    double min = 0.;
    double max = curve1.getUpperBoundDistance(curve2);

    while (max - min >= epsilon) {
        auto split = (max + min)/2.;
        if (lessThan(split, curve1, curve2)) {
            max = split;
        }
        else {
            min = split;
        }
    }

    return std::make_pair(min, max);
}

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL
