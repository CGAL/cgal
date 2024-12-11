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
#include <CGAL/Frechet_distance/internal/certificate.h>
#include <CGAL/Frechet_distance/internal/curve.h>
#include <CGAL/Frechet_distance/internal/filter.h>
#include <CGAL/Frechet_distance/internal/frechet_light_types.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>
#include <CGAL/Frechet_distance/internal/id.h>
#include <CGAL/Frechet_distance/internal/high_level_predicates.h>

#include <optional>
#include <vector>
#include <limits>

namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a
*/
  template <typename C>
class FrechetNaive
{
    // TODO: clean up
    using Curve = C;
    using K = typename Curve::K;
    using Point = typename C::Point;
    using PointID = typename Curve::PointID;
    using distance_t = typename Curve::distance_t;
    using CPoint = CGAL::Frechet_distance_::internal::CPoint<C>;
    using CInterval = CGAL::Frechet_distance_::internal::CInterval<C>;
    using CIntervals = CGAL::Frechet_distance_::internal::CIntervals<C>;
    using CIntervalsID = CGAL::Frechet_distance_::internal::CIntervalsID<C>;
    using CPosition = CGAL::Frechet_distance_::internal::CPosition<C>;
    using QSimpleInterval = CGAL::Frechet_distance_::internal::QSimpleInterval<C>;
    using QSimpleIntervals = CGAL::Frechet_distance_::internal::QSimpleIntervals<C>;
    using QSimpleOutputs = CGAL::Frechet_distance_::internal::QSimpleOutputs<C>;
    using Certificate = CGAL::Frechet_distance_::internal::Certificate<C>;
    using Filter = CGAL::Frechet_distance_::internal::Filter<C>;
    using Inputs = CGAL::Frechet_distance_::internal::Inputs<C>;
    using Outputs = CGAL::Frechet_distance_::internal::Outputs<C>;
    using BoxData = CGAL::Frechet_distance_::internal::BoxData<C>;
    using Lambda = CGAL::Frechet_distance_::internal::Lambda<C>;
    using CurvePair = std::array<Curve const*, 2>;

    // this we actually need???
    using Interval = Interval<C>;

public:
    FrechetNaive() = default;

    bool lessThan(distance_t const& distance, Curve const& curve1,
                  Curve const& curve2);
    std::pair<double,double> calcDistance(Curve const& curve1, Curve const& curve2, double epsilon);
};


template <typename C>
bool FrechetNaive<C>::lessThan(distance_t const& distance, Curve const& curve1,
                               Curve const& curve2)
{
    using OptLambda = std::optional<Lambda>;

        assert(curve1.size() >= 2);
        assert(curve2.size() >= 2);
        distance_t dist_sqr = distance * distance;

        if (CGAL::squared_distance(curve1[0],curve2[0]) > dist_sqr || CGAL::squared_distance(curve1.back(), curve2.back()) > dist_sqr) { return false; }

        std::vector<std::vector<OptLambda>> reachable1(curve1.size()-1, std::vector<OptLambda>(curve2.size(), std::nullopt));
        std::vector<std::vector<OptLambda>> reachable2(curve1.size(), std::vector<OptLambda>(curve2.size()-1, std::nullopt));
        for (size_t i = 0; i < curve1.size() - 1; ++i) {
                reachable1[i][0] = 0.;
                if (CGAL::squared_distance(curve2[0], curve1[i+1]) > dist_sqr) { break; }
        }
        for (size_t j = 0; j < curve2.size() - 1; ++j) {
                reachable2[0][j] = 0.;
                if (CGAL::squared_distance(curve1[0], curve2[j+1]) > dist_sqr) { break; }
        }

        for (size_t i = 0; i < curve1.size(); ++i) {
                for (size_t j = 0; j < curve2.size(); ++j) {
                        if (i < curve1.size() - 1 && j > 0) {
                                Interval free_int = HLPred::intersection_interval(curve2, j, curve1, i, distance);
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
                                Interval free_int = HLPred::intersection_interval(curve1, i, curve2, j, distance);
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

        assert((bool)reachable1.back().back() == (bool)reachable2.back().back());

        return (bool)reachable1.back().back();
}

template <typename C>
std::pair<double,double> FrechetNaive<C>::calcDistance(Curve const& curve1, Curve const& curve2, double epsilon)
{
        double min = 0.;
        double max = curve1.getUpperBoundDistance(curve2).sup();

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
