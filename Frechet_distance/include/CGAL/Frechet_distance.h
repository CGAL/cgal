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

#ifndef CGAL_FRECHET_DISTANCE_H
#define CGAL_FRECHET_DISTANCE_H

#include <CGAL/license/Frechet_distance.h>
#include <CGAL/basic.h>
#include <CGAL/Frechet_distance/internal/Frechet_distance.h>

#include <CGAL/Named_function_parameters.h>

#include <iterator>

namespace CGAL
{

/**
 * \ingroup PkgFrechetDistanceFunctions
 * determines if the Frechet distance between two polylines is larger than a given distance.
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer` with `Traits::Point_d` as value type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param polyline1 the first polyline defined by a sequence of consecutive points
 * \param polyline2 the second polyline defined by a sequence of consecutive points
 * \param distance the distance to compare against
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{Traits}
 *     \cgalParamDefault{`Traits()`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre the polylines must not be empty
 */
template <class Traits, class PointRange, class NamedParameters =  parameters::Default_named_parameters>
bool is_Frechet_distance_larger(const PointRange& polyline1,
                                const PointRange& polyline2,
                                const double distance,
                                const NamedParameters& np = parameters::default_values())
{
    constexpr bool force_filtering =
      internal_np::Lookup_named_param_def<internal_np::force_filtering_t, NamedParameters, std::false_type>::type::value;

    Traits traits =  parameters::choose_parameter<Traits>(
        parameters::get_parameter(np, internal_np::geom_traits));

    constexpr bool filtered = force_filtering ||
      std::is_same_v<typename decltype(Frechet_distance_::internal::toCurve<force_filtering>(polyline1, traits))::IFT,
                     Interval_nt<false>>;
    Protect_FPU_rounding<filtered> p;

    auto icurve1 = Frechet_distance_::internal::toCurve<force_filtering>(polyline1, traits);
    auto icurve2 = Frechet_distance_::internal::toCurve<force_filtering>(polyline2, traits);

    using distance_t = const typename decltype(icurve1)::distance_t;

    return ! Frechet_distance_::internal::lessThan(icurve1, icurve2, distance_t(distance));
}

/**
 * \ingroup PkgFrechetDistanceFunctions
 * approximates the Fréchet distance between two polylines up to an additive error
 * of `precision`.
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer` with `Traits::Point_d` as value type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param polyline1 the first polyline defined by a sequence of consecutive points
 * \param polyline2 the second polyline defined by a sequence of consecutive points
 * \param precision the precision of the approximation
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{Traits}
 *     \cgalParamDefault{`Traits()`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre the polylines must not be empty
 *
 * @return an interval enclosing the exact result, the difference between the upper and
 * the lower bound being less than `precision`.
 */
template <class Traits, class PointRange, class NamedParameters =  parameters::Default_named_parameters>
std::pair<double,double> approximate_Frechet_distance(const PointRange& polyline1,
                                                      const PointRange& polyline2,
                                                      const double precision,
                                                      const NamedParameters& np = parameters::default_values())
{
    constexpr bool force_filtering =
      internal_np::Lookup_named_param_def<internal_np::force_filtering_t, NamedParameters, std::false_type>::type::value;

    Traits traits =  parameters::choose_parameter<Traits>(
        parameters::get_parameter(np, internal_np::geom_traits));


    constexpr bool filtered = force_filtering ||
      std::is_same_v<typename decltype(Frechet_distance_::internal::toCurve<force_filtering>(polyline1, traits))::IFT,
                     Interval_nt<false>>;
    Protect_FPU_rounding<filtered> p;

    auto icurve1 = Frechet_distance_::internal::toCurve<force_filtering>(polyline1, traits);
    auto icurve2 = Frechet_distance_::internal::toCurve<force_filtering>(polyline2, traits);

    return Frechet_distance_::internal::calcDistance(icurve1, icurve2, precision);
}

}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_H
