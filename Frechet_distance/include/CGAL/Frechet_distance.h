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
#include <CGAL/Dimension.h>

#include <CGAL/Named_function_parameters.h>

#include <iterator>

namespace CGAL
{

#ifdef CGAL_EIGEN3_ENABLED
template<class Dim> struct Epick_d;
template<class Dim> struct Epeck_d;

template < typename P, int D>
struct Ambient_dimension<P, Epick_d<Dimension_tag<D>>>
{
  static constexpr int value = D;
};

template < typename P, int D>
struct Ambient_dimension<P, Epeck_d<Dimension_tag<D>>>
{
  static constexpr int value = D;
};
#endif

namespace internal
{
  template <typename P,
            int Dimension = ::CGAL::Ambient_dimension<P>::value >
  struct Get_default_traits_base
  #ifdef CGAL_EIGEN3_ENABLED
  {
    using type = Frechet_distance_traits_d<typename Kernel_traits<P>::Kernel>;
  };
  #else
  ;
  #endif

  template <typename P>
  struct Get_default_traits_base<P, 2>
  {
    using type = Frechet_distance_traits_2<typename Kernel_traits<P>::Kernel>;
  };

  template <typename P>
  struct Get_default_traits_base<P, 3>
  {
    using type = Frechet_distance_traits_3<typename Kernel_traits<P>::Kernel>;
  };

  template <class P, class K=typename Kernel_traits<P>::Kernel>
  struct Get_default_traits
  {
    using type = typename Get_default_traits_base<P>::type;
  };

  template <class P>
  struct Get_default_traits<P, ::CGAL::internal_kernel_traits::Dummy_kernel<P>>
  {
    using type = ::CGAL::internal_kernel_traits::Dummy_kernel<P>;
  };
} // end of internal namespace

/**
 * \ingroup PkgFrechetDistanceFunctions
 * determines if the Frechet distance between two polylines is larger than a given distance bound.
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer` with `Traits::Point_d` as value type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param polyline1 the first polyline defined by a sequence of consecutive points
 * \param polyline2 the second polyline defined by a sequence of consecutive points
 * \param distance_bound the distance to compare against
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a model of `FrechetDistanceTraits`}
 *     \cgalParamDefault{`Frechet_distance_traits_2`, `Frechet_distance_traits_3`, or`Frechet_distance_traits_d`, depending on the dimension of the point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre the polylines must not be empty
 */
template <class PointRange, class NamedParameters =  parameters::Default_named_parameters>
bool is_Frechet_distance_larger(const PointRange& polyline1,
                                const PointRange& polyline2,
                                const double distance_bound,
                                const NamedParameters& np = parameters::default_values())
{
    constexpr bool force_filtering =
      internal_np::Lookup_named_param_def<internal_np::force_filtering_t, NamedParameters, std::false_type>::type::value;

    using Point_d = std::remove_cv_t<std::remove_reference_t<decltype(*polyline1.begin())>>;
    using Default_traits = typename internal::Get_default_traits<Point_d>::type;
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                NamedParameters,
                                                                Default_traits>::type;
    Traits traits =  parameters::choose_parameter<Traits>(
        parameters::get_parameter(np, internal_np::geom_traits));

    constexpr bool filtered = force_filtering ||
      std::is_same_v<typename decltype(Frechet_distance::internal::toCurve<force_filtering>(polyline1, traits))::IFT,
                     Interval_nt<false>>;
    Protect_FPU_rounding<filtered> p;

    auto icurve1 = Frechet_distance::internal::toCurve<force_filtering>(polyline1, traits);
    auto icurve2 = Frechet_distance::internal::toCurve<force_filtering>(polyline2, traits);

    using distance_t = const typename decltype(icurve1)::distance_t;

    return ! Frechet_distance::internal::lessThan(icurve1, icurve2, distance_t(distance_bound));
}

/**
 * \ingroup PkgFrechetDistanceFunctions
 * returns an estimate of the Fréchet distance between the two polylines that is at most `error_bound`.
 * away from the actual Fréchet distance between the two polylines.
 *
 * \tparam Traits a model of `FrechetDistanceTraits`
 * \tparam PointRange  a model of the concept `RandomAccessContainer` with `Traits::Point_d` as value type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param polyline1 the first polyline defined by a sequence of consecutive points
 * \param polyline2 the second polyline defined by a sequence of consecutive points
 * \param error_bound a maximum bound by which the Fréchet distance estimate is allowed to deviate from the actual Fréchet distance.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a model of `FrechetDistanceTraits`}
 *     \cgalParamDefault{`Frechet_distance_traits_2`, `Frechet_distance_traits_3`, or`Frechet_distance_traits_d`, depending on the dimension of the point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre the polylines must not be empty
 *
 * @return an interval enclosing the exact result, the difference between the upper and
 * the lower bound being less than `error_bound`.
 */
template <class PointRange, class NamedParameters =  parameters::Default_named_parameters>
std::pair<double,double> approximate_Frechet_distance(const PointRange& polyline1,
                                                      const PointRange& polyline2,
                                                      const double error_bound,
                                                      const NamedParameters& np = parameters::default_values())
{
    constexpr bool force_filtering =
      internal_np::Lookup_named_param_def<internal_np::force_filtering_t, NamedParameters, std::false_type>::type::value;

    using Point_d = std::remove_cv_t<std::remove_reference_t<decltype(*polyline1.begin())>>;
    using Default_traits = typename internal::Get_default_traits<Point_d>::type;
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                NamedParameters,
                                                                Default_traits>::type;

    Traits traits =  parameters::choose_parameter<Traits>(
      parameters::get_parameter(np, internal_np::geom_traits));

    constexpr bool filtered = force_filtering ||
      std::is_same_v<typename decltype(Frechet_distance::internal::toCurve<force_filtering>(polyline1, traits))::IFT,
                     Interval_nt<false>>;
    Protect_FPU_rounding<filtered> p;

    auto icurve1 = Frechet_distance::internal::toCurve<force_filtering>(polyline1, traits);
    auto icurve2 = Frechet_distance::internal::toCurve<force_filtering>(polyline2, traits);

    return Frechet_distance::internal::calcDistance(icurve1, icurve2, error_bound);
}

}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_H
