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

#ifndef CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
#define CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance/internal/curve.h>
#include <CGAL/Frechet_distance/internal/frechet_light.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>
#include <CGAL/Frechet_distance/internal/Frechet_distance_traits.h>
#include <CGAL/STL_Extension/internal/Has_nested_type_Has_filtered_predicates_tag.h>
#include <CGAL/Frechet_distance_traits_2.h>
#include <CGAL/Frechet_distance_traits_3.h>
#include <CGAL/Frechet_distance_traits_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Lazy.h>

namespace CGAL {
namespace Frechet_distance {
namespace internal {

template <class K, bool hfp = ::CGAL::internal::Has_nested_type_Has_filtered_predicates_tag<K>::value>
struct Is_filtered_kernel
{
  static constexpr bool value = K::Has_filtered_predicates_tag::value;
};

template <class K>
struct Is_filtered_kernel<K,false>
{
  static constexpr bool value = false;
};

template <bool force_filtering, class PointRange, class Traits>
auto toCurve(const PointRange& point_range, const Traits& traits)
{
  using IPoint = typename Traits::Point_d;

  if constexpr (Is_filtered_kernel<typename Kernel_traits<IPoint>::Kernel>::value)
  {
    if constexpr (std::is_floating_point_v<typename Traits::FT>)
    {
      if constexpr (Traits::Dimension::value==2)
      {
        using AK = CGAL::Simple_cartesian<Interval_nt_advanced>;
        using EK = CGAL::Simple_cartesian<Exact_rational>;

        using Filtered_traits = std::pair<Frechet_distance_traits_2<AK>, Frechet_distance_traits_2<EK>>;

        return Curve<Filtered_traits, true>(point_range, traits);
      }
      else if constexpr (Traits::Dimension::value==3)
      {
        using AK = CGAL::Simple_cartesian<Interval_nt_advanced>;
        using EK = CGAL::Simple_cartesian<Exact_rational>;

        using Filtered_traits = std::pair<Frechet_distance_traits_3<AK>, Frechet_distance_traits_3<EK>>;

        return Curve<Filtered_traits, true>(point_range, traits);
      }
      else
      {
        using AT = Frechet_distance::internal::Frechet_distance_traits<CGAL::Interval_nt_advanced, Traits::Dimension::value>;
        using ET = Frechet_distance::internal::Frechet_distance_traits<CGAL::Exact_rational, Traits::Dimension::value>;
        using Filtered_traits = std::pair<AT,ET>;

        return Curve<Filtered_traits, true>(point_range, traits);
        //return Curve<Traits, false>(point_range, traits);
      }
    }
    else
    {
      if constexpr (Traits::Dimension::value==2)
      {
        using AK = CGAL::Simple_cartesian<Interval_nt_advanced>;
        using EK = typename Kernel_traits<typename Traits::Point_d>::Kernel::Exact_kernel;

        using Filtered_traits = std::tuple<typename Traits::Point_d, Frechet_distance_traits_2<AK>, Frechet_distance_traits_2<EK>>;

        return Curve<Filtered_traits, true>(point_range, traits);
      }
      else if constexpr (Traits::Dimension::value==3)
      {
        using AK = CGAL::Simple_cartesian<Interval_nt_advanced>;
        using EK = typename Kernel_traits<typename Traits::Point_d>::Kernel::Exact_kernel;

        using Filtered_traits = std::tuple<typename Traits::Point_d, Frechet_distance_traits_3<AK>, Frechet_distance_traits_3<EK>>;

        return Curve<Filtered_traits, true>(point_range, traits);
      }
      else
      {
        //TODO: not implemented
        return Curve<Traits, false>(point_range, traits);
      }
    }
  }
  else
  {
    if constexpr (force_filtering)
    {
      if constexpr (std::is_floating_point_v<typename Traits::FT>)
      {
        if constexpr (Traits::Dimension::value==2)
        {
          using AT = Frechet_distance_traits_2<Simple_cartesian<Interval_nt_advanced>>;
          using ET = Frechet_distance_traits_2<Simple_cartesian<Exact_rational>>;

          using Filtered_traits = std::pair<AT,ET>;

          return Curve<Filtered_traits, true>(point_range, traits);
        }
        else if constexpr (Traits::Dimension::value==3)
        {
          using AT = Frechet_distance_traits_3<Simple_cartesian<Interval_nt_advanced>>;
          using ET = Frechet_distance_traits_3<Simple_cartesian<Exact_rational>>;

          using Filtered_traits = std::pair<AT,ET>;

          return Curve<Filtered_traits, true>(point_range, traits);
        }
        else
        {
          //TODO: not implemented
          return Curve<Traits, false>(point_range, traits);
        }
      }
      else
      {
        if constexpr (Traits::Dimension::value==2)
        {
          using AT = Frechet_distance_traits_2<Simple_cartesian<Interval_nt_advanced>>;
          using ET = Traits;

          using Filtered_traits = std::tuple<typename ET::Point_d,AT,ET>;

          return Curve<Filtered_traits, true>(point_range, traits, traits);
        }
        else if constexpr (Traits::Dimension::value==3)
        {
          using AT = Frechet_distance_traits_3<Simple_cartesian<Interval_nt_advanced>>;
          using ET = Traits;

          using Filtered_traits = std::tuple<typename ET::Point_d,AT,ET>;

          return Curve<Filtered_traits, true>(point_range, traits, traits);
        }
        else
        {
          //TODO: not implemented
          return Curve<Traits, false>(point_range, traits);
        }
      }
    }
    else
      return Curve<Traits, false>(point_range, traits);
  }
}

template <class Traits, bool is_filtered>
bool lessThan(Curve<Traits, is_filtered> const& curve1, Curve<Traits, is_filtered> const& curve2,
              const typename Curve<Traits, is_filtered>::distance_t& distance /*, const Traits& traits */)
{
  FrechetLight<Curve<Traits, is_filtered>> frechet;
  return frechet.lessThanWithFilters(distance, curve1, curve2);
}

template <typename Traits, bool is_filtered>
std::pair<double,double> calcDistance(Curve<Traits, is_filtered> const& curve1,
                                      Curve<Traits, is_filtered> const& curve2,
                                      double precision)
{
    FrechetLight<Curve<Traits, is_filtered>> frechet;
    return frechet.calcDistance(curve1, curve2, precision);
}

} } }  // end of namespace CGAL::Frechet_distance::internal

#endif  // CGAL_INTERNAL_Frechet_distance_FRECHET_DISTANCE_H
