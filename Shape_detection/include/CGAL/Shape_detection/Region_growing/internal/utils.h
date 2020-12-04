// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<typename GeomTraits>
  class Default_sqrt {

  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const {

      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits,
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {

  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) {
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {

  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) {
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  struct Compare_scores {
    const std::vector<FT>& m_scores;

    Compare_scores(const std::vector<FT>& scores) :
    m_scores(scores)
    { }

    bool operator()(const std::size_t i, const std::size_t j) const {

      CGAL_precondition(i < m_scores.size());
      CGAL_precondition(j < m_scores.size());

      return m_scores[i] > m_scores[j];
    }
  };

  template<
  typename InputRange,
  typename PointMap,
  typename Plane_3>
  void create_planes_from_points(
    const InputRange& input_range,
    const PointMap point_map,
    std::vector< std::vector<std::size_t> >& regions,
    std::vector<Plane_3>& planes) {

    using Traits = typename Kernel_traits<Plane_3>::Kernel;
    using FT = typename Traits::FT;

    using Local_traits =
    CGAL::Exact_predicates_inexact_constructions_kernel;
    using To_local_converter =
    Cartesian_converter<Traits, Local_traits>;

                using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
                using Local_plane_3 = typename Local_traits::Plane_3;

    planes.clear();
    planes.reserve(regions.size());

    std::vector<Local_point_3> points;
    const To_local_converter to_local_converter = To_local_converter();

    for (const auto& region : regions) {
      CGAL_assertion(region.size() > 0);

      points.clear();
      for (std::size_t i = 0; i < region.size(); ++i) {

        CGAL_precondition(region[i] < input_range.size());

        const auto& key = *(input_range.begin() + region[i]);
        points.push_back(to_local_converter(get(point_map, key)));
      }
      CGAL_postcondition(points.size() == region.size());

      Local_plane_3 fitted_plane;
      Local_point_3 fitted_centroid;

      CGAL::linear_least_squares_fitting_3(
        points.begin(), points.end(),
        fitted_plane, fitted_centroid,
        CGAL::Dimension_tag<0>(),
        Local_traits(),
        CGAL::Eigen_diagonalize_traits<Local_FT, 3>());

                  const Plane_3 plane = Plane_3(
        static_cast<FT>(fitted_plane.a()),
        static_cast<FT>(fitted_plane.b()),
        static_cast<FT>(fitted_plane.c()),
        static_cast<FT>(fitted_plane.d()));
      planes.push_back(plane);
    }
    CGAL_postcondition(planes.size() == regions.size());
  }

} // namespace internal

namespace RG {

  template<typename GeomTraits>
  class Plane {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Vector_3 = typename Traits::Vector_3;

    template<
    typename Input_range,
    typename Point_map>
    Plane(
      const Input_range& input_range,
      const Point_map point_map,
      const std::vector<std::size_t>& region,
      const Plane_3& plane) {

      FT x = FT(0), y = FT(0), z = FT(0);
      for (const std::size_t idx : region) {
        const auto& p = get(point_map, *(input_range.begin() + idx));
        x += p.x();
        y += p.y();
        z += p.z();
      }
      x /= static_cast<FT>(region.size());
      y /= static_cast<FT>(region.size());
      z /= static_cast<FT>(region.size());

      m_centroid = Point_3(x, y, z);
      m_base1 = plane.base1() / static_cast<FT>(CGAL::sqrt(
        CGAL::to_double(plane.base1() * plane.base1())));
      m_base2 = plane.base2() / static_cast<FT>(CGAL::sqrt(
        CGAL::to_double(plane.base2() * plane.base2())));
    }

    Point_2 to_2d(const Point_3& query) const {
      const Vector_3 v(m_centroid, query);
      return Point_2(v * m_base1, v * m_base2);
    }

    Point_3 to_3d(const Point_2& query) const {
      return m_centroid + query.x() * m_base1 + query.y() * m_base2;
    }

  private:
    Point_3 m_centroid;
    Vector_3 m_base1, m_base2;
  };

} // namespace RG
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_UTILS_H
