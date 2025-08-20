// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shepard Liu	 <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_COORDINATE_CONVERTER_H
#define CGAL_DRAW_AOS_ARR_COORDINATE_CONVERTER_H
#include <cmath>

#include <CGAL/number_type_config.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/*!
 * \brief class handling coordinate conversion between 2D parameterized surface coordinates and cartesian coordinates.
 *
 * \tparam GeomTraits
 */
template <typename GeomTraits>
class Arr_coordinate_converter
{
  using Geom_traits = GeomTraits;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point = typename Approx_traits::Point;

public:
  Arr_coordinate_converter(const GeomTraits& traits)
      : m_traits(traits) {}

  /*!
   * \brief Converts a point in cartesian coordinates to parameterized surface coordinates.
   *
   * \param pt
   * \return Point
   */
  Point to_uv(Approx_point pt) const { return pt; }

  /*!
   * \brief Converts a point in parameterized surface coordinates to cartesian coordinates.
   *
   * \param pt
   * \return Approx_point
   */
  Approx_point to_cartesian(Point pt) const { return pt; }

private:
  const GeomTraits& m_traits;
};

/*!
 * \brief Converter specialization for geodesic arc on sphere traits.
 *
 * Provides conversions between spherical coordinates and right-handed Cartesian coordinates. Sphercial coordinates are
 * represented as azimuth ( [0, 2 Pi) ) and polar ( [0, Pi] ) angle in radians. Points on the identification curve have
 * azimuth == 0. The south pole has polar == 0.
 *
 * \tparam Kernel
 * \tparam atanX
 * \tparam atanY
 */
template <typename Kernel, int atanX, int atanY>
class Arr_coordinate_converter<Arr_geodesic_arc_on_sphere_traits_2<Kernel, atanX, atanY>>
{
  using Geom_traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Point = typename Approx_traits::Point;

public:
  Arr_coordinate_converter(const Geom_traits& traits)
      : m_traits(traits) {}

  Point to_uv(Approx_point point) const {
    if(point.location() == Approx_point::MAX_BOUNDARY_LOC) return Point(0, CGAL_PI);
    if(point.location() == Approx_point::MIN_BOUNDARY_LOC) return Point(0, 0);
    Approx_nt azimuth_from_id =
        std::fmod(std::atan2(point.dy(), point.dx()) - std::atan2(atanY, atanX) + 2 * CGAL_PI, 2 * CGAL_PI);
    return Point(azimuth_from_id, std::acos(-point.dz()));
  }

  Approx_point to_cartesian(Point point) const {
    using Direction_3 = typename Geom_traits::Approximate_kernel::Direction_3;

    Approx_nt polar = point.y();
    if(point.y() == CGAL_PI) return Approx_point(Direction_3(0, 0, 1), Approx_point::MAX_BOUNDARY_LOC);
    if(point.y() == 0) return Approx_point(Direction_3(0, 0, -1), Approx_point::MIN_BOUNDARY_LOC);
    Approx_nt azimuth = point.x() + std::atan2(atanY, atanX);
    Approx_nt x = std::sin(polar) * std::cos(azimuth);
    Approx_nt y = std::sin(polar) * std::sin(azimuth);
    Approx_nt z = -std::cos(polar);
    Direction_3 dir(x, y, z);
    return Approx_point(dir, azimuth == 0 ? Approx_point::MID_BOUNDARY_LOC : Approx_point::NO_BOUNDARY_LOC);
  }

private:
  const Geom_traits& m_traits;
};

} // namespace draw_aos
} // namespace CGAL

#endif