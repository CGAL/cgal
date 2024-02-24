// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl
//                 Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_EDGE_INTERSECTION_ORACLES_3_H
#define CGAL_ISOSURFACING_3_EDGE_INTERSECTION_ORACLES_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/assertions.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{EdgeIntersectionOracle_3}
 *
 * \brief The class `Dichotomy_edge_intersection` uses a dichotomy to find the intersection point
 * between an edge and the isosurface.
 *
 * This class is for example suitable to be used as the `EdgeIntersectionOracle` template parameter of isosurfacing
 * domain classes when values are computed using an implicit function.
 * It is however not optimal when the values are interpolated from discrete values
 * since the intersection can be computed analytically in this case.
 *
 * \sa `CGAL::Isosurfacing::Linear_interpolation_edge_intersection`
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3`
 * \sa `CGAL::Isosurfacing::Dual_contouring_domain_3`
 * \sa `CGAL::Isosurfacing::Value_field_3`
 */
struct Dichotomy_edge_intersection
{
  /*!
   * \brief computes the intersection point between an edge and the isosurface.
   *
   * \tparam Domain must be a model of `IsosurfacingDomain_3`
   *
   * \param p_0 the geometric position of the first vertex of the edge
   * \param p_1 the geometric position of the second vertex of the edge
   * \param val_0 the value at the first vertex of the edge
   * \param val_1 the value at the second vertex of the edge
   * \param domain the isosurfacing domain
   * \param isovalue the isovalue defining the isosurfacing with which we seek an intersection
   * \param p the intersection point, if it exists
   *
   * \return `true` if the intersection point exists, `false` otherwise
   */
  template <typename Domain> // == Isosurfacing_domain_3 or similar
  bool operator()(const typename Domain::Geom_traits::Point_3& p_0,
                  const typename Domain::Geom_traits::Point_3& p_1,
                  const typename Domain::Geom_traits::FT val_0,
                  const typename Domain::Geom_traits::FT val_1,
                  const Domain& domain,
                  const typename Domain::Geom_traits::FT isovalue,
                  typename Domain::Geom_traits::Point_3& p) const
  {
    using Geom_traits = typename Domain::Geom_traits;
    using FT = typename Geom_traits::FT;
    using Point_3 = typename Geom_traits::Point_3;

    typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

    const bool sl = (val_0 <= isovalue);
    const bool sr = (val_1 <= isovalue);

    if(sl == sr)
      return false;

    Point_3 pl = p_0;
    Point_3 pr = p_1;

    unsigned int dichotomy_iterations = 10, iter = 0;
    const FT eps = (std::max)(FT(1e-7), std::abs(isovalue) * FT(1e-7));
    do
    {
      p = point((x_coord(pl) + x_coord(pr)) / FT(2),
                (y_coord(pl) + y_coord(pr)) / FT(2),
                (z_coord(pl) + z_coord(pr)) / FT(2));

      const FT val_p = domain.value(p);
      const bool sp = (val_p <= isovalue);

      if(sl == sp)
        pl = p;
      else if(sp == sr)
        pr = p;
      else
        break;

      if(std::abs(val_p - isovalue) < eps)
        return true;
    }
    while(++iter < dichotomy_iterations);


    return true;
  }
};

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{EdgeIntersectionOracle_3}
 *
 * \brief The class `Linear_interpolation_edge_intersection` uses linear interpolation
 * to find the intersection point between an edge and the isosurface.
 *
 * This class is for example suitable when interpolated discrete values are being used.
 *
 * \sa `CGAL::Isosurfacing::Dichotomy_edge_intersection`
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3`
 * \sa `CGAL::Isosurfacing::Dual_contouring_domain_3`
 * \sa `CGAL::Isosurfacing::Interpolated_discrete_values_3`
 */
struct Linear_interpolation_edge_intersection
{
  /*!
   * \brief computes the intersection point between an edge and the isosurface.
   *
   * \tparam Domain must be a model of `IsosurfacingDomain_3`
   *
   * \param p_0 the geometric position of the first vertex of the edge
   * \param p_1 the geometric position of the second vertex of the edge
   * \param val_0 the value at the first vertex of the edge
   * \param val_1 the value at the second vertex of the edge
   * \param domain the isosurfacing domain
   * \param isovalue the isovalue defining the isosurfacing with which we seek an intersection
   * \param p the intersection point, if it exists
   *
   * \return `true` if the intersection point exists, `false` otherwise
   */
  template <typename Domain>
  bool operator()(const typename Domain::Geom_traits::Point_3& p_0,
                  const typename Domain::Geom_traits::Point_3& p_1,
                  const typename Domain::Geom_traits::FT val_0,
                  const typename Domain::Geom_traits::FT val_1,
                  const Domain& domain,
                  const typename Domain::Geom_traits::FT isovalue,
                  typename Domain::Geom_traits::Point_3& p) const
  {
    using Geom_traits = typename Domain::Geom_traits;
    using FT = typename Geom_traits::FT;

    typename Geom_traits::Compute_x_3 x_coord = domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = domain.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = domain.geom_traits().construct_point_3_object();

    if((val_0 <= isovalue) == (val_1 <= isovalue))
      return false;

    const FT den = val_0 - val_1;
    const FT u = is_zero(den) ? 0.5 : (val_0 - isovalue) / den;
    p = point((FT(1) - u) * x_coord(p_0) + u * x_coord(p_1),
              (FT(1) - u) * y_coord(p_0) + u * y_coord(p_1),
              (FT(1) - u) * z_coord(p_0) + u * z_coord(p_1));

    return true;
  }
};

/*
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{EdgeIntersectionOracle_3}
 *
 * \brief The class `Ray_marching_edge_intersection` uses ray marching to find the intersection point
 * between an edge and the isosurface.
 *
 * This class is suitable when the values stem from a signed distance function.
 */
struct Ray_marching_edge_intersection
{
  template <typename Domain>
  bool operator()(const typename Domain::Edge_descriptor& e,
                  const Domain& domain,
                  const typename Domain::Geom_traits::FT isovalue,
                  typename Domain::Geom_traits::Point_3& p) const
  {
    // @todo this is for the case where we know domain.value is an SDF
    // then we can do better than a dichotomy
    // Take code from the AW3 sharp branch
    CGAL_assertion(false);
    return false;
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EDGE_INTERSECTION_ORACLES_3_H
