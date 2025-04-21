// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France), GeometryFactory (France).
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
 * \cgalModels{IsosurfacingEdgeIntersectionOracle_3}
 *
 * \brief The class `Dichotomy_edge_intersection` uses a dichotomy to find the intersection point
 * between an edge and the isosurface.
 *
 * This class is for example suitable to be used as the `EdgeIntersectionOracle` template
 * parameter of isosurfacing domain classes when values are computed using an implicit function.
 *
 * \warning It is not optimal to use this class when values are interpolated from discrete values
 * since the intersection can be computed analytically in this case.
 *
 * \sa `CGAL::Isosurfacing::Linear_interpolation_edge_intersection`
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3`
 * \sa `CGAL::Isosurfacing::Dual_contouring_domain_3`
 */
struct Dichotomy_edge_intersection
{
  unsigned int m_max_iterations;
  double m_relative_eps;

public:
  /*!
  * Constructor, enabling setting up the two criteria which can stop the dichotomy: either a
  * threshold on the value (i.e., the difference between the isovalue and the value at the current
  * point is smaller than `relative_eps * isovalue`), or a maximum number of iterations.
  */
  Dichotomy_edge_intersection(unsigned int max_iterations = 10,
                              double relative_eps = 1e-7)
    : m_max_iterations(max_iterations),
      m_relative_eps(relative_eps)
  { }

  /*!
   * \brief computes the intersection point between an edge and the isosurface.
   *
   * The result (if it exists) is stored in `p`.
   *
   * \tparam Domain must be a model of `IsosurfacingDomain_3`
   *
   * \param p_0 the location of the first vertex of the edge
   * \param p_1 the location of the second vertex of the edge
   * \param val_0 the value at the first vertex of the edge
   * \param val_1 the value at the second vertex of the edge
   * \param domain the isosurfacing domain
   * \param isovalue the isovalue defining the isosurface with which we seek an intersection
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

    unsigned int dichotomy_iterations = m_max_iterations, iter = 0;
    const FT eps = (std::max)(FT(m_relative_eps), std::abs(isovalue) * FT(m_relative_eps));
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
 * \cgalModels{IsosurfacingEdgeIntersectionOracle_3}
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
   * \param p_0 the location of the first vertex of the edge
   * \param p_1 the location of the second vertex of the edge
   * \param val_0 the value at the first vertex of the edge
   * \param val_1 the value at the second vertex of the edge
   * \param domain the isosurfacing domain
   * \param isovalue the isovalue defining the isosurface with which we seek an intersection
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

#ifndef DOXYGEN_RUNNING
/*
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{IsosurfacingEdgeIntersectionOracle_3}
 *
 * \brief The class `Ray_marching_edge_intersection` uses ray marching to find the intersection point
 * between an edge and the isosurface.
 *
 * This class is suitable when the values stem from a signed distance function.
 */
//
// @todo this is for the case where we know domain.value is an SDF
// then we could do better than a dichotomy
// see https://github.com/MaelRL/cgal/blob/AW3-Sharp_and_sparse-GF/Alpha_wrap_3/include/CGAL/Alpha_wrap_3/internal/offset_intersection.h
struct Ray_marching_edge_intersection
{
  template <typename Domain> // == Isosurfacing_domain_3 or similar
  bool operator()(const typename Domain::Geom_traits::Point_3& p_0,
                  const typename Domain::Geom_traits::Point_3& p_1,
                  const typename Domain::Geom_traits::FT val_0,
                  const typename Domain::Geom_traits::FT val_1,
                  const Domain& domain,
                  const typename Domain::Geom_traits::FT isovalue,
                  typename Domain::Geom_traits::Point_3& p) const;
};
#endif // DOXYGEN_RUNNING

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EDGE_INTERSECTION_ORACLES_3_H
