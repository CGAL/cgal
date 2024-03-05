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

#ifndef CGAL_ISOSURFACING_3_DUAL_CONTOURING_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_DUAL_CONTOURING_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/edge_intersection_oracles_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \cgalModels{IsosurfacingDomainWithGradient_3}
 *
 * \brief A domain that can be used as input in the %Dual Contouring algorithm.
 *
 * \details This class is essentially wrapper around the different bricks provided by its
 * template parameters: `Partition` provides the spacial partitioning, `ValueField` and `GradientField`
 * the values and gradients that define the isosurface. The optional template parameter
 * `EdgeIntersectionOracle` gives control over the method used to computate edge-isosurface intersection points.
 *
 * \tparam Partition must be a model of `IsosurfacingPartition_3`
 * \tparam ValueField must be a model of `IsosurfacingValueField_3`
 * \tparam GradientField must be a model of `IsosurfacingGradientField_3`
 * \tparam EdgeIntersectionOracle must be a model of `IsosurfacingEdgeIntersectionOracle_3`
 *
 * \sa `CGAL::Isosurfacing::dual_contouring()`
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3()`
 */
template <typename Partition,
          typename ValueField,
          typename GradientField,
          typename EdgeIntersectionOracle = Dichotomy_edge_intersection>
class Dual_contouring_domain_3
#ifndef DOXYGEN_RUNNING
  : public internal::Isosurfacing_domain_3<Partition, ValueField, GradientField, EdgeIntersectionOracle>
#endif
{
private:
  using Base = internal::Isosurfacing_domain_3<Partition, ValueField, GradientField, EdgeIntersectionOracle>;

public:
  /**
   * \brief constructs a domain that can be used with the %Dual Contouring algorithm.
   *
   * \param partition the space partitioning data structure
   * \param values a continuous field of scalar values, defined over the geometric span of `partition`
   * \param gradients a continuous field of normalized vectors, defined over the geometric span of `partition`
   * \param intersection_oracle the oracle for edge-isosurface intersection computation
   *
   * \warning the domain class keeps a reference to the `partition`, `values` and `gradients` objects.
   * As such, users must ensure that the lifetime of these objects exceeds that of the domain object.
   */
  Dual_contouring_domain_3(const Partition& partition,
                           const ValueField& values,
                           const GradientField& gradients,
                           const EdgeIntersectionOracle& intersection_oracle = EdgeIntersectionOracle())
    : Base(partition, values, gradients, intersection_oracle)
  { }
};

/**
 * \ingroup IS_Domains_grp
 *
 * \brief creates a new instance of a domain that can be used with the %Dual Contouring algorithm.
 *
 * \tparam Partition must be a model of `IsosurfacingPartition_3`
 * \tparam ValueField must be a model of `IsosurfacingValueField_3`
 * \tparam GradientField must be a model of `IsosurfacingGradientField_3`
 * \tparam EdgeIntersectionOracle must be a model of `IsosurfacingEdgeIntersectionOracle_3`
 *
 * \param partition the space partitioning data structure
 * \param values a continuous field of scalar values, defined over the geometric span of `partition`
 * \param gradients a continuous field of normalized vectors, defined over the geometric span of `partition`
 * \param intersection_oracle the oracle for edge-isosurface intersection computation
 *
 * \warning the domain class keeps a reference to the `partition`, `values` and `gradients` objects.
 * As such, users must ensure that the lifetime of these objects exceeds that of the domain object.
 */
template <typename Partition,
          typename ValueField,
          typename GradientField,
          typename EdgeIntersectionOracle = Dichotomy_edge_intersection>
Dual_contouring_domain_3<Partition, ValueField, GradientField, EdgeIntersectionOracle>
create_dual_contouring_domain_3(const Partition& partition,
                                const ValueField& values,
                                const GradientField& gradients,
                                const EdgeIntersectionOracle& intersection_oracle = EdgeIntersectionOracle())
{
  return { partition, values, gradients, intersection_oracle };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_DUAL_CONTOURING_DOMAIN_3_H
