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

#ifndef CGAL_ISOSURFACING_3_VALUE_FUNCTION_3_H
#define CGAL_ISOSURFACING_3_VALUE_FUNCTION_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/partition_traits.h>

#include <functional>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domain_helpers_grp
 *
 * \cgalModels{IsosurfacingValueField_3}
 *
 * \brief The class `Value_function_3` represents a field of scalars computed
 * using a user-provided unary function.
 *
 * \tparam Partition must be a model of `IsosurfacingPartition_3`
 *
 * \sa `CGAL::Isosurfacing::Marching_cubes_domain_3`
 * \sa `CGAL::Isosurfacing::Dual_contouring_domain_3`
 */
template <typename Partition>
class Value_function_3
{
public:
  using Geom_traits = typename Partition::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  using PT = partition_traits<Partition>;
  using vertex_descriptor = typename PT::vertex_descriptor;

private:
  std::function<FT(const Point_3&)> m_fn;
  const Partition& m_partition;

public:
  /**
   * \brief constructs a field of values using a value function and a partition.
   *
   * \tparam Function must provide the following function signature:
   *                  `FT operator()(const %Point_3&) const`
   *
   * \param fn the function providing values
   * \param partition the space partitioning data structure
   */
  template <typename Function>
  Value_function_3(const Function& fn,
                   const Partition& partition)
    : m_fn{fn},
      m_partition{partition}
  { }

public:
  /**
   * \brief evaluates the function at the point `p`
   */
  FT operator()(const Point_3& p) const
  {
    return m_fn(p);
  }

  /**
   * \brief evaluates the function at the vertex `v`
   */
  FT operator()(const vertex_descriptor& v) const
  {
    return this->operator()(PT::point(v, m_partition));
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_VALUE_FUNCTION_3_H
