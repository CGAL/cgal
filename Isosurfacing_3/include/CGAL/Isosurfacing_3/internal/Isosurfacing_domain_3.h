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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/partition_traits.h>
#include <CGAL/Isosurfacing_3/edge_intersection_oracles_3.h>

#include <CGAL/tags.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// This class is pretty much just the concatenation of the following classes:
// - Partition: Space partitioning data structure, e.g. Cartesian grid, octree, ...
// - Values: values over the 3D space
// - Gradients: gradients over the 3D space
// - Oracle: edge-isosurface intersection computation
template <typename Partition,
          typename ValueField,
          typename GradientField,
          typename IntersectionOracle = CGAL::Isosurfacing::Dichotomy_edge_intersection>
class Isosurfacing_domain_3
{
public:
  using Edge_intersection_oracle = IntersectionOracle;

  using Geom_traits = typename Partition::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using PT = CGAL::Isosurfacing::partition_traits<Partition>;

  using vertex_descriptor = typename PT::vertex_descriptor;
  using edge_descriptor = typename PT::edge_descriptor;
  using cell_descriptor = typename PT::cell_descriptor;

  using Edge_vertices = typename PT::Edge_vertices;
  using Cells_incident_to_edge = typename PT::Cells_incident_to_edge;
  using Cell_vertices = typename PT::Cell_vertices;
  using Cell_edges = typename PT::Cell_edges;

  static constexpr Cell_type CELL_TYPE = PT::CELL_TYPE;
  static constexpr std::size_t VERTICES_PER_CELL = PT::VERTICES_PER_CELL;
  static constexpr std::size_t EDGES_PER_CELL = PT::EDGES_PER_CELL;

private:
  const Partition& m_partition;
  const ValueField& m_values;
  const GradientField& m_gradients;
  const IntersectionOracle m_intersection_oracle;

public:
  Isosurfacing_domain_3(const Partition& partition,
                        const ValueField& values,
                        const GradientField& gradients,
                        const IntersectionOracle& intersection_oracle = IntersectionOracle())
    : m_partition{partition},
      m_values{values},
      m_gradients{gradients},
      m_intersection_oracle{intersection_oracle}
  { }

  const Geom_traits& geom_traits() const
  {
    return m_partition.geom_traits();
  }

  const Edge_intersection_oracle& intersection_oracle() const
  {
    return m_intersection_oracle;
  }

public:
  // The following functions are dispatching to the partition_traits' static functions.

  // returns the position of vertex `v`
  decltype(auto) /*Point_3*/ point(const vertex_descriptor& v) const
  {
    return PT::point(v, m_partition);
  }

  // returns the value of the function at vertex `v`
  decltype(auto) /*FT*/ value(const vertex_descriptor& v) const
  {
    return m_values(v);
  }

  // returns the value of the function at point `p`
  decltype(auto) /*FT*/ value(const Point_3& p) const
  {
    return m_values(p);
  }

  // returns the gradient at point `p`
  decltype(auto) /*Vector_3*/ gradient(const Point_3& p) const
  {
    return m_gradients(p);
  }

  // returns a container with the two vertices incident to the edge `e`
  decltype(auto) /*Edge_vertices*/ incident_vertices(const edge_descriptor& e) const
  {
    return PT::incident_vertices(e, m_partition);
  }

  // returns a container with all cells incident to the edge `e`
  decltype(auto) /*Cells_incident_to_edge*/ incident_cells(const edge_descriptor& e) const
  {
    return PT::incident_cells(e, m_partition);
  }

  // returns a container with all vertices of the cell `c`
  decltype(auto) /*Cell_vertices*/ cell_vertices(const cell_descriptor& c) const
  {
    return PT::cell_vertices(c, m_partition);
  }

  // returns a container with all edges of the cell `c`
  decltype(auto) /*Cell_edges*/ cell_edges(const cell_descriptor& c) const
  {
    return PT::cell_edges(c, m_partition);
  }

  // iterates over all vertices `v`, calling `f(v)` on each of them
  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  void for_each_vertex(Functor& f) const
  {
    PT::for_each_vertex(f, m_partition, ConcurrencyTag{});
  }

  // iterates over all edges `e`, calling `f(e)` on each of them
  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  void for_each_edge(Functor& f) const
  {
    PT::for_each_edge(f, m_partition, ConcurrencyTag{});
  }

  // iterates over all cells `c`, calling `f(c)` on each of them
  template <typename ConcurrencyTag = CGAL::Sequential_tag, typename Functor>
  void for_each_cell(Functor& f) const
  {
    PT::for_each_cell(f, m_partition, ConcurrencyTag{});
  }

  // finds the intersection of the isosurface with the edge `e` (if any)
  bool construct_intersection(const Point_3& p_0, const Point_3& p_1,
                              const FT val_0, const FT val_1,
                              const FT isovalue,
                              Point_3& p) const
  {
    return m_intersection_oracle(p_0, p_1, val_0, val_1, *this, isovalue, p);
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H
