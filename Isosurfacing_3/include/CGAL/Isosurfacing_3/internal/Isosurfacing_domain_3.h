// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// A wrapper class to puzzle a domain together from different combinations of topology,
// geometry, function, and gradient.
template <typename GeomTraits,
          typename Topology_,
          typename Geometry_,
          typename Function_,
          typename Gradient_>
class Isosurfacing_domain_3
{
public:
  using Geom_traits = GeomTraits;

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Topology = Topology_;
  using Vertex_descriptor = typename Topology_::Vertex_descriptor;
  using Edge_descriptor = typename Topology_::Edge_descriptor;
  using Cell_descriptor = typename Topology_::Cell_descriptor;

  using Vertices_incident_to_edge = typename Topology_::Vertices_incident_to_edge;
  using Cells_incident_to_edge = typename Topology_::Cells_incident_to_edge;
  using Cell_vertices = typename Topology_::Cell_vertices;
  using Cell_edges = typename Topology_::Cell_edges;

  using Geometry = Geometry_;
  using Function = Function_;
  using Gradient = Gradient_;

  static constexpr Cell_type CELL_TYPE = Topology_::CELL_TYPE;
  static constexpr std::size_t VERTICES_PER_CELL = Topology_::VERTICES_PER_CELL;
  static constexpr std::size_t EDGES_PER_CELL = Topology_::EDGES_PER_CELL;

private:
  const Topology m_topo;
  const Geometry m_geom;
  const Function m_func;
  const Gradient m_grad;
  const Geom_traits m_gt;

public:
  // creates a base domain from topology, geometry, implicit function, gradient, and geometric traits
  Isosurfacing_domain_3(const Topology& topo,
                        const Geometry& geom,
                        const Function& func,
                        const Gradient& grad,
                        const Geom_traits& gt)
    : m_topo{topo},
      m_geom{geom},
      m_func{func},
      m_grad{grad},
      m_gt(gt)
  { }

  const Geom_traits& geom_traits() const
  {
    return m_gt;
  }

public:
  // gets the position of vertex `v`
  decltype(auto) /*Point_3*/ point(const Vertex_descriptor& v) const
  {
    return m_geom.operator()(v);
  }

  // gets the value of the function at vertex `v`
  decltype(auto) /*FT*/ value(const Vertex_descriptor& v) const
  {
    return m_func.operator()(v);
  }

  // gets the gradient at vertex `v`
  decltype(auto) /*Vector_3*/ gradient(const Point_3& p) const
  {
    return m_grad.operator()(p);
  }

  // gets a container with the two vertices incident to the edge `e`
  decltype(auto) /*Vertices_incident_to_edge*/ incident_vertices(const Edge_descriptor& e) const
  {
    return m_topo.incident_vertices(e);
  }

  // gets a container with all cells incident to the edge `e`
  decltype(auto) /*Cells_incident_to_edge*/ incident_cells(const Edge_descriptor& e) const
  {
    return m_topo.incident_cells(e);
  }

  // gets a container with all vertices of the cell `c`
  decltype(auto) /*Cell_vertices*/ cell_vertices(const Cell_descriptor& c) const
  {
    return m_topo.cell_vertices(c);
  }

  // gets a container with all edges of the cell `c`
  decltype(auto) /*Cell_edges*/ cell_edges(const Cell_descriptor& c) const
  {
    return m_topo.cell_edges(c);
  }

  // iterates over all vertices `v`, calling `f(v)` on each of them
  template <typename ConcurrencyTag, typename Functor>
  void for_each_vertex(Functor& f) const
  {
    m_topo.for_each_vertex(f, ConcurrencyTag{});
  }

  // iterates over all edges `e`, calling `f(e)` on each of them
  template <typename ConcurrencyTag, typename Functor>
  void for_each_edge(Functor& f) const
  {
    m_topo.for_each_edge(f, ConcurrencyTag{});
  }

  // iterates over all cells `c`, calling `f(c)` on each of them
  template <typename ConcurrencyTag, typename Functor>
  void for_each_cell(Functor& f) const
  {
    m_topo.for_each_cell(f, ConcurrencyTag{});
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_ISOSURFACING_DOMAIN_3_H
