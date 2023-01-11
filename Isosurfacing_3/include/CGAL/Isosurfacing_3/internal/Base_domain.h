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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_BASE_DOMAIN_H
#define CGAL_ISOSURFACING_3_INTERNAL_BASE_DOMAIN_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>

#include <memory>

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
class Base_domain
{
public:
  using Geom_traits = GeomTraits;

  using FT = typename Geom_traits::FT;
  using Point = typename Geom_traits::Point_3;
  using Vector = typename Geom_traits::Vector_3;

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

public:
  // creates a base_domain from a topology, geometry, input function, and gradient
  Base_domain(const Topology& topo,
              const Geometry& geom,
              const Function& func,
              const Gradient& grad)
      : topo(topo),
        geom(geom),
        func(func),
        grad(grad)
  { }

  // gets the position of vertex `v`
  Point position(const Vertex_descriptor& v) const
  {
    return geom.operator()(v);
  }

  // gets the value of the function at vertex `v`
  FT value(const Vertex_descriptor& v) const
  {
    return func.operator()(v);
  }

  // gets the gradient at vertex `v`
  Vector gradient(const Point& p) const
  {
    return grad.operator()(p);
  }

  // gets a container with the two vertices incident to the edge `e`
  Vertices_incident_to_edge edge_vertices(const Edge_descriptor& e) const
  {
    return topo.edge_vertices(e);
  }

  // gets a container with all cells incident to the edge `e`
  Cells_incident_to_edge cells_incident_to_edge(const Edge_descriptor& e) const
  {
    return topo.cells_incident_to_edge(e);
  }

  // gets a container with all vertices of the cell `c`
  Cell_vertices cell_vertices(const Cell_descriptor& c) const
  {
    return topo.cell_vertices(c);
  }

  // gets a container with all edges of the cell `c`
  Cell_edges cell_edges(const Cell_descriptor& c) const
  {
    return topo.cell_edges(c);
  }

  // iterates over all vertices `v`, calling `f(v)` on each of them
  template <typename Concurrency_tag, typename Functor>
  void iterate_vertices(Functor& f) const
  {
    topo.iterate_vertices(f, Concurrency_tag());
  }

  // iterates over all edges `e`, calling `f(e)` on each of them
  template <typename Concurrency_tag, typename Functor>
  void iterate_edges(Functor& f) const
  {
    topo.iterate_edges(f, Concurrency_tag());
  }

  // iterates over all cells `c`, calling `f(c)` on each of them
  template <typename Concurrency_tag, typename Functor>
  void iterate_cells(Functor& f) const
  {
    topo.iterate_cells(f, Concurrency_tag());
  }

private:
  const Topology topo;
  const Geometry geom;
  const Function func;
  const Gradient grad;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_BASE_DOMAIN_H
