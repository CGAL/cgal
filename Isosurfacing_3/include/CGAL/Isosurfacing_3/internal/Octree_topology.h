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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H
#define CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Octree>
class Octree_topology
{
public:
  using Vertex_descriptor = typename Octree::Vertex_handle;
  using Edge_descriptor = typename Octree::Edge_handle;
  using Cell_descriptor = typename Octree::Voxel_handle;

  static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
  static constexpr std::size_t VERTICES_PER_CELL = 8;
  static constexpr std::size_t EDGES_PER_CELL = 12;

  using Vertices_incident_to_edge = std::array<Vertex_descriptor, 2>;
  using Cells_incident_to_edge = std::array<Cell_descriptor, 4>;  // @todo: not always 4
  using Cell_vertices = std::array<Vertex_descriptor, 8>;
  using Cell_edges = std::array<Edge_descriptor, 12>;

public:
  Octree_topology(const Octree& octree)
    : m_octree(octree)
  { }

  Vertices_incident_to_edge incident_vertices(const Edge_descriptor& e) const
  {
    return m_octree.edge_vertices(e);
  }

  Cells_incident_to_edge incident_cells(const Edge_descriptor& e) const
  {
    return m_octree.edge_voxels(e);
  }

  Cell_vertices cell_vertices(const Cell_descriptor& c) const
  {
    return m_octree.voxel_vertices(c);
  }

  Cell_edges cell_edges(const Cell_descriptor& c) const
  {
    return m_octree.voxel_edges(c);
  }

  template <typename Functor>
  void for_each_vertex(Functor& f, Sequential_tag) const
  {
    for(const Vertex_descriptor& v : m_octree.leaf_vertices())
      f(v);
  }

  template <typename Functor>
  void for_each_edge(Functor& f, Sequential_tag) const
  {
    for(const Edge_descriptor& e : m_octree.leaf_edges())
      f(e);
  }

  template <typename Functor>
  void for_each_cell(Functor& f, Sequential_tag) const
  {
    for(const Cell_descriptor& v : m_octree.leaf_voxels())
      f(v);
  }

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Functor>
  void for_each_vertex(Functor& f, Parallel_tag) const
  {
    const auto& vertices = m_octree.leaf_vertices();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(vertices[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
  }

  template <typename Functor>
  void for_each_edge(Functor& f, Parallel_tag) const
  {
    const auto& edges = m_octree.leaf_edges();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(edges[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
  }

  template <typename Functor>
  void for_each_cell(Functor& f, Parallel_tag) const
  {
    const auto& cells = m_octree.leaf_voxels();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(cells[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
  }
#endif // CGAL_LINKED_WITH_TBB

private:
  const Octree& m_octree;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_OCTREE_TOPOLOGY_H
