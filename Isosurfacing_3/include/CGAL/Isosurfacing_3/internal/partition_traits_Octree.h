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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H
#define CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>
#include <CGAL/Isosurfacing_3/internal/Cell_type.h>

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif // CGAL_LINKED_WITH_TBB

#include <array>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename GeomTraits>
class Octree_wrapper;

} // namespace internal

template <typename Partition>
struct partition_traits;

template <typename GeomTraits>
struct partition_traits<internal::Octree_wrapper<GeomTraits> >
{
  using Octree = internal::Octree_wrapper<GeomTraits>;

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
  static decltype(auto) /*Point_3*/ point(const Vertex_descriptor& v,
                                          const Octree& o)
  {
    return o.point(v);
  }

  static Vertices_incident_to_edge incident_vertices(const Edge_descriptor& e,
                                                     const Octree& o)
  {
    return o.edge_vertices(e);
  }

  static Cells_incident_to_edge incident_cells(const Edge_descriptor& e,
                                               const Octree& o)
  {
    return o.edge_voxels(e);
  }

  static Cell_vertices cell_vertices(const Cell_descriptor& c,
                                     const Octree& o)
  {
    return o.voxel_vertices(c);
  }

  static Cell_edges cell_edges(const Cell_descriptor& c,
                               const Octree& o)
  {
    return o.voxel_edges(c);
  }

  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              const Octree& o,
                              const CGAL::Sequential_tag)
  {
    for(const Vertex_descriptor& v : o.leaf_vertices())
      f(v);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Octree& o,
                            Sequential_tag)
  {
    for(const Edge_descriptor& e : o.leaf_edges())
      f(e);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            const Octree& o,
                            CGAL::Sequential_tag)
  {
    for(const Cell_descriptor& v : o.leaf_voxels())
      f(v);
  }

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              const Octree& o,
                              Parallel_tag)
  {
    const auto& vertices = o.leaf_vertices();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(vertices[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Octree& o,
                            Parallel_tag)
  {
    const auto& edges = o.leaf_edges();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(edges[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            const Octree& o,
                            Parallel_tag)
  {
    const auto& cells = o.leaf_voxels();

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(cells[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
  }
#endif // CGAL_LINKED_WITH_TBB
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H
