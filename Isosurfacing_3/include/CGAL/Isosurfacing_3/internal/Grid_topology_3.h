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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_GRID_TOPOLOGY_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_GRID_TOPOLOGY_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/tables.h>
#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#endif // CGAL_LINKED_WITH_TBB

#include <array>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// The topology of a Cartesian grid.
// All elements are created with the help of the cube tables.
class Grid_topology_3
{
public:
  // identifies a vertex by its (i, j, k) indices
  using Vertex_descriptor = std::array<std::size_t, 3>;

  // identifies an edge by its starting vertex (i, j, k) and the direction x -> 0, y -> 1, z -> 2
  using Edge_descriptor = std::array<std::size_t, 4>;

  // identifies a cell by its corner vertex with the smallest (i, j, k) index
  using Cell_descriptor = std::array<std::size_t, 3>;

  static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
  static constexpr std::size_t VERTICES_PER_CELL = 8;
  static constexpr std::size_t EDGES_PER_CELL = 12;

  using Vertices_incident_to_edge = std::array<Vertex_descriptor, 2>;
  using Cells_incident_to_edge = std::array<Cell_descriptor, 4>;
  using Cell_vertices = std::array<Vertex_descriptor, VERTICES_PER_CELL>;
  using Cell_edges = std::array<Edge_descriptor, EDGES_PER_CELL>;

public:
  // creates the topology of a grid with size_i, size_j, and size_k vertices in the respective dimensions.
  Grid_topology_3(const std::size_t size_i,
                  const std::size_t size_j,
                  const std::size_t size_k)
    : size_i{size_i},
      size_j{size_j},
      size_k{size_k}
  { }

  // gets a container with the two vertices incident to edge e
  Vertices_incident_to_edge incident_vertices(const Edge_descriptor& e) const
  {
    Vertices_incident_to_edge ev;
    ev[0] = { e[0], e[1], e[2] };  // start vertex
    ev[1] = { e[0], e[1], e[2] };  // end vertex
    ev[1][e[3]] += 1;              // one position further in the direction of the edge
    return ev;
  }

  // gets a container with all cells incident to edge e
  Cells_incident_to_edge incident_cells(const Edge_descriptor& e) const
  {
    // lookup the neighbor cells relative to the edge
    const int local = internal::Cube_table::edge_store_index[e[3]];
    auto neighbors = internal::Cube_table::edge_to_voxel_neighbor[local];

    Cells_incident_to_edge cite;
    for(std::size_t i=0; i<cite.size(); ++i) {
      for(std::size_t j=0; j<cite[i].size(); ++j)
      {
        // offset the relative indices by the edge position
        cite[i][j] = e[j] + neighbors[i][j];
      }
    }

    return cite;
  }

  // gets a container with all vertices of cell c
  Cell_vertices cell_vertices(const Cell_descriptor& c) const
  {
    Cell_vertices cv;
    for(std::size_t i=0; i<cv.size(); ++i) {
      for(std::size_t j=0; j<c.size(); ++j)
      {
        // lookup the relative vertex indices and offset them by the cell position
        cv[i][j] = c[j] + internal::Cube_table::local_vertex_position[i][j];
      }
    }

    return cv;
  }

  // gets a container with all edges of cell c
  Cell_edges cell_edges(const Cell_descriptor& c) const
  {
    Cell_edges ce;
    for(std::size_t i=0; i<ce.size(); ++i) {
      for(std::size_t j=0; j<c.size(); ++j)
      {
        // lookup the relative edge indices and offset them by the cell position
        ce[i][j] = c[j] + internal::Cube_table::global_edge_id[i][j];
      }

      // set the direction of the edge
      ce[i][3] = internal::Cube_table::global_edge_id[i][3];
    }

    return ce;
  }

  // iterates sequentially over all vertices v calling f(v) on every one
  template <typename Functor>
  void for_each_vertex(Functor& f, Sequential_tag) const
  {
    for(std::size_t i=0; i<size_i; ++i)
      for(std::size_t j=0; j<size_j; ++j)
        for(std::size_t k=0; k<size_k; ++k)
          f({i, j, k});
  }

  // iterates sequentially over all edges e calling f(e) on every one
  template <typename Functor>
  void for_each_edge(Functor& f, Sequential_tag) const
  {
    for(std::size_t i=0; i<size_i-1; ++i) {
      for(std::size_t j=0; j<size_j-1; ++j) {
        for(std::size_t k=0; k<size_k-1; ++k)
        {
          // all three edges starting at vertex (i, j, k)
          f({i, j, k, 0});
          f({i, j, k, 1});
          f({i, j, k, 2});
        }
      }
    }
  }

  // iterates sequentially over all cells c calling f(c) on every one
  template <typename Functor>
  void for_each_cell(Functor& f, Sequential_tag) const
  {
    for(std::size_t i=0; i<size_i-1; ++i)
      for(std::size_t j=0; j<size_j-1; ++j)
        for(std::size_t k=0; k<size_k-1; ++k)
          f({i, j, k});
  }

#ifdef CGAL_LINKED_WITH_TBB
  // iterates in parallel over all vertices v calling f(v) on every one
  template <typename Functor>
  void for_each_vertex(Functor& f, Parallel_tag) const
  {
    const std::size_t sj = size_j;
    const std::size_t sk = size_k;

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        for(std::size_t j=0; j<sj; ++j)
          for(std::size_t k=0; k<sk; ++k)
            f({i, j, k});
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_i), iterator);
  }

  // iterates in parallel over all edges e calling f(e) on every one
  template <typename Functor>
  void for_each_edge(Functor& f, Parallel_tag) const
  {
    const std::size_t sj = size_j;
    const std::size_t sk = size_k;

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i) {
        for(std::size_t j=0; j<sj - 1; ++j) {
          for(std::size_t k=0; k<sk - 1; ++k)
          {
            f({i, j, k, 0});
            f({i, j, k, 1});
            f({i, j, k, 2});
          }
        }
      }
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size_i - 1), iterator);
  }

  // iterates in parallel over all cells c calling f(c) on every one
  template <typename Functor>
  void for_each_cell(Functor& f, Parallel_tag) const
  {
    const std::size_t sj = size_j;
    const std::size_t sk = size_k;

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range3d<std::size_t>& r) {
      const std::size_t i_begin = r.pages().begin();
      const std::size_t i_end = r.pages().end();
      const std::size_t j_begin = r.rows().begin();
      const std::size_t j_end = r.rows().end();
      const std::size_t k_begin = r.cols().begin();
      const std::size_t k_end = r.cols().end();

      for(std::size_t i = i_begin; i != i_end; ++i)
        for(std::size_t j = j_begin; j != j_end; ++j)
          for(std::size_t k = k_begin; k != k_end; ++k)
            f({i, j, k});
    };

    tbb::blocked_range3d<std::size_t> range(0, size_i - 1, 0, size_j - 1, 0, size_k - 1);
    tbb::parallel_for(range, iterator);
  }
#endif // CGAL_LINKED_WITH_TBB

private:
  std::size_t size_i;
  std::size_t size_j;
  std::size_t size_k;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_GRID_TOPOLOGY_3_H
