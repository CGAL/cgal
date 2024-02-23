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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_CARTESIAN_GRID_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
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

template <typename GeomTraits, typename MemoryPolicy>
class Cartesian_grid_3;

template <typename Partition>
struct partition_traits;

template <typename GeomTraits, typename MemoryPolicy>
struct partition_traits<Cartesian_grid_3<GeomTraits, MemoryPolicy> >
{
  using Self = Cartesian_grid_3<GeomTraits, MemoryPolicy>;

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

  static decltype(auto) /*Point_3*/ point(const Vertex_descriptor& v,
                                          const Grid& g)
  {
    return g.point(v[0], v[1], v[2]);
  }

  // returns a container with the two vertices incident to edge e
  static Vertices_incident_to_edge incident_vertices(const Edge_descriptor& e,
                                                     const Grid&)
  {
    Vertices_incident_to_edge ev;
    ev[0] = { e[0], e[1], e[2] };  // start vertex
    ev[1] = { e[0], e[1], e[2] };  // end vertex
    ev[1][e[3]] += 1;              // one position further in the direction of the edge
    return ev;
  }

  // returns a container with all cells incident to edge e
  static Cells_incident_to_edge incident_cells(const Edge_descriptor& e,
                                               const Grid&)
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

  // returns a container with all vertices of cell c
  static Cell_vertices cell_vertices(const Cell_descriptor& c,
                                     const Grid&)
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

  // returns a container with all edges of cell c
  static Cell_edges cell_edges(const Cell_descriptor& c,
                               const Grid&)
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
  static void for_each_vertex(Functor& f,
                              const Grid& g,
                              const CGAL::Sequential_tag)
  {
    for(std::size_t i=0; i<g.xdim(); ++i)
      for(std::size_t j=0; j<g.ydim(); ++j)
        for(std::size_t k=0; k<g.zdim(); ++k)
          f({i, j, k});
  }

  // iterates sequentially over all edges e calling f(e) on every one
  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Grid& g,
                            const CGAL::Sequential_tag)
  {
    for(std::size_t i=0; i<g.xdim()-1; ++i) {
      for(std::size_t j=0; j<g.ydim()-1; ++j) {
        for(std::size_t k=0; k<g.zdim()-1; ++k)
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
  static void for_each_cell(Functor& f,
                            const Grid& g,
                            const CGAL::Sequential_tag)
  {
    for(std::size_t i=0; i<g.xdim()-1; ++i)
      for(std::size_t j=0; j<g.ydim()-1; ++j)
        for(std::size_t k=0; k<g.zdim()-1; ++k)
          f({i, j, k});
  }

  #ifdef CGAL_LINKED_WITH_TBB
  // iterates in parallel over all vertices v calling f(v) on every one
  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              const Grid& g,
                              const CGAL::Parallel_tag)
  {
    const std::size_t sj = g.ydim();
    const std::size_t sk = g.zdim();

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i=r.begin(); i!=r.end(); ++i)
        for(std::size_t j=0; j<sj; ++j)
          for(std::size_t k=0; k<sk; ++k)
            f({i, j, k});
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, g.xdim()), iterator);
  }

  // iterates in parallel over all edges e calling f(e) on every one
  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Grid& g,
                            const CGAL::Parallel_tag)
  {
    const std::size_t sj = g.ydim();
    const std::size_t sk = g.zdim();

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i=r.begin(); i != r.end(); ++i) {
        for(std::size_t j=0; j<sj-1; ++j) {
          for(std::size_t k=0; k<sk-1; ++k)
          {
            f({i, j, k, 0});
            f({i, j, k, 1});
            f({i, j, k, 2});
          }
        }
      }
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, g.xdim() - 1), iterator);
  }

  // iterates in parallel over all cells c calling f(c) on every one
  template <typename Functor>
  static void for_each_cell(Functor& f,
                            const Grid& g,
                            const CGAL::Parallel_tag)
  {
    const std::size_t sj = g.ydim();
    const std::size_t sk = g.zdim();

    // for now only parallelize outer loop
    auto iterator = [&f, sj, sk](const tbb::blocked_range3d<std::size_t>& r)
    {
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

    tbb::blocked_range3d<std::size_t> range(0, g.xdim() - 1, 0, g.ydim() - 1, 0, g.zdim() - 1);
    tbb::parallel_for(range, iterator);
  }
  #endif // CGAL_LINKED_WITH_TBB
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_CARTESIAN_GRID_3_H
