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
//                 Sven Oesau

#ifndef CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H
#define CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Isosurfacing_3/internal/Cell_type.h>
#include <CGAL/Isosurfacing_3/internal/tables.h>

#include <CGAL/tags.h>
#include <boost/functional/hash.hpp>
#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif // CGAL_LINKED_WITH_TBB

#include <array>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// this is to be able to specialize std::hash
struct OW_Edge_handle : public std::tuple<std::size_t, std::size_t>
{
  using std::tuple<std::size_t, std::size_t>::tuple; // inherit constructors
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

namespace std {

template <>
struct hash<CGAL::Isosurfacing::internal::OW_Edge_handle>
{
  std::size_t operator()(const CGAL::Isosurfacing::internal::OW_Edge_handle& e) const
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, std::get<0>(e));
    boost::hash_combine(seed, std::get<1>(e));
    return seed;
  }
};

} // namespace std

namespace CGAL {
namespace Isosurfacing {

template <typename Partition>
struct partition_traits;

template<typename GeomTraits>
struct partition_traits<CGAL::Octree<GeomTraits, std::vector<typename GeomTraits::Point_3> > >
{
  using Orthtree = CGAL::Octree<GeomTraits, std::vector<typename GeomTraits::Point_3> >;

public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Node_index = typename Orthtree::Node_index;
  using Uniform_coords = typename Orthtree::Global_coordinates;  // coordinates on max depth level

  using vertex_descriptor = std::size_t;
  using edge_descriptor = internal::OW_Edge_handle;
  using cell_descriptor = std::size_t;

  static constexpr Cell_type CELL_TYPE = CUBICAL_CELL;
  static constexpr std::size_t VERTICES_PER_CELL = 8;
  static constexpr std::size_t EDGES_PER_CELL = 12;

  using Edge_vertices = std::array<vertex_descriptor, 2>;
  using Cells_incident_to_edge = std::array<cell_descriptor, 4>;  // @todo: not always 4
  using Cell_vertices = std::array<vertex_descriptor, 8>;
  using Cell_edges = std::array<edge_descriptor, 12>;

public:
  static std::set<edge_descriptor> get_leaf_edges(const Orthtree& o) {
    std::set<edge_descriptor> leaf_edge_set;
    std::size_t dim = std::size_t(1) << o.depth();
    for (Node_index node_index : o.traverse(CGAL::Orthtrees::Leaves_traversal<Orthtree>(o)))
    {
      const Uniform_coords& coords_uniform = uniform_coordinates(node_index, o);

      // write all leaf edges in a set
      const Uniform_coords& coords_global = o.global_coordinates(node_index);
      const std::size_t depth = o.depth(node_index);
      const std::size_t df = std::size_t(1) << (o.depth() - depth);
      for (const auto& edge_voxels : internal::Cube_table::edge_to_voxel_neighbor)
      {
        bool are_all_voxels_leafs = true;
        for (const auto& node_ijk : edge_voxels)
        {
          const std::size_t x = coords_uniform[0] + df * node_ijk[0];
          const std::size_t y = coords_uniform[1] + df * node_ijk[1];
          const std::size_t z = coords_uniform[2] + df * node_ijk[2];
          // check for overflow / ignore edges on boundary
          if (x >= dim || y >= dim || z >= dim)
          {
            are_all_voxels_leafs = false;
            break;
          }

          const Node_index n = get_node(x, y, z, o);
          if (o.depth(n) > depth)
          {
            are_all_voxels_leafs = false;
            break;
          }
        }

        if (are_all_voxels_leafs)
        {
          // add to leaf edge set
          std::size_t e_gl = e_glIndex(edge_voxels[0][3],
            coords_global[0], coords_global[1], coords_global[2],
            depth);
          leaf_edge_set.insert({ e_gl, depth });
        }
      }
    }

    return leaf_edge_set;
  }

  static std::set<vertex_descriptor> get_leaf_vertices(const Orthtree& o) {
    std::set<vertex_descriptor> leaf_vertex_set;
    std::size_t dim = std::size_t(1) << o.depth();
    for (Node_index node_index : o.traverse(CGAL::Orthtrees::Leaves_traversal<Orthtree>(o)))
    {
      const Uniform_coords& coords_uniform = uniform_coordinates(node_index, o);

      // write all vertices edges in a set
      const Uniform_coords& coords_global = o.global_coordinates(node_index);
      const std::size_t depth = o.depth(node_index);
      const std::size_t df = std::size_t(1) << (o.depth() - depth);

      for (int i = 0; i < internal::Cube_table::N_VERTICES; ++i)
      {
        Uniform_coords v_coords = coords_global;
        typename Orthtree::Local_coordinates local_coords(i);

        for (int j = 0; j < 3; ++j)
          v_coords[j] += std::size_t(local_coords[j]);

        for (int j = 0; j < 3; ++j)
          v_coords[j] *= static_cast<uint32_t>(df);

        const std::size_t lex = lex_index(v_coords[0], v_coords[1], v_coords[2], o.depth());
        leaf_vertex_set.insert(lex);
      }
    }

    return leaf_vertex_set;
  }

  static void get_leaves(const Orthtree& o, std::vector<cell_descriptor> &cells, std::vector<edge_descriptor>& edges, std::vector<vertex_descriptor>& vertices) {
    std::set<edge_descriptor> leaf_edge_set;
    std::set<vertex_descriptor> leaf_vertex_set;
    std::size_t dim = std::size_t(1) << o.depth();
    cells.clear();
    for (Node_index node_index : o.traverse(CGAL::Orthtrees::Leaves_traversal<Orthtree>(o)))
    {
      const Uniform_coords& coords_uniform = uniform_coordinates(node_index, o);
      cells.push_back(node_index);

      // write all leaf edges in a set
      const Uniform_coords& coords_global = o.global_coordinates(node_index);
      const std::size_t depth = o.depth(node_index);
      const std::size_t df = std::size_t(1) << (o.depth() - depth);

      for (int i = 0; i < internal::Cube_table::N_VERTICES; ++i)
      {
        Uniform_coords v_coords = coords_global;
        typename Orthtree::Local_coordinates local_coords(i);

        for (int j = 0; j < 3; ++j)
          v_coords[j] += std::size_t(local_coords[j]);

        for (int j = 0; j < 3; ++j)
          v_coords[j] *= static_cast<uint32_t>(df);

        const std::size_t lex = lex_index(v_coords[0], v_coords[1], v_coords[2], o.depth());
        leaf_vertex_set.insert(lex);
      }

      for (const auto& edge_voxels : internal::Cube_table::edge_to_voxel_neighbor)
      {
        bool are_all_voxels_leafs = true;
        for (const auto& node_ijk : edge_voxels)
        {
          const std::size_t x = coords_uniform[0] + df * node_ijk[0];
          const std::size_t y = coords_uniform[1] + df * node_ijk[1];
          const std::size_t z = coords_uniform[2] + df * node_ijk[2];
          // check for overflow / ignore edges on boundary
          if (x >= dim || y >= dim || z >= dim)
          {
            are_all_voxels_leafs = false;
            break;
          }

          const Node_index n = get_node(x, y, z, o);
          if (o.depth(n) > depth)
          {
            are_all_voxels_leafs = false;
            break;
          }
        }

        if (are_all_voxels_leafs)
        {
          // add to leaf edge set
          std::size_t e_gl = e_glIndex(edge_voxels[0][3],
            coords_global[0], coords_global[1], coords_global[2],
            depth);
          leaf_edge_set.insert({ e_gl, depth });
        }
      }
    }

    edges.clear();
    vertices.clear();

    std::copy(leaf_edge_set.begin(), leaf_edge_set.end(), std::back_inserter(edges));
    std::copy(leaf_vertex_set.begin(), leaf_vertex_set.end(), std::back_inserter(vertices));
  }

  static Point_3 point(const vertex_descriptor& v,
                                          const Orthtree& o)
  {
    std::size_t dim_ = std::size_t(1) << o.depth();
    const auto bbox = o.bbox(0);
    FT hx_ = (bbox.xmax() - bbox.xmin()) / dim_;
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(v, o.depth());

    const FT x0 = bbox.xmin() + i * hx_;
    const FT y0 = bbox.ymin() + j * hx_;
    const FT z0 = bbox.zmin() + k * hx_;

    return { x0, y0, z0 };
  }


  static Edge_vertices incident_vertices(const edge_descriptor& e,
                                         const Orthtree& o)
  {
    return edge_vertices(e, o);
  }

  static Cells_incident_to_edge incident_cells(const edge_descriptor& e_id,
    const Orthtree& o)
  {
    namespace Tables = internal::Cube_table;

    std::size_t e_global_id, depth;
    std::tie(e_global_id, depth) = static_cast<const std::tuple<std::size_t, std::size_t>&>(e_id);
    const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];

    const std::size_t max_depth = o.depth();
    const std::size_t df = std::size_t(1) << (max_depth - depth);

    const size_t v0_lex_index = e_global_id / 3;
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(v0_lex_index, depth);
    i *= df;
    j *= df;
    k *= df;

    const auto& voxel_neighbors = Tables::edge_to_voxel_neighbor[e_local_index];
    Node_index n0 = get_node(i + voxel_neighbors[0][0], j + voxel_neighbors[0][1], k + voxel_neighbors[0][2], o);
    Node_index n1 = get_node(i + voxel_neighbors[1][0], j + voxel_neighbors[1][1], k + voxel_neighbors[1][2], o);
    Node_index n2 = get_node(i + voxel_neighbors[2][0], j + voxel_neighbors[2][1], k + voxel_neighbors[2][2], o);
    Node_index n3 = get_node(i + voxel_neighbors[3][0], j + voxel_neighbors[3][1], k + voxel_neighbors[3][2], o);

    const Uniform_coords n0_uniform_coords = uniform_coordinates(n0, o);
    const Uniform_coords n1_uniform_coords = uniform_coordinates(n1, o);
    const Uniform_coords n2_uniform_coords = uniform_coordinates(n2, o);
    const Uniform_coords n3_uniform_coords = uniform_coordinates(n3, o);

    std::size_t n0_lex = lex_index(n0_uniform_coords[0], n0_uniform_coords[1], n0_uniform_coords[2], max_depth);
    std::size_t n1_lex = lex_index(n1_uniform_coords[0], n1_uniform_coords[1], n1_uniform_coords[2], max_depth);
    std::size_t n2_lex = lex_index(n2_uniform_coords[0], n2_uniform_coords[1], n2_uniform_coords[2], max_depth);
    std::size_t n3_lex = lex_index(n3_uniform_coords[0], n3_uniform_coords[1], n3_uniform_coords[2], max_depth);

    return { n0_lex, n1_lex, n2_lex, n3_lex };
  }

  static std::size_t depth(const cell_descriptor& c,
    const Orthtree& o) {
    return o.depth(c);
  }

  static Cell_vertices cell_vertices(const cell_descriptor& c,
                                     const Orthtree& o)
  {
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(c, o.depth());
    Node_index node_index = get_node(i, j, k, o);
    const std::size_t df = std::size_t(1) << (o.depth() - o.depth(node_index));

    std::array<vertex_descriptor, 8> v;
    for (int v_id = 0; v_id < internal::Cube_table::N_VERTICES; ++v_id)
    {
      const int* l = internal::Cube_table::local_vertex_location[v_id];
      const std::size_t lex = lex_index(i + df * l[0], j + df * l[1], k + df * l[2], o.depth());
      v[v_id] = lex;
    }

    return v;
  }

  static Cell_edges cell_edges(const cell_descriptor& c,
                               const Orthtree& o)
  {
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(c, o.depth());
    Node_index node_index = get_node(i, j, k, o);

    const Uniform_coords& coords_global = o.global_coordinates(node_index);
    const std::size_t depth = o.depth(node_index);

    std::array<edge_descriptor, internal::Cube_table::N_EDGES> edges;
    for (std::size_t e_id = 0; e_id < edges.size(); ++e_id)
    {
      const std::size_t e_gl = e_glIndex(e_id, coords_global[0], coords_global[1], coords_global[2], depth);
      edges[e_id] = { e_gl, depth };
    }

    return edges;
  }

  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              const Orthtree& o,
                              const CGAL::Sequential_tag)
  {
    for(const vertex_descriptor& v : get_leaf_vertices(o))
      f(v);
  }

  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              std::vector<vertex_descriptor>& vertices,
                              const Orthtree& o,
                              const CGAL::Sequential_tag)
  {
    for (const vertex_descriptor& v : vertices)
      f(v);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Orthtree& o,
                            Sequential_tag)
  {
    for(const edge_descriptor& e : get_leaf_edges(o))
       f(e);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            std::vector<edge_descriptor>& edges,
                            const Orthtree& o,
                            Sequential_tag)
  {
    if (edges.empty()) {
      std::set<edge_descriptor> edge_set = get_leaf_edges(o, edges);
      std::copy(edge_set.begin(), edge_set.end(), std::back_inserter(edges));
    }

    for (const edge_descriptor& e : edges)
      f(e);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            const Orthtree& o,
                            CGAL::Sequential_tag)
  {
    for (const cell_descriptor& v : o.traverse(CGAL::Orthtrees::Leaves_traversal<typename Orthtree::Orthtree>(o)))
      f(v);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            std::vector<cell_descriptor>& cells,
                            const Orthtree& o,
                            CGAL::Sequential_tag)
  {
    if (cells.empty()) {
      auto cell_range = o.traverse(CGAL::Orthtrees::Leaves_traversal<typename Orthtree::Orthtree>(o));
      std::copy(cell_range.begin(), cell_range.end(), std::back_inserter(cells));
    }

    for (const cell_descriptor& v : cells)
      f(v);
  }

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              const Orthtree& o,
                              Parallel_tag)
  {
    auto edges = get_leaf_edges(o);

    tbb::parallel_for_each(edges.begin(), edges.end(), f);
  }

  template <typename Functor>
  static void for_each_vertex(Functor& f,
                              std::vector<vertex_descriptor>& vertices,
                              const Orthtree& o,
                              Parallel_tag)
  {
    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
    {
      for(std::size_t i = r.begin(); i != r.end(); ++i)
        f(vertices[i]);
    };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, vertices.size()), iterator);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            const Orthtree& o,
                            Parallel_tag)
  {
    std::set<edge_descriptor> edges = get_leaf_edges(o);

    tbb::parallel_for_each(edges.begin(), edges.end(), f);
  }

  template <typename Functor>
  static void for_each_edge(Functor& f,
                            std::vector<edge_descriptor>& edges,
                            const Orthtree& o,
                            Parallel_tag)
  {
    if (edges.empty()) {
      std::set<edge_descriptor> edge_set = get_leaf_edges(o);
      std::copy(edge_set.begin(), edge_set.end(), std::back_inserter(edges));
    }

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          f(edges[i]);
      };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, edges.size()), iterator);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            const Orthtree& o,
                            Parallel_tag)
  {
    auto func = [&](const typename Orthtree::Node_index& n)
    {
      const Uniform_coords& uc = uniform_coordinates(n, o);
      f(lex_index(uc[0], uc[1], uc[2], o.depth()));
    };

    const auto& cells = o.traverse(CGAL::Orthtrees::Leaves_traversal<Orthtree>(o));
    tbb::parallel_for_each(cells.begin(), cells.end(), func);
  }

  template <typename Functor>
  static void for_each_cell(Functor& f,
                            std::vector<cell_descriptor>& cells,
                            const Orthtree& o,
                            Parallel_tag)
  {
    if (cells.empty()) {
      auto cell_range = o.traverse(CGAL::Orthtrees::Leaves_traversal<typename Orthtree::Orthtree>(o));
      std::copy(cell_range.begin(), cell_range.end(), std::back_inserter(cells));
    }

    auto iterator = [&](const tbb::blocked_range<std::size_t>& r)
      {
        for (std::size_t i = r.begin(); i != r.end(); ++i) {
          const Uniform_coords& uc = uniform_coordinates(cells[i], o);
          f(lex_index(uc[0], uc[1], uc[2], o.depth()));
        }
      };

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()), iterator);
  }
#endif // CGAL_LINKED_WITH_TBB

  template <typename ConcurrencyTag, typename Functor>
  static void for_each_vertex(Functor& f, const Orthtree& o) { return for_each_vertex(f, o, ConcurrencyTag{}); }
  template <typename ConcurrencyTag, typename Functor>
  static void for_each_edge(Functor& f, std::vector<edge_descriptor>& edges, const Orthtree& o) { return for_each_edge(f, edges, o, ConcurrencyTag{}); }
  template <typename ConcurrencyTag, typename Functor>
  static void for_each_edge(Functor& f, const Orthtree& o) { return for_each_edge(f, o, ConcurrencyTag{}); }
  template <typename ConcurrencyTag, typename Functor>
  static void for_each_cell(Functor& f, std::vector<cell_descriptor>& cells, const Orthtree& o) { return for_each_cell(f, cells, o, ConcurrencyTag{}); }
  template <typename ConcurrencyTag, typename Functor>
  static void for_each_cell(Functor& f, const Orthtree& o) { return for_each_cell(f, o, ConcurrencyTag{}); }

private:
  static Uniform_coords uniform_coordinates(Node_index node_index, const Orthtree &o)
  {
    Uniform_coords coords = o.global_coordinates(node_index);
    const std::size_t df = std::size_t(1) << (o.depth() - o.depth(node_index));
    for (int i = 0; i < 3; ++i)
      coords[i] *= static_cast<uint32_t>(df);

    return coords;
  }

  static std::tuple<std::size_t, std::size_t, std::size_t> ijk_index(const std::size_t lex_index,
    const std::size_t depth)
  {
    const std::size_t dim = (std::size_t(1) << depth) + 1;
    return std::make_tuple(lex_index % dim, (lex_index / dim) % dim, lex_index / (dim * dim));
  }
  // computes unique edge global index.
// \param e local edge index
// \param i_idx i-index of cell
// \param j_idx j-index of cell
// \param k_idx k-index of cell
// \param depth of cell
  static std::size_t e_glIndex(const std::size_t e,
    const std::size_t i_idx,
    const std::size_t j_idx,
    const std::size_t k_idx,
    const std::size_t depth)
  {
    const unsigned long long gei_pattern_ = 670526590282893600ull;
    const size_t i = i_idx + (size_t)((gei_pattern_ >> 5 * e) & 1);        // global_edge_id[eg][0];
    const size_t j = j_idx + (size_t)((gei_pattern_ >> (5 * e + 1)) & 1);  // global_edge_id[eg][1];
    const size_t k = k_idx + (size_t)((gei_pattern_ >> (5 * e + 2)) & 1);  // global_edge_id[eg][2];
    const size_t offs = (size_t)((gei_pattern_ >> (5 * e + 3)) & 3);

    return (3 * lex_index(i, j, k, depth) + offs);
  }

  static std::size_t lex_index(const std::size_t i,
    const std::size_t j,
    const std::size_t k,
    const std::size_t depth)
  {
    std::size_t dim = (std::size_t(1) << depth) + 1;
    return k * dim * dim + j * dim + i;
  }

  static std::array<vertex_descriptor, 2> edge_vertices(const edge_descriptor& e_id, const Orthtree& o)
  {
    namespace Tables = internal::Cube_table;

    std::size_t e_global_id, depth, max_depth = o.depth();
    std::tie(e_global_id, depth) = static_cast<const std::tuple<std::size_t, std::size_t>&>(e_id);
    const std::size_t df = std::size_t(1) << (max_depth - depth);

    const size_t v0_lex_index = e_global_id / 3;
    std::size_t i0, j0, k0;
    std::tie(i0, j0, k0) = ijk_index(v0_lex_index, depth);

    // v1
    const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];
    const int* v1_local = Tables::local_vertex_location[Tables::edge_to_vertex[e_local_index][1]];

    const std::size_t i1 = i0 + v1_local[0];
    const std::size_t j1 = j0 + v1_local[1];
    const std::size_t k1 = k0 + v1_local[2];

    const std::size_t v0 = lex_index(df * i0, df * j0, df * k0, max_depth);
    const std::size_t v1 = lex_index(df * i1, df * j1, df * k1, max_depth);

    return { v0, v1 };
  }

  static Node_index get_node(const std::size_t i,
    const std::size_t j,
    const std::size_t k,
    const Orthtree& o)
  {
    Node_index node_index = o.root();
    const std::size_t x = i;
    const std::size_t y = j;
    const std::size_t z = k;

    while (!o.is_leaf(node_index))
    {
      std::size_t dist_to_max = o.depth() - o.depth(node_index) - 1;
      typename Orthtree::Local_coordinates loc;
      if (x & (std::size_t(1) << dist_to_max))
        loc[0] = true;

      if (y & (std::size_t(1) << dist_to_max))
        loc[1] = true;

      if (z & (std::size_t(1) << dist_to_max))
        loc[2] = true;

      node_index = o.child(node_index, loc.to_ulong());
    }

    return node_index;
  }

  Node_index get_node(const std::size_t lex_index, const Orthtree &o) const
  {
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(lex_index, o.depth());
    return get_node(i, j, k, o);
  }
};

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_PARTITION_TRAITS_OCTREE_H
