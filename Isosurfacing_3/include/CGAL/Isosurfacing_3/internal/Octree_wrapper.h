// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France), GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Daniel Zint
//                 Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_OCTREE_WRAPPER_H
#define CGAL_ISOSURFACING_3_INTERNAL_OCTREE_WRAPPER_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/partition_traits_Octree.h>
#include <CGAL/Isosurfacing_3/internal/tables.h>

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>

#include <boost/functional/hash.hpp>

#include <array>
#include <map>
#include <tuple>
#include <vector>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename GeomTraits>
class Octree_wrapper
{
  /*
    * Naming convention from "A parallel dual marching cubes approach to quad
    * only surface reconstruction - Grosso & Zint"
    *
    *        ^ y
    *        |
    *       v2------e2------v3
    *       /|             /|
    *     e11|           e10|
    *     /  e3          /  e1
    *   v6------e6------v7  |
    *    |   |          |   |
    *    |  v0------e0--|---v1 --> x
    *    e7 /           e5 /
    *    | e8           | e9
    *    |/             |/
    *   v4------e4------v5
    *   /
    *  < z
    */

public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Octree = CGAL::Octree<Geom_traits, std::vector<Point_3> >;

  using Vertex_handle = std::size_t;
  using Voxel_handle = std::size_t;

  using Node_index = typename Octree::Node_index;
  using Uniform_coords = typename Octree::Global_coordinates;  // coordinates on max depth level

private:
  std::size_t max_depth_ = 0;

  FT offset_x_;
  FT offset_y_;
  FT offset_z_;

  // @todo should be iso_cuboid_3 if we mirrored the other domains, but we'll see
  // once CGAL::Octree has been updated and this wrapper isn't required anymore.
  CGAL::Bbox_3 bbox_;

  std::size_t dim_ = 1;

  FT hx_ = 0;

  std::vector<Point_3> point_range_;
  Octree octree_;
  GeomTraits gt_;

  // std::set<Uniform_coords> leaf_node_uniform_coordinates_;
  std::vector<Voxel_handle> leaf_voxels_;
  std::vector<Edge_handle> leaf_edges_;
  std::vector<Vertex_handle> leaf_vertices_;

public:
  Octree_wrapper(const CGAL::Bbox_3& bbox)
    : offset_x_(bbox.xmin()),
      offset_y_(bbox.ymin()),
      offset_z_(bbox.zmin()),
      bbox_(bbox),
      point_range_({{bbox.xmin(), bbox.ymin(), bbox.zmin()},
                    {bbox.xmax(), bbox.ymax(), bbox.zmax()}}),
      octree_(point_range_)
  { }

  template <typename Split_predicate>
  void refine(const Split_predicate& split_predicate)
  {
    namespace Tables = internal::Cube_table;

    octree_.refine(split_predicate);

    max_depth_ = octree_.depth();
    dim_ = std::size_t(1) << max_depth_;
    hx_ = bbox_.x_span() / dim_;

    // store leaf elements in sets + initialize value maps
    /*std::set<Voxel_handle> leaf_voxels_set;
    std::set<Edge_handle> leaf_edges_set;
    std::set<Vertex_handle> leaf_vertices_set;
    std::size_t cell_count = 0;
    for(Node_index node_index : octree_.traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(octree_)))
    {
      const Uniform_coords& coords_uniform = uniform_coordinates(node_index);

      // write all leaf nodes in a set
      leaf_voxels_set.insert(lex_index(coords_uniform[0], coords_uniform[1], coords_uniform[2], max_depth_));
      cell_count++;

      // init vertex values
      for(int i=0; i<Tables::N_VERTICES; ++i)
      {
        Uniform_coords vuc = vertex_uniform_coordinates(node_index, i);
        const std::size_t lex = lex_index(vuc[0], vuc[1], vuc[2], max_depth_);
        leaf_vertices_set.insert(lex);
      }

      // write all leaf edges in a set
      const Uniform_coords& coords_global = octree_.global_coordinates(node_index);
      const std::size_t depth = octree_.depth(node_index);
      const std::size_t df = depth_factor(depth);
      for(const auto& edge_voxels : Tables::edge_to_voxel_neighbor)
      {
        bool are_all_voxels_leafs = true;
        for(const auto& node_ijk : edge_voxels)
        {
          const std::size_t x = coords_uniform[0] + df * node_ijk[0];
          const std::size_t y = coords_uniform[1] + df * node_ijk[1];
          const std::size_t z = coords_uniform[2] + df * node_ijk[2];
          // check for overflow / ignore edges on boundary
          if(x >= dim_ || y >= dim_ || z >= dim_)
          {
            are_all_voxels_leafs = false;
            break;
          }

          const Node_index n = get_node(x, y, z);
          if(octree_.depth(n) > depth)
          {
            are_all_voxels_leafs = false;
            break;
          }
        }

        if(are_all_voxels_leafs)
        {
          // add to leaf edge set
          std::size_t e_gl = e_glIndex(edge_voxels[0][3],
                                       coords_global[0], coords_global[1], coords_global[2],
                                       depth);
          leaf_edges_set.insert({e_gl, depth});
        }
      }
    }

    leaf_voxels_ = std::vector<Voxel_handle>{leaf_voxels_set.begin(), leaf_voxels_set.end()};
    leaf_edges_ = std::vector<Edge_handle>{leaf_edges_set.begin(), leaf_edges_set.end()};
    leaf_vertices_ = std::vector<Vertex_handle>{leaf_vertices_set.begin(), leaf_vertices_set.end()};*/
  }

  const Octree& octree() const
  {
    return octree_;
  }

  const Geom_traits& geom_traits() const
  {
    return gt_;
  }

  std::size_t dim() const
  {
    return dim_;
  }

  FT hx() const
  {
    return hx_;
  }

  FT offset_x() const
  {
    return offset_x_;
  }

  FT offset_y() const
  {
    return offset_y_;
  }

  FT offset_z() const
  {
    return offset_z_;
  }

  std::size_t max_depth() const
  {
    return max_depth_;
  }

  const std::vector<Edge_handle>& leaf_edges() const
  {
    return leaf_edges_;
  }

  const std::vector<Vertex_handle>& leaf_vertices() const
  {
    return leaf_vertices_;
  }

  const std::vector<Voxel_handle>& leaf_voxels() const
  {
    return leaf_voxels_;
  }

  std::size_t depth(const std::size_t lex) const {
    Node_index n = get_node(lex);
    return octree_.depth(n);
  }

  std::size_t depth_factor(const std::size_t depth) const
  {
    return std::size_t(1) << (max_depth_ - depth);
  }

  Uniform_coords uniform_coordinates(Node_index node_index) const
  {
    Uniform_coords coords = octree_.global_coordinates(node_index);
    const std::size_t df = depth_factor(octree_.depth(node_index));
    for(int i=0; i<3; ++i)
      coords[i] *= static_cast<uint32_t>(df);

    return coords;
  }

  Point_3 point(const Uniform_coords& vertex_coordinates) const
  {
    const FT x0 = offset_x_ + vertex_coordinates[0] * hx_;
    const FT y0 = offset_y_ + vertex_coordinates[1] * hx_;
    const FT z0 = offset_z_ + vertex_coordinates[2] * hx_;

    return { x0, y0, z0 };
  }

  Point_3 point(const Vertex_handle& v) const
  {
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(v, max_depth_);

    const FT x0 = offset_x_ + i * hx_;
    const FT y0 = offset_y_ + j * hx_;
    const FT z0 = offset_z_ + k * hx_;

    return { x0, y0, z0 };
  }

  Uniform_coords vertex_uniform_coordinates(Node_index node_index,
                                            const typename Octree::Local_coordinates local_coords) const
  {
    const Uniform_coords& node_coords = octree_.global_coordinates(node_index);
    Uniform_coords v_coords = node_coords;
    for(int i=0; i<3; ++i)
      v_coords[i] += std::size_t(local_coords[i]);

    const std::size_t df = depth_factor(octree_.depth(node_index));
    for(int i=0; i<3; ++i)
      v_coords[i] *= static_cast<uint32_t>(df);

    return v_coords;
  }

  Node_index get_node(const std::size_t i,
                      const std::size_t j,
                      const std::size_t k) const
  {
    Node_index node_index = octree_.root();
    const std::size_t x = i;
    const std::size_t y = j;
    const std::size_t z = k;

    while(!octree_.is_leaf(node_index))
    {
      std::size_t dist_to_max = max_depth_ - octree_.depth(node_index) - 1;
      typename Octree::Local_coordinates loc;
      if(x & (std::size_t(1) << dist_to_max))
          loc[0] = true;

      if(y & (std::size_t(1) << dist_to_max))
          loc[1] = true;

      if(z & (std::size_t(1) << dist_to_max))
          loc[2] = true;

      node_index = octree_.child(node_index, loc.to_ulong());
    }

    return node_index;
  }

  Node_index get_node(const std::size_t lex_index) const
  {
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(lex_index, max_depth_);
    return get_node(i, j, k);
  }

  std::size_t lex_index(const std::size_t i,
                        const std::size_t j,
                        const std::size_t k,
                        const std::size_t depth) const
  {
    std::size_t dim = (std::size_t(1) << depth) + 1;
    return k * dim * dim + j * dim + i;
  }

  std::tuple<std::size_t, std::size_t, std::size_t> ijk_index(const std::size_t lex_index,
                                                              const std::size_t depth) const
  {
    const std::size_t dim = (std::size_t(1) << depth) + 1;
    return std::make_tuple(lex_index % dim, (lex_index / dim) % dim, lex_index / (dim * dim));
  }

  std::array<Vertex_handle, 2> edge_vertices(const Edge_handle& e_id) const
  {
    namespace Tables = internal::Cube_table;

    std::size_t e_global_id, depth;
    std::tie(e_global_id, depth) = static_cast<const std::tuple<std::size_t,std::size_t>&>(e_id);
    const std::size_t df = depth_factor(depth);

    const size_t v0_lex_index = e_global_id / 3;
    std::size_t i0, j0, k0;
    std::tie(i0, j0, k0) = ijk_index(v0_lex_index, depth);

    // v1
    const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];
    const int* v1_local = Tables::local_vertex_location[Tables::edge_to_vertex[e_local_index][1]];

    const std::size_t i1 = i0 + v1_local[0];
    const std::size_t j1 = j0 + v1_local[1];
    const std::size_t k1 = k0 + v1_local[2];

    const std::size_t v0 = lex_index(df * i0, df * j0, df * k0, max_depth_);
    const std::size_t v1 = lex_index(df * i1, df * j1, df * k1, max_depth_);

    return { v0, v1 };
  }

  /// returns the 4 voxels incident to an edge. If an edge has only three incident
  /// voxels, one will appear twice. The voxels are given with the uniform
  /// lexicographical index.
  std::array<std::size_t, 4> edge_voxels(const Edge_handle& e_id) const
  {
    namespace Tables = internal::Cube_table;

    std::size_t e_global_id, depth;
    std::tie(e_global_id, depth) = static_cast<const std::tuple<std::size_t,std::size_t>&>(e_id);
    const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];

    const std::size_t df = depth_factor(depth);

    const size_t v0_lex_index = e_global_id / 3;
    std::size_t i, j, k;
    std::tie(i, j, k) = ijk_index(v0_lex_index, depth);
    i *= df;
    j *= df;
    k *= df;

    const auto& voxel_neighbors = Tables::edge_to_voxel_neighbor[e_local_index];
    Node_index n0 = get_node(i + voxel_neighbors[0][0], j + voxel_neighbors[0][1], k + voxel_neighbors[0][2]);
    Node_index n1 = get_node(i + voxel_neighbors[1][0], j + voxel_neighbors[1][1], k + voxel_neighbors[1][2]);
    Node_index n2 = get_node(i + voxel_neighbors[2][0], j + voxel_neighbors[2][1], k + voxel_neighbors[2][2]);
    Node_index n3 = get_node(i + voxel_neighbors[3][0], j + voxel_neighbors[3][1], k + voxel_neighbors[3][2]);

    const Uniform_coords n0_uniform_coords = uniform_coordinates(n0);
    const Uniform_coords n1_uniform_coords = uniform_coordinates(n1);
    const Uniform_coords n2_uniform_coords = uniform_coordinates(n2);
    const Uniform_coords n3_uniform_coords = uniform_coordinates(n3);

    std::size_t n0_lex = lex_index(n0_uniform_coords[0], n0_uniform_coords[1], n0_uniform_coords[2], max_depth_);
    std::size_t n1_lex = lex_index(n1_uniform_coords[0], n1_uniform_coords[1], n1_uniform_coords[2], max_depth_);
    std::size_t n2_lex = lex_index(n2_uniform_coords[0], n2_uniform_coords[1], n2_uniform_coords[2], max_depth_);
    std::size_t n3_lex = lex_index(n3_uniform_coords[0], n3_uniform_coords[1], n3_uniform_coords[2], max_depth_);

    return { n0_lex, n1_lex, n2_lex, n3_lex };
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_OCTREE_WRAPPER_H
