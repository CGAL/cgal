// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Tong Zhao, Cédric Portaneri

#ifndef CGAL_OCTREE_3_H
#define CGAL_OCTREE_3_H

/*
 * Not present or relevant for benchmarking
 */
//#include <CGAL/license/Implicit_surface_reconstruction_3.h>

#include <CGAL/Octree/Octree_node.h>

#include <CGAL/bounding_box.h>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

/*
 * These headers were not included here originally
 * Adding them was necessary to make this header self sufficient
 */
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include <fstream>

#include <stack>
#include <queue>
#include <vector>
#include <math.h>

namespace CGAL {

  struct HashIntPoint_3 {
    uint64_t operator()(const IntPoint_3 &pt) const {
      return ((std::hash<uint64_t>()(static_cast<uint64_t>(pt.x())) ^
               (std::hash<uint64_t>()(static_cast<uint64_t>(pt.y())) << 1)) >> 1)
             ^ (std::hash<uint64_t>()(static_cast<uint64_t>(pt.z())) << 1);
    }
  };

  enum DebugOctreeVisuType {
    SHOW_ALL_LEAFS = 0,
    SHOW_NON_EMPTY_LEAFS = 1,
    SHOW_NON_EMPTY_NODES = 2
  };

  template<class Kernel,
          class PointRange,
          class PointMap,
          class NormalMap>
  class Octree {
  public: // types :
    typedef Octree_node<Kernel, PointRange> Node;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef IntPoint_3 IntPoint;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;
    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

  private: // data members :
    Node m_root;                      /* root node of the octree */
    uint8_t m_max_depth_reached = 0;  /* octree actual highest depth reached */
    PointRange &m_ranges;              /* input point range */
    PointMap m_points_map;          /* property map: `value_type of InputIterator` -> `Point` (Position) */
    Point m_bbox_min;                  /* input bounding box min value */
    FT m_bbox_side;              /* input bounding box side length (cube) */
    std::vector<FT> m_side_per_depth;      /* side length per node's depth */
    std::vector<size_t> m_unit_per_depth; /* number of unit node (smallest) inside one node for each depth for one axis */

  public: // functions :
    Octree(
            PointRange &pwn,
            PointMap &point_map,
            const FT enlarge_ratio = 1.2) :
            m_ranges(pwn),
            m_points_map(point_map) {
      // compute bbox
      Iso_cuboid bbox = CGAL::bounding_box(boost::make_transform_iterator
                                                   (m_ranges.begin(),
                                                    CGAL::Property_map_to_unary_function<PointMap>(
                                                            m_points_map)),
                                           boost::make_transform_iterator
                                                   (m_ranges.end(),
                                                    CGAL::Property_map_to_unary_function<PointMap>(
                                                            m_points_map)));

      Point bbox_centroid = midpoint(bbox.min(), bbox.max());

      // scale bbox
      bbox = bbox.transform(Aff_transformation_3<Kernel>(SCALING, enlarge_ratio));

      // use isotrope bbox
      FT x_len = bbox.xmax() - bbox.xmin();
      FT y_len = bbox.ymax() - bbox.ymin();
      FT z_len = bbox.zmax() - bbox.zmin();
      FT max_len = (x_len < y_len) ? y_len : x_len;
      max_len = (max_len < z_len) ? z_len : max_len;
      bbox = Iso_cuboid(bbox.min(), bbox.min() + max_len * Vector(1.0, 1.0, 1.0));

      // translate bbox back to initial centroid
      Point bbox_transformed_centroid = midpoint(bbox.min(), bbox.max());
      Vector diff_centroid = bbox_centroid - bbox_transformed_centroid;
      bbox = bbox.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));

      // save octree attributes
      m_bbox_min = bbox.min();
      m_bbox_side = bbox.max()[0] - m_bbox_min[0];
      for (InputIterator it = pwn.begin(); it != pwn.end(); it++)
        m_root.add_point(it);
    }

    ~Octree() {
      m_root.unsplit();
    }

    // template < typename CellCriteria, typename NormalCriteria > // or other useful criterion
    void refine(size_t max_depth, size_t max_pts_num) {
      if (max_depth < 0 || max_pts_num < 1) {
        CGAL_TRACE_STREAM << "wrong octree refinement criteria\n";
        return;
      }
      for (int i = 0; i <= (int) max_depth; i++)
        m_side_per_depth.push_back(m_bbox_side / (FT) (1 << i));

      refine_recurse(&m_root, max_depth, max_pts_num);

      for (int i = 0; i <= (int) m_max_depth_reached; i++)
        m_unit_per_depth.push_back(1 << (m_max_depth_reached - i));
    }

    void grade() {
      std::queue<Node *> leaf_nodes;
      fill_leaf_queue(&m_root, leaf_nodes);

      while (!leaf_nodes.empty()) {
        Node *node = leaf_nodes.front();
        leaf_nodes.pop();
        if (!node->is_leaf()) continue;

        std::list<Node *> neighbors_to_split = node->find_unbalanced_neighbors_to_split();
        if (!neighbors_to_split.empty()) leaf_nodes.push(node);
        for (Node *neighbor : neighbors_to_split) {
          neighbor->split();
          reassign_points(neighbor);
          for (int child_id = 0; child_id < 8; child_id++) {
            Node *neighbor_child = neighbor->child(child_id);
            leaf_nodes.push(neighbor_child);
          }
        }
      }
    }

    template<class OutputPwnIterator, class OutputPointIterator>
    void generate_points(OutputPwnIterator output_pwn, OutputPointIterator output_steiner) {
      std::set<IntPoint> all_corner_locations;
      std::map<IntPoint, Vector> corner_normals;

      // 1. for each leaf, insert barycenter pwn and get corner unique location
      std::queue<Node *> leaf_nodes;
      fill_leaf_queue(&m_root, leaf_nodes);
      while (!leaf_nodes.empty()) {
        Node *node = leaf_nodes.front();
        leaf_nodes.pop();

        IntPoint node_corners_location[8];
        for (int child_id = 0; child_id < 8; child_id++) {
          node_corners_location[child_id] = get_corner_location(node, child_id);
          all_corner_locations.insert(node_corners_location[child_id]);
        }

        Point bary_position = compute_barycenter_position(node);
        Vector bary_normal = CGAL::NULL_VECTOR;
        for (const InputIterator &pwn_it : node->points()) {
          bary_normal += compute_weighted_normal(bary_position, pwn_it);

          for (int child_id = 0; child_id < 8; child_id++) {
            IntPoint corner_loc = node_corners_location[child_id];
            Vector weighted_normal = compute_weighted_normal(compute_corner_position(corner_loc),
                                                             pwn_it);

            auto find = corner_normals.find(corner_loc);
            if (find == corner_normals.end()) {
              corner_normals.insert(std::make_pair(corner_loc, weighted_normal));
            } else {
              find->second += weighted_normal;
            }
          }
        }

        if (bary_normal == CGAL::NULL_VECTOR) {
          *output_steiner = bary_position;
          output_steiner++;
        } else {
          bary_normal /= CGAL::sqrt(
                  bary_normal.squared_length()); // will be optional with weighted normal computation fonctor
          *output_pwn = std::make_pair(bary_position, bary_normal);
          output_pwn++;
        }
      }

      // 2. for each corner location, insert corner pwn
      for (const IntPoint &corner_location : all_corner_locations) {
        Point corner_position = compute_corner_position(corner_location);

        auto find = corner_normals.find(corner_location);
        if (find == corner_normals.end() || find->second == CGAL::NULL_VECTOR) {
          *output_steiner = corner_position;
          output_steiner++;
        } else {
          Vector corner_normal = find->second;
          corner_normal /= CGAL::sqrt(
                  corner_normal.squared_length()); // will be optional with weighted normal computation fonctor
          *output_pwn = std::make_pair(corner_position, corner_normal);
          output_pwn++;
        }
      }
    }

    void fill_leaf_queue(Node *node, std::queue<Node *> &queue) {
      if (node->is_leaf()) {
        queue.push(node);
      } else {
        for (int child_id = 0; child_id < 8; child_id++) {
          fill_leaf_queue(node->child(child_id), queue);
        }
      }
    }

    Node *root() { return &m_root; }

    const Node *root() const { return &m_root; }

    size_t num_corner() {
      std::set<IntPoint> all_corner_locations;
      std::queue<Node *> leaf_nodes;
      fill_leaf_queue(&m_root, leaf_nodes);
      while (!leaf_nodes.empty()) {
        Node *node = leaf_nodes.front();
        leaf_nodes.pop();

        IntPoint node_corners_location[8];
        for (int child_id = 0; child_id < 8; child_id++) {
          node_corners_location[child_id] = get_corner_location(node, child_id);
          all_corner_locations.insert(node_corners_location[child_id]);
        }
      }

      return all_corner_locations.size();
    }

  private: // functions :

    Point compute_barycenter_position(Node *node) const {
      FT size = m_side_per_depth[node->depth()];
      FT bary[3];
      for (int i = 0; i < 3; i++)
        bary[i] = node->location()[i] * size + (size / 2.0) + m_bbox_min[i];
      return Point(bary[0], bary[1], bary[2]);
    }

    IntPoint get_corner_location(Node *node, size_t corner_idx) const {
      size_t size = m_unit_per_depth[node->depth()];
      IntPoint loc;
      for (int i = 0; i < 3; i++)
        loc[i] = (node->location()[i] + ((corner_idx >> i) & 1)) * size;
      return loc;
    }

    Point compute_corner_position(const IntPoint &corner_loc) const {
      FT size = m_side_per_depth[m_max_depth_reached];
      FT corner_pos[3];
      for (int i = 0; i < 3; i++)
        corner_pos[i] = (FT) (corner_loc[i] * size) + m_bbox_min[i];
      return Point(corner_pos[0], corner_pos[1], corner_pos[2]);
    }

    Vector compute_weighted_normal(const Point &corner_loc, const InputIterator &pwn_it) const {
      return pwn_it->second; // simple case for now, will be replaced by a generic functor
    }

    void refine_recurse(Node *node, size_t dist_to_max_depth, size_t max_pts_num) {
      if (dist_to_max_depth == 0 || node->num_points() <= max_pts_num) {
        if (m_max_depth_reached < node->depth()) m_max_depth_reached = node->depth();
        return;
      }

      node->split();
      reassign_points(node);
      for (int child_id = 0; child_id < 8; child_id++) {
        refine_recurse(node->child(child_id), dist_to_max_depth - 1, max_pts_num);
      }
    }

    void reassign_points(Node *node) {
      Point barycenter = compute_barycenter_position(node);

      for (const InputIterator &pwn_it : node->points()) {
        const Point &point = get(m_points_map, *pwn_it);

        int is_right = barycenter[0] < point[0];
        int is_up = barycenter[1] < point[1];
        int is_front = barycenter[2] < point[2];

        bool equal_right = std::abs(barycenter[0] - point[0]) < 1e-6;
        bool equal_up = std::abs(barycenter[1] - point[1]) < 1e-6;
        bool equal_front = std::abs(barycenter[2] - point[2]) < 1e-6;

        int child_id = (is_front << 2) | (is_up << 1) | is_right;
        node->child(child_id)->add_point(pwn_it);

        if (equal_right) {
          int sym_child_id = (is_front << 2) | (is_up << 1) | (!is_right);
          node->child(sym_child_id)->add_point(pwn_it);
        }

        if (equal_up) {
          int sym_child_id = (is_front << 2) | (!is_up << 1) | is_right;
          node->child(sym_child_id)->add_point(pwn_it);
        }

        if (equal_front) {
          int sym_child_id = (!is_front << 2) | (is_up << 1) | (!is_right);
          node->child(sym_child_id)->add_point(pwn_it);
        }

      }
    }

  public: // debugging :

    void dump_header(int num_cuboids, std::ofstream &out_file) const {
      out_file << "ply\n" // header
                  "format ascii 1.0\n"
                  "element vertex "
               << num_cuboids * 8 << "\n"
                                     "property float x\n"
                                     "property float y\n"
                                     "property float z\n"
                                     "element edge "
               << num_cuboids * 12 << "\n"
                                      "property int vertex1\n"
                                      "property int vertex2\n"
                                      "end_header\n";
    }

    void dump_cuboid_vertices(const Point &min, const Point &max, std::ofstream &out_file) const {
      for (int child_id = 0; child_id < 8; child_id++) {
        for (int j = 0; j < 3; j++) {
          if (((child_id >> j) & 1) == 0)
            out_file << min[j] << " ";
          else
            out_file << max[j] << " ";
        }
        out_file << "\n";
      }
    }

    void dump_cuboid_vertices(const Point &barycenter, const FT &half_size, std::ofstream &out_file) const {
      Vector half_size_vec = half_size * Vector(1., 1., 1.);
      dump_cuboid_vertices(barycenter - half_size_vec, barycenter + half_size_vec, out_file);
    }

    void dump_cuboid_edges(int cuboid_id, std::ofstream &out_file) const {
      for (int i = 0; i < 8; i++) {
        int v1 = (cuboid_id * 8) + i;
        if ((i & 3) == 3) v1 -= 3;
        int v2 = v1 + 1 + (i & 1);
        out_file << v1 << " " << v2 << "\n";
      }
      for (int i = 0; i < 4; i++) {
        int v1 = (cuboid_id * 8) + i;
        int v2 = v1 + 4;
        out_file << v1 << " " << v2 << "\n";
      }
    }

    void dump_bbox(const std::string &filename) const {
      std::cout << " dump bbox " + filename + "\n";
      std::ofstream out_file(filename + ".ply");
      dump_header(1, out_file);
      dump_cuboid_vertices(m_bbox_min, m_bbox_min + m_bbox_side * Vector(1, 1, 1), out_file);
      dump_cuboid_edges(0, out_file);
      out_file.close();
    }

    void dump_octree(const std::string &filename, DebugOctreeVisuType visu_type) {
      std::cout << " dump octree " + filename + "\n";
      std::ofstream out_file(filename + ".ply");
      int num_cuboid_to_draw = 0;
      get_num_cuboid_to_draw(&m_root, num_cuboid_to_draw, visu_type);

      dump_header(num_cuboid_to_draw, out_file);
      dump_octree_vertices_recursive(&m_root, out_file, visu_type);
      for (int i = 0; i < num_cuboid_to_draw; i++) {
        dump_cuboid_edges(i, out_file);
      }
      out_file.close();
    }

    void get_num_cuboid_to_draw(Node *node, int &num_cuboid_to_draw, DebugOctreeVisuType visu_type) {
      bool is_leaf = node->is_leaf();
      bool is_non_empty = !node->is_empty();

      if ((visu_type == SHOW_ALL_LEAFS && is_leaf) ||
          (visu_type == SHOW_NON_EMPTY_NODES && is_non_empty) ||
          (visu_type == SHOW_NON_EMPTY_LEAFS && is_non_empty && is_leaf)) {
        num_cuboid_to_draw++;
      }
      if (!is_leaf) {
        for (int child_id = 0; child_id < 8; child_id++) {
          get_num_cuboid_to_draw(node->child(child_id), num_cuboid_to_draw, visu_type);
        }
      }
    }

    void dump_octree_vertices_recursive(Node *node, std::ofstream &out_file, DebugOctreeVisuType visu_type) {
      bool is_leaf = node->is_leaf();
      bool is_non_empty = !node->is_empty();

      if ((visu_type == SHOW_ALL_LEAFS && is_leaf) ||
          (visu_type == SHOW_NON_EMPTY_NODES && is_non_empty) ||
          (visu_type == SHOW_NON_EMPTY_LEAFS && is_non_empty && is_leaf)) {
        dump_cuboid_vertices(compute_barycenter_position(node), m_side_per_depth[node->depth()] / 2.0,
                             out_file);
      }
      if (!is_leaf) {
        for (int child_id = 0; child_id < 8; child_id++) {
          dump_octree_vertices_recursive(node->child(child_id), out_file, visu_type);
        }
      }
    }

    bool debug_grading() {
      return debug_grading_recursive(&m_root);
    }

    bool debug_grading_recursive(Node *node) {
      if (node->is_leaf()) {
        return node->is_balanced();
      }
      for (int child_id = 0; child_id < 8; child_id++) {
        if (!debug_grading_recursive(node->child(child_id)))
          return false;
      }
      return true;
    }

    void dump_octree_pwn(const std::string &filename, const PointRange &ranges) const {
      std::cout << " dump octree points with normal " + filename + "\n";
      std::ofstream out_file(filename + ".ply");
      out_file << "ply\n" // header
                  "format ascii 1.0\n"
                  "element vertex "
               << ranges.size() << "\n"
                                   "property float x\n"
                                   "property float y\n"
                                   "property float z\n"
                                   "property float nx\n"
                                   "property float ny\n"
                                   "property float nz\n"
                                   "end_header\n";

      for (InputIterator pwn_it = ranges.cbegin(); pwn_it != ranges.cend(); pwn_it++) {
        const Point &point = get(m_points_map, *pwn_it);

        out_file << point[0] << " " << point[1] << " " << point[2] << "\n";
      }
    }

    void dump_octree_point(const std::string &filename, const std::vector<Point> &ranges) const {
      std::cout << " dump octree points " + filename + "\n";
      std::ofstream out_file(filename + ".ply");
      out_file << "ply\n" // header
                  "format ascii 1.0\n"
                  "element vertex "
               << ranges.size() << "\n"
                                   "property float x\n"
                                   "property float y\n"
                                   "property float z\n"
                                   "end_header\n";

      for (const Point &point : ranges) {
        out_file << point[0] << " " << point[1] << " " << point[2] << "\n";
      }
    }

  }; // end class Octree

} // namespace CGAL

#endif // CGAL_OCTREE_3_H
