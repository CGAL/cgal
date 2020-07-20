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

#include <CGAL/Octree/Node.h>
#include <CGAL/Octree/Split_criterion.h>
#include <CGAL/Octree/Walker_criterion.h>
#include <CGAL/Octree/Walker_iterator.h>

#include <CGAL/bounding_box.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>

#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <fstream>
#include <ostream>
#include <functional>

#include <stack>
#include <queue>
#include <vector>
#include <math.h>

namespace CGAL {

namespace Octree {

/*!
 * \ingroup PkgOctreeClasses
 *
 * \brief Class Octree is a data structure for efficient computations in 3D space.
 *
 * \details It builds a heirarchy of nodes which subdivide the space based on a collection of points.
 * Each node represents an axis aligned cubic region of space.
 * A node contains the range of points that are present in the region it defines,
 * and it may contain eight other nodes which further subdivide the region.
 *
 * \tparam PointRange is a range type that provides random access iterators over the indices of a set of points.
 * \tparam PointMap is a type that maps items in the PointRange to Point data
 */
template<class PointRange, class PointMap>
class Octree {

public:

  /// \name Public Types
  /// @{

  /*!
   * \brief The point type is deduced from the type of the property map used
   */
  typedef typename boost::property_traits<PointMap>::value_type Point;

  /*!
   * \brief The Kernel used is deduced from the point type
   */
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;

  /*!
   * \brief The floating point type is decided by the Kernel
   */
  typedef typename Kernel::FT FT;

  /*!
   * \brief
   */
  typedef boost::iterator_range<typename PointRange::iterator> Points_iterator_range;

  /*!
   * \brief The Sub-tree / Octant type
   */
  typedef Node::Node<typename PointRange::iterator> Node;

  /*!
   * \brief A function that determines whether a node needs to be split when refining a tree
   */
  typedef std::function<bool(const Node &)> Split_criterion_function;

  /*!
   * \brief A range that provides input-iterator access to the nodes of a tree
   */
  typedef boost::iterator_range<Walker_iterator<const Node>> Node_range;

  /*!
   * \brief A function that determines the next node in a traversal given the current one
   */
  typedef std::function<const Node *(const Node *)> Node_walker;

  /// @}

private: // Private types

  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;
  typedef typename Kernel::Sphere_3 Sphere;
  typedef typename CGAL::Bbox_3 Bbox;
  typedef typename PointRange::iterator Range_iterator;
  typedef typename std::iterator_traits<Range_iterator>::value_type Range_type;

private: // data members :

  Node m_root;                      /* root node of the octree */
  uint8_t m_max_depth_reached = 0;  /* octree actual highest depth reached */

  PointRange &m_ranges;              /* input point range */
  PointMap m_points_map;          /* property map: `value_type of InputIterator` -> `Point` (Position) */

  Point m_bbox_min;                  /* input bounding box min value */
  FT m_bbox_side;              /* input bounding box side length (cube) */

  std::vector<FT> m_side_per_depth;      /* side length per node's depth */
  std::vector<size_t> m_unit_per_depth; /* number of unit node (smallest) inside one node for each depth for one axis */

public:

  /// \name Construction, Destruction
  /// @{

  /*!
   * \brief Create an octree from a collection of points
   *
   * The resulting octree will have a root node with no children that contains the points passed.
   * That root node will have a bounding box that encloses all of the points passed,
   * with padding according to the enlarge_ratio
   *
   * \param point_range random access iterator over the indices of the points
   * \param point_map maps the point indices to their coordinate locations
   * \param enlarge_ratio the degree to which the bounding box should be enlarged
   */
  Octree(
          PointRange &point_range,
          PointMap &point_map,
          const FT enlarge_ratio = 1.2) :
          m_ranges(point_range),
          m_points_map(point_map) {


    // compute bounding box that encloses all points
    Iso_cuboid bbox = CGAL::bounding_box(boost::make_transform_iterator
                                                 (m_ranges.begin(),
                                                  CGAL::Property_map_to_unary_function<PointMap>(
                                                          m_points_map)),
                                         boost::make_transform_iterator
                                                 (m_ranges.end(),
                                                  CGAL::Property_map_to_unary_function<PointMap>(
                                                          m_points_map)));

    // Find the center point of the box
    Point bbox_centroid = midpoint(bbox.min(), bbox.max());

    // scale bounding box to add padding
    bbox = bbox.transform(Aff_transformation_3<Kernel>(SCALING, enlarge_ratio));

    // Convert the bounding box into a cube
    FT x_len = bbox.xmax() - bbox.xmin();
    FT y_len = bbox.ymax() - bbox.ymin();
    FT z_len = bbox.zmax() - bbox.zmin();
    FT max_len = std::max({x_len, y_len, z_len});
    bbox = Iso_cuboid(bbox.min(), bbox.min() + max_len * Vector(1.0, 1.0, 1.0));

    // Shift the squared box to make sure it's centered in the original place
    Point bbox_transformed_centroid = midpoint(bbox.min(), bbox.max());
    Vector diff_centroid = bbox_centroid - bbox_transformed_centroid;
    bbox = bbox.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));

    // save octree attributes
    m_bbox_min = bbox.min();
    m_bbox_side = bbox.max()[0] - m_bbox_min[0];
    m_root.points() = {point_range.begin(), point_range.end()};
  }

  /// @}

  /// \name Tree Building
  /// @{

  /*!
   * \brief Subdivide an octree's nodes and sub-nodes until it meets the given criteria
   *
   * \todo
   *
   * \param split_criterion rule to use when determining whether or not a node needs to be subdivided
   */
  void refine(const Split_criterion_function &split_criterion) {

    // create a side length map
    for (int i = 0; i <= (int) 32; i++)
      m_side_per_depth.push_back(m_bbox_side / (FT) (1 << i));

    // Initialize a queue of nodes that need to be refined
    std::queue<Node *> todo;
    todo.push(&m_root);

    // Process items in the queue until it's consumed fully
    while (!todo.empty()) {

      // Get the next element
      auto current = todo.front();
      todo.pop();
      int depth = current->depth();

      // Check if this node needs to be processed
      if (split_criterion(*current)) {

        // Split this node
        current->split();

        // Redistribute its points
        reassign_points((*current));

        // Process each of its children
        for (int i = 0; i < 8; ++i)
          todo.push(&(*current)[i]);

      }
    }
  }

  /*!
   * \brief Refine an octree using a max depth and max number of points in a node as split criterion
   *
   * \todo
   *
   * \param max_depth deepest a tree is allowed to be (nodes at this depth will not be split)
   * \param bucket_size maximum points a node is allowed to contain
   */
  void refine(size_t max_depth, size_t bucket_size) {
    refine(Split_criterion::Max_depth_or_bucket_size(max_depth, bucket_size));
  }

  /// @}

  /// \name Accessors
  /// @{

  /*!
   * \brief Provides read and write access to the root node, and by extension the rest of the tree
   *
   * \todo
   *
   * \return a reference to the root node of the tree
   */
  Node &root() { return m_root; }

  /*!
   * \brief Provides read-only access to the root node, and by extension the rest of the tree
   *
   * \todo
   *
   * \return a const reference to the root node of the tree
   */
  const Node &root() const { return m_root; }

  /*!
   * \brief Constructs an input range of nodes using a tree walker function
   *
   * \todo
   *
   * \tparam Walker type of the walker rule
   * \param walker the rule to use when determining the order of the sequence of points produced
   * \return a forward input iterator over the nodes of the tree
   */
  template<class Walker>
  Node_range walk(const Walker &walker = Walker()) const {

    const Node *first = walker.first(&m_root);

    Node_walker next = std::bind(&Walker::template next<Points_iterator_range>,
                                 walker, std::placeholders::_1);

    return boost::make_iterator_range(Walker_iterator<const Node>(first, next),
                                      Walker_iterator<const Node>());
  }

  /*!
   * \brief Find the leaf node which would contain a point
   *
   * Traverses the octree and finds the deepest cell that has a domain enclosing the point passed.
   *
   * \param p The point to find a node for
   * \return A reference to the node which would contain the point
   */
  const Node &locate(const Point &p) const {

    // Start at the root node
    auto *node_for_point = &m_root;

    // Descend the tree until reaching a leaf node
    while (!node_for_point->is_leaf()) {

      // Find the point to split around
      Point center = compute_barycenter_position(*node_for_point);

      // Find the index of the correct sub-node
      typename Node::Index index;
      for (int dimension = 0; dimension < 3; ++dimension) {

        index[dimension] = center[dimension] < p[dimension];
      }

      // Find the correct sub-node of the current node
      node_for_point = &(*node_for_point)[index.to_ulong()];
    }

    // Return the result
    return *node_for_point;
  }

  /*!
   * \brief Find the bounding box of a node
   *
   * \todo
   *
   * \param node the node to determine the bounding box of
   * \return the bounding box defined by that node's relationship to the tree
   */
  Bbox bbox(const Node &node) const {

    // Determine the side length of this node
    FT size = m_side_per_depth[node.depth()];

    // Determine the location this node should be split
    FT min_corner[3];
    FT max_corner[3];
    for (int i = 0; i < 3; i++) {

      min_corner[i] = m_bbox_min[i] + (node.location()[i] * size);
      max_corner[i] = min_corner[i] + size;
    }

    // Create the cube
    return {min_corner[0], min_corner[1], min_corner[2],
            max_corner[0], max_corner[1], max_corner[2]};
  }

  /*!
   * \brief Find the K points in a tree that are nearest to the search point
   *
   * \todo
   *
   * \tparam Point_output_iterator an output iterator type that will accept points
   * \param search_point the location to find points near
   * \param k the number of points to find
   * \param output the output iterator to add the found points to
   */
  template<typename Point_output_iterator>
  void nearest_k_neighbours(const Point &search_point, std::size_t k, Point_output_iterator output) const {

    // Create an empty list of points
    std::vector<Point_with_distance> points_list;
    points_list.reserve(k);

    // Invoking the recursive function adds those points to the vector (passed by reference)
//    nearest_k_neighbours_recursive(search_point, points_list, m_root, std::numeric_limits<FT>::max(), k);
    auto search_bounds = Sphere(search_point, std::numeric_limits<FT>::max());
    _nearest_k_neighbours_recursive(search_bounds, m_root, points_list);

    // Add all the points found to the output
    for (auto &item : points_list)
      *output++ = item.point;
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
   * \brief Compares the topology of a pair of Octrees
   *
   * \todo
   *
   * \param rhs tree to compare with
   * \return whether the trees have the same topology
   */
  bool operator==(const Octree<PointRange, PointMap> &rhs) const {

    // Identical trees should have the same bounding box
    if (rhs.m_bbox_min != m_bbox_min || rhs.m_bbox_side != m_bbox_side)
      return false;

    // Identical trees should have the same depth
    if (rhs.m_max_depth_reached != m_max_depth_reached)
      return false;

    // If all else is equal, recursively compare the trees themselves
    return rhs.m_root == m_root;
  }

  bool operator!=(const Octree<PointRange, PointMap> &rhs) const {
    return !operator==(rhs);
  }

  /// @}


private: // functions :

  Point compute_barycenter_position(const Node &node) const {

    // Determine the side length of this node
    FT size = m_side_per_depth[node.depth()];

    // Determine the location this node should be split
    FT bary[3];
    for (int i = 0; i < 3; i++)
      bary[i] = node.location()[i] * size + (size / 2.0) + m_bbox_min[i];

    // Convert that location into a point
    return {bary[0], bary[1], bary[2]};
  }

  void reassign_points(Node &node, Range_iterator begin, Range_iterator end, const Point &center,
                       std::bitset<3> coord = {},
                       std::size_t dimension = 0) {

    // Root case: reached the last dimension
    if (dimension == 3) {

      node[coord.to_ulong()].points() = {begin, end};

      return;
    }

    // Split the point collection around the center point on this dimension
    Range_iterator split_point = std::partition(begin, end,
                                                [&](const Range_type &a) -> bool {
                                                  return (get(m_points_map, a)[dimension] < center[dimension]);
                                                });

    // Further subdivide the first side of the split
    std::bitset<3> coord_left = coord;
    coord_left[dimension] = false;
    reassign_points(node, begin, split_point, center, coord_left, dimension + 1);

    // Further subdivide the second side of the split
    std::bitset<3> coord_right = coord;
    coord_right[dimension] = true;
    reassign_points(node, split_point, end, center, coord_right, dimension + 1);

  }

  void reassign_points(Node &node) {

    Point center = compute_barycenter_position(node);
    reassign_points(node, node.points().begin(), node.points().end(), center);
  }

  bool do_intersect(const Node &node, const Sphere &sphere) const {

    // Create a cubic bounding box from the node
    Bbox node_cube = bbox(node);

    // Check for overlap between the node's box and the sphere as a box, to quickly catch some cases
    // FIXME: Activating this causes slower times!
//    if (!do_overlap(node_cube, sphere.bbox()))
//      return false;

    // Check for intersection between the node and the sphere
    return CGAL::do_intersect(node_cube, sphere);
  }

  FT nearest_k_neighbours_recursive(const Point &p, std::vector<std::pair<Point, FT>> &out, const Node &node,
                                    FT search_bounds_radius_squared, std::size_t k) const {

    FT largest_radius_squared_found = search_bounds_radius_squared;

    // Check whether we've reached the bottom of the tree
    if (node.is_leaf()) {

      // Base case: the node has no children

      // Check if the node contains any points
      if (!node.is_empty()) {

        // If it does, loop through each point
        for (auto i : node.points()) {
          auto point = get(m_points_map, i);

          // Find the distance of the point
          FT new_distance_squared = CGAL::squared_distance(point, p);

          // Check whether the new distance is an improvement
          if (new_distance_squared < largest_radius_squared_found) {

            // If it is, add it to the list
            out.push_back({point, new_distance_squared});

            // Check whether the list is full
            if (out.size() >= k) {

              // If the list has at least K points, sort them
              std::sort(out.begin(), out.end(), [=](auto &left, auto &right) {
                return left.second < right.second;
              });

              // Make sure the list has only k points (discarding the furthest ones, if there are more
              out.resize(k);

              // Update the largest distance
              largest_radius_squared_found = out.back().second;
            }

          }
        }
      }

    } else {

      // If the node has children

      // Create a list of the child nodes
      std::vector<std::pair<typename Node::Index, FT>> children_with_distances;
      children_with_distances.reserve(8);

      // Check each node
      for (int index = 0; index < 8; ++index) {
        auto &n = node[index];

        // Find the distance of the node's center
        auto node_center_distance = CGAL::squared_distance(p, compute_barycenter_position(n));

        // Add this node to the list
        children_with_distances.emplace_back(index, node_center_distance);
      }

      // Sort the children by their distance
      std::sort(children_with_distances.begin(), children_with_distances.end(), [=](auto &left, auto &right) {
        return left.second < right.second;
      });

      // Search each of them
      for (auto child : children_with_distances) {
        auto &n = node[child.first.to_ulong()];

        // Check whether this node is capable of containing closer points
        if (do_intersect(n, Sphere{p, largest_radius_squared_found + 0.01 /*TODO: This is my epsilon*/})) {

          // Recursive case
          largest_radius_squared_found =
                  nearest_k_neighbours_recursive(p, out, n, largest_radius_squared_found, k);

        }
      }
    }

    return largest_radius_squared_found;
  }


  // TODO: There has to be a better way than using structs like these!
  struct Point_with_distance {
    Point point;
    FT distance;
  };
  struct Node_index_with_distance {
    typename Node::Index index;
    FT distance;
  };

  void _nearest_k_neighbours_recursive(Sphere &search_bounds, const Node &node,
                                       std::vector<Point_with_distance> &results) const {

    // Check whether the node has children
    if (node.is_leaf()) {

      // Base case: the node has no children

      // Loop through each of the points contained by the node
      // Note: there might be none, and that should be fine!
      for (auto point_index : node.points()) {

        // Retrieve each point from the octree's point map
        auto point = get(m_points_map, point_index);

        // Pair that point with its distance from the search point
        Point_with_distance current_point_with_distance =
                {point, squared_distance(point, search_bounds.center())};

        // Check if the new point is within the bounds
        if (current_point_with_distance.distance < search_bounds.squared_radius()) {

          // Check if the results list is full
          if (results.size() == results.capacity()) {

            // Delete a point if we need to make room
            results.pop_back();
          }

          // Add the new point
          results.push_back(current_point_with_distance);

          // Sort the list
          std::sort(results.begin(), results.end(), [=](auto &left, auto &right) {
            return left.distance < right.distance;
          });

          // Check if the results list is full
          if (results.size() == results.capacity()) {

            // Set the search radius
            search_bounds = Sphere(search_bounds.center(), results.back().distance);
          }
        }
      }
    } else {

      // Recursive case: the node has children

      // Create a list to map children to their distances
      std::vector<Node_index_with_distance> children_with_distances;
      children_with_distances.reserve(8);

      // Fill the list with child nodes
      for (int index = 0; index < 8; ++index) {
        auto &child_node = node[index];

        // Add a child to the list, with its distance
        children_with_distances.push_back(
                {typename Node::Index(index),
                 CGAL::squared_distance(search_bounds.center(), compute_barycenter_position(child_node))}
        );
      }

      // Sort the children by their distance from the search point
      std::sort(children_with_distances.begin(), children_with_distances.end(), [=](auto &left, auto &right) {
        return left.distance < right.distance;
      });

      // Loop over the children
      for (auto child_with_distance : children_with_distances) {
        auto &child_node = node[child_with_distance.index.to_ulong()];

        // Check whether the bounding box of the child intersects with the search bounds
        if (do_intersect(child_node, search_bounds)) {

          // Recursively invoke this function
          _nearest_k_neighbours_recursive(search_bounds, child_node, results);
        }
      }
    }
  }

}; // end class Octree

} // namespace Octree

} // namespace CGAL

#endif // CGAL_OCTREE_3_H
