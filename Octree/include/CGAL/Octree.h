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

#include <CGAL/Octree/Octree_node.h>
#include <CGAL/Octree/Split_criterion.h>
#include <CGAL/Octree/Tree_walker_criterion.h>

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
#include <ostream>

#include <stack>
#include <queue>
#include <vector>
#include <math.h>

namespace CGAL {

  template<class PointRange,
          class PointMap>
  class Octree {
  public: // types

    // Deduce the kernel
    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;

    // Define the Node based on this kernel
    typedef Octree_node<Kernel, PointRange> Node;

    // Some other useful types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;

    // New Types :
    typedef typename PointRange::iterator Range_iterator;
    typedef typename std::iterator_traits<Range_iterator>::value_type Range_type;
    // TODO: Kernel can be deduced from the point map

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
      FT max_len = (x_len < y_len) ? y_len : x_len;
      max_len = (max_len < z_len) ? z_len : max_len;
      bbox = Iso_cuboid(bbox.min(), bbox.min() + max_len * Vector(1.0, 1.0, 1.0));

      // Shift the squared box to make sure it's centered in the original place
      Point bbox_transformed_centroid = midpoint(bbox.min(), bbox.max());
      Vector diff_centroid = bbox_centroid - bbox_transformed_centroid;
      bbox = bbox.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));

      // save octree attributes
      m_bbox_min = bbox.min();
      m_bbox_side = bbox.max()[0] - m_bbox_min[0];
      m_root.begin() = pwn.begin();
      m_root.end() = pwn.end();
    }

    ~Octree() {
      m_root.unsplit();
    }

    void refine(std::function<bool(const Node &)> split_criterion) {

      // create a side length map
      for (int i = 0; i <= (int) 10; i++)
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

    void refine(size_t max_depth, size_t bucket_size) {
      refine(Split_to_max_depth_or_bucket_size(max_depth, bucket_size));
    }

    Node &root() { return m_root; }

    const Node &root() const { return m_root; }

    bool operator==(Octree<PointRange, PointMap> &rhs) {

      // Identical trees should have the same bounding box
      if (rhs.m_bbox_min != m_bbox_min || rhs.m_bbox_side != m_bbox_side)
        return false;

      // Identical trees should have the same depth
      if (rhs.m_max_depth_reached != m_max_depth_reached)
        return false;

      // If all else is equal, recursively compare the trees themselves
      return rhs.m_root == m_root;
    }


  private: // functions :

    Point compute_barycenter_position(Node &node) const {

      // Determine the side length of this node
      FT size = m_side_per_depth[node.depth()];

      // Determine the location this node should be split
      FT bary[3];
      for (int i = 0; i < 3; i++)
        bary[i] = node.location()[i] * size + (size / 2.0) + m_bbox_min[i];

      // Convert that location into a point
      return {bary[0], bary[1], bary[2]};
    }

    std::array<Range_iterator, 9> partition_around_point(Range_iterator begin, Range_iterator end, Point point) {

      auto partitions = std::array<Range_iterator, 9>();

      partitions[0] = begin;
      partitions[8] = end;

      // Split on x
      partitions[4] = std::partition(partitions[0], partitions[8],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[0] < point[0]);
                                     });

      // Split on y, to the left of x
      partitions[2] = std::partition(partitions[0], partitions[4],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[1] < point[1]);
                                     });

      // Split on y, to the right of x
      partitions[6] = std::partition(partitions[4], partitions[8],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[1] < point[1]);
                                     });

      // Split on z, to the left of y and the left of x
      partitions[1] = std::partition(partitions[0], partitions[2],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[2] < point[2]);
                                     });

      // Split on z, to the right of y and the left of x
      partitions[3] = std::partition(partitions[2], partitions[4],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[2] < point[2]);
                                     });

      // Split on z, to the left of y and the right of x
      partitions[5] = std::partition(partitions[4], partitions[6],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[2] < point[2]);
                                     });

      // Split on z, to the right of y and the right of x
      partitions[7] = std::partition(partitions[6], partitions[8],
                                     [&](const Range_type &a) -> bool {
                                       return (get(m_points_map, a)[2] < point[2]);
                                     });

      return partitions;
    }

    void reassign_points(Node &node) {

      Point center = compute_barycenter_position(node);

      auto partitions = partition_around_point(node.begin(), node.end(), center);

      for (int i = 0; i < 8; ++i) {

        node[i].begin() = partitions[i];
        node[i].end() = partitions[i + 1];
      }
    }

    void _split_node(Node &node,
                     Range_iterator begin, Range_iterator end,
                     std::size_t dimension, std::bitset<3> coords,
                     Point center) {

      // Split the point collection around the center point on this dimension
      Range_iterator split_point = std::partition(begin, end,
                                                  [&](const Range_type &a) -> bool {
                                                    return (get(m_points_map, a)[dimension] < split_point[dimension]);
                                                  });
    }

    void _reassign_points(Node &node) {

    }
  }; // end class Octree

} // namespace CGAL

#endif // CGAL_OCTREE_3_H
