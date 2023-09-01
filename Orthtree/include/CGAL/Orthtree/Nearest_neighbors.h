// Copyright (c) 2023 INRIA
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro

#ifndef ORTHTREE_EXAMPLES_NEAREST_NEIGHBORS_H
#define ORTHTREE_EXAMPLES_NEAREST_NEIGHBORS_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree.h>

namespace CGAL {

namespace internal {

template <typename Tree, typename Result>
void nearest_k_neighbors_recursive(const Tree& orthtree,
                                   typename Tree::Sphere& search_bounds,
                                   typename Tree::Node_index node,
                                   std::vector<Result>& results,
                                   typename Tree::FT epsilon = 0) {

  // Check whether the node has children
  if (orthtree.is_leaf(node)) {

    // Base case: the node has no children

    // Loop through each of the points contained by the node
    // Note: there might be none, and that should be fine!
    for (auto& p: orthtree.data(node)) {

      // Pair that point with its distance from the search point
      Result current_point_with_distance =
        {p, squared_distance(orthtree.traits().get_geometric_object_for_element_object()(p), search_bounds.center())};

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
        std::sort(results.begin(), results.end(), [=](auto& left, auto& right) {
          return left.distance < right.distance;
        });

        // Check if the results list is full
        if (results.size() == results.capacity()) {

          // Set the search radius
          search_bounds = typename Tree::Sphere(search_bounds.center(), results.back().distance + epsilon);
        }
      }
    }
  } else {

    struct Node_index_with_distance {
      typename Tree::Node_index index;
      typename Tree::FT distance;

      Node_index_with_distance(const typename Tree::Node_index& index, const typename Tree::FT& distance) :
        index(index), distance(distance) {}
    };

    // Recursive case: the node has children

    // Create a list to map children to their distances
    std::vector<Node_index_with_distance> children_with_distances;
    children_with_distances.reserve(Tree::Degree::value);

    // Fill the list with child nodes
    for (int i = 0; i < Tree::Degree::value; ++i) {
      auto child_node = orthtree.child(node, i);

      // Add a child to the list, with its distance
      children_with_distances.emplace_back(
        child_node,
        CGAL::squared_distance(search_bounds.center(), orthtree.barycenter(child_node))
      );
    }

    // Sort the children by their distance from the search point
    std::sort(children_with_distances.begin(), children_with_distances.end(), [=](auto& left, auto& right) {
      return left.distance < right.distance;
    });

    // Loop over the children
    for (auto child_with_distance: children_with_distances) {

      // Check whether the bounding box of the child intersects with the search bounds
      if (CGAL::do_intersect(orthtree.bbox(child_with_distance.index), search_bounds)) {

        // Recursively invoke this function
        CGAL::internal::nearest_k_neighbors_recursive(orthtree, search_bounds, child_with_distance.index, results);
      }
    }
  }
}

}

namespace Orthtrees {

/*!
  \brief finds the `k` points within a specific radius that are
  nearest to the center of `query_sphere`.

  This function guarantees that there are no closer points than the ones returned,
  but it does not guarantee that it will return at least `k` points.
  For a query where the search radius encloses `k` or fewer points, all enclosed points will be returned.
  If the search radius is too small, no points may be returned.
  This function is useful when the user already knows how sparse the points are,
  or if they do not care about points that are too far away.
  Setting a small radius may have performance benefits.

  \tparam Tree must be an orthtree with traits which are a model of CollectionPartitioningOrthtreeTraits
  \tparam OutputIterator must be a model of `OutputIterator` that accepts points
  \param orthtree the tree to search within
  \param query_sphere the region to search within
  \param k the number of points to find
  \param output the output iterator to add the found points to (in order of increasing distance)
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_k_neighbors_in_radius(
  const Tree& orthtree,
  typename Tree::Sphere& query_sphere,
  std::size_t k,
  OutputIterator output
) {

  // todo: this type is over-constrained, this must be made more generic
  struct Node_element_with_distance {
    typename Tree::Traits::Node_data_element point;
    typename Tree::FT distance;
  };

  // Create an empty list of points
  std::vector<Node_element_with_distance> points_list;
  if (k != (std::numeric_limits<std::size_t>::max)())
    points_list.reserve(k);

  // Invoking the recursive function adds those points to the vector (passed by reference)
  CGAL::internal::nearest_k_neighbors_recursive(orthtree, query_sphere, orthtree.root(), points_list);

  // Add all the points found to the output
  for (auto& item: points_list)
    *output++ = item.point;

  return output;
}


/*!
  \brief finds the `k` nearest neighbors of `query`.

  Nearest neighbors are outputted in order of increasing distance to
  `query`.

  \tparam Tree must be an orthtree with traits which are a model of CollectionPartitioningOrthtreeTraits
  \tparam OutputIterator a model of `OutputIterator` that accept `Point_d` objects.
  \param orthtree the tree to search within
  \param query query point.
  \param k number of neighbors.
  \param output output iterator.
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_neighbors(const Tree& orthtree, const typename Tree::Point& query,
                                 std::size_t k,
                                 OutputIterator output) {
  typename Tree::Sphere query_sphere(query, (std::numeric_limits<typename Tree::FT>::max)());
  return nearest_k_neighbors_in_radius(orthtree, query_sphere, k, output);
}

/*!
  \brief finds the points in `sphere`.

  Nearest neighbors are outputted in order of increasing distance to
  the center of `sphere`.

  \tparam Tree must be an orthtree with traits which are a model of CollectionPartitioningOrthtreeTraits
  \tparam OutputIterator a model of `OutputIterator` that accept `Point_d` objects.
  \param orthtree the tree to search within
  \param query query sphere.
  \param output output iterator.
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_neighbors(const Tree& orthtree, const typename Tree::Sphere& query, OutputIterator output) {
  typename Tree::Sphere query_sphere = query;
  return nearest_k_neighbors_in_radius(orthtree, query_sphere, (std::numeric_limits<std::size_t>::max)(), output);
}

}
}

#endif //ORTHTREE_EXAMPLES_NEAREST_NEIGHBORS_H
