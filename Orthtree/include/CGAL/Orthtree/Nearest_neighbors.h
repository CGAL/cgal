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

#ifndef ORTHTREE_NEAREST_NEIGHBORS_H
#define ORTHTREE_NEAREST_NEIGHBORS_H

#include <CGAL/license/Orthtree.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>
#include <boost/function_output_iterator.hpp>

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

    // Loop through each of the elements contained by the node
    // Note: there might be none, and that should be fine!
    for (auto& e: orthtree.data(node)) {

      // Pair that element with its distance from the search point
      Result current_element_with_distance =
        {e, orthtree.traits().squared_distance_of_element_object()(e, orthtree.traits().construct_center_3_object()(search_bounds))};

      // Check if the new element is within the bounds
      if (current_element_with_distance.distance < orthtree.traits().compute_squared_radius_3_object()(search_bounds)) {

        // Check if the results list is full
        if (results.size() == results.capacity()) {

          // Delete a element if we need to make room
          results.pop_back();
        }

        // Add the new element
        results.push_back(current_element_with_distance);

        // Sort the list
        std::sort(results.begin(), results.end(), [=](auto& left, auto& right) {
          return left.distance < right.distance;
        });

        // Check if the results list is full
        if (results.size() == results.capacity()) {

          // Set the search radius
          search_bounds = orthtree.traits().construct_sphere_3_object()(orthtree.traits().construct_center_3_object()(search_bounds), results.back().distance + epsilon);
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
    children_with_distances.reserve(Tree::degree);

    // Fill the list with child nodes
    for (int i = 0; i < Tree::degree; ++i) {
      auto child_node = orthtree.child(node, i);

      // Add a child to the list, with its distance
      children_with_distances.emplace_back(
        child_node,
        CGAL::squared_distance(orthtree.traits().construct_center_3_object()(search_bounds), orthtree.barycenter(child_node))
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
  \ingroup PkgOrthtreeNeighbors
  \brief finds at most `k` elements within a specific radius that are
  nearest to the center of the sphere `query`: if `query` does not contain
  at least `k` elements, only contained elements will be returned.

  This function is useful when the user already knows how sparse the elements are,
  or if they do not care about elements that are too far away.
  Setting a small radius may have performance benefits.

  \tparam Tree must be `Orthtree<GT>` with `GT` being a model of `CollectionPartitioningOrthtreeTraits`
  \tparam OutputIterator must be a model of `OutputIterator` that accepts `GT::Node_data_element`

  \param orthtree the tree to search within
  \param query the region to search within
  \param k the number of elements to find
  \param output the output iterator to add the found elements to (in order of increasing distance)
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_k_neighbors_in_radius(
  const Tree& orthtree,
  typename Tree::Sphere& query,
  std::size_t k,
  OutputIterator output
) {

  // todo: this type is over-constrained, this must be made more generic
  struct Node_element_with_distance {
    typename Tree::Traits::Node_data_element element;
    typename Tree::FT distance;
  };

  // Create an empty list of elements
  std::vector<Node_element_with_distance> element_list;
  if (k != (std::numeric_limits<std::size_t>::max)())
    element_list.reserve(k);

  // Invoking the recursive function adds those elements to the vector (passed by reference)
  CGAL::internal::nearest_k_neighbors_recursive(orthtree, query, orthtree.root(), element_list);

  // Add all the points found to the output
  for (auto& item: element_list)
    *output++ = item.element;

  return output;
}


/*!
  \ingroup PkgOrthtreeNeighbors
  \brief finds the `k` nearest neighbors of the point `query`.

  Nearest neighbors are outputted in order of increasing distance to
  `query`.

  \tparam Tree must be `Orthtree<GT>` with `GT` being a model of `CollectionPartitioningOrthtreeTraits`
  \tparam OutputIterator a model of `OutputIterator` that accepts `GT::Node_data_element` objects.

  \param orthtree the tree to search within
  \param query query point
  \param k number of neighbors to find
  \param output output iterator
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_neighbors(const Tree& orthtree, const typename Tree::Point& query,
                                 std::size_t k,
                                 OutputIterator output) {
  typename Tree::Sphere query_sphere = orthtree.traits().construct_sphere_3_object()(query, (std::numeric_limits<typename Tree::FT>::max)());
  return nearest_k_neighbors_in_radius(orthtree, query_sphere, k, output);
}

/*!
  \ingroup PkgOrthtreeNeighbors
  \brief finds the elements in the sphere `query`.

  Elements are outputted in order of increasing distance to
  the center of the sphere.

  \tparam Tree must be `Orthtree<GT>` with `GT` being a model of `CollectionPartitioningOrthtreeTraits`
  \tparam OutputIterator a model of `OutputIterator` that accepts `GT::Node_data_element` objects.

  \param orthtree the tree to search within
  \param query query sphere
  \param output output iterator
 */
template <typename Tree, typename OutputIterator>
OutputIterator nearest_neighbors(const Tree& orthtree, const typename Tree::Sphere& query, OutputIterator output) {
  typename Tree::Sphere query_sphere = query;
  return nearest_k_neighbors_in_radius(orthtree, query_sphere, (std::numeric_limits<std::size_t>::max)(), output);
}

}
}

#endif //ORTHTREE_NEAREST_NEIGHBORS_H
