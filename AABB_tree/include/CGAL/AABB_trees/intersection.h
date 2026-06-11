// Copyright (c) 2026  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_AABB_TREE_INTERSECTIONS_H
#define CGAL_AABB_TREE_INTERSECTIONS_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal.h>
#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal_traits.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#endif

namespace CGAL::AABB_trees {

  /// \brief Tests whether two AABB trees contain intersecting primitives.
  ///
  /// Returns `true` if at least one primitive of `tree1` intersects
  /// a primitive of `tree2`, and `false` otherwise.
  template< typename Concurrency_tag = Sequential_tag,
            typename AABBTree1,
            typename AABBTree2 >
  bool do_intersect(const AABBTree1 &tree1,
                    const AABBTree2 &tree2)
  {
    CGAL::internal::AABB_tree::Two_tree_do_intersect_traits traversal_traits(tree1.traits(), tree2.traits());
    CGAL::internal::AABB_tree::two_tree_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
    return traversal_traits.is_intersection_found();
  }

  /// \brief Computes all intersecting primitive pairs between two AABB trees.
  ///
  /// Traverses both trees and outputs all pairs of primitives that intersect.
  /// Each output element is a pair `(id1, id2)` where:
  ///       - `id1` is the ID of a primitive from `tree1`
  ///       - `id2` is the ID of a primitive from `tree2`
  ///
  /// \tparam Concurrency_tag  enables sequential versus parallel algorithm. Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
  /// \tparam AABBTree1        Type of the first AABB tree.
  /// \tparam AABBTree2        Type of the second AABB tree.
  /// \tparam OutputIterator   Output iterator storing pairs of primitive IDs.
  ///
  /// \see do_intersect()
  template< typename Concurrency_tag = Sequential_tag,
            typename AABBTree1,
            typename AABBTree2,
            typename OutputIterator >
  void all_pairs_of_intersecting_primitives(const AABBTree1 &tree1,
                                            const AABBTree2 &tree2,
                                            OutputIterator out)
  {
    CGAL::internal::AABB_tree::Two_tree_listing_intersecting_primitives_traits traversal_traits(tree1.traits(), tree2.traits(), out);
    CGAL::internal::AABB_tree::two_tree_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
  }

  /// \brief Computes all self-intersecting primitive pairs in a single AABB tree.
  ///
  /// \note This function is equivalent to calling:
  ///       `all_intersected_primitives(tree, tree, out)`.
  template< typename Concurrency_tag = Sequential_tag,
            typename AABBTree,
            typename OutputIterator>
  void all_pairs_of_intersecting_primitives(const AABBTree &tree,
                                            OutputIterator out)
  {
    all_pairs_of_intersecting_primitives<Concurrency_tag>(tree, tree, out);
  }

} // end namespace CGAL::AABB_trees

#endif
