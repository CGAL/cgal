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

#include <memory>

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal.h>
#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal_traits.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#include <CGAL/mutex.h>
#endif

namespace CGAL {

  /// returns `true`, iff at least one primitive of `tree1` does intersect a primitive of `tree2`.
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

/// puts in `out` the ids of all pairs of a primitive of `tree1` and a primitive of `tree2` that do intersect.
template< typename Concurrency_tag = Sequential_tag,
          typename AABBTree1,
          typename AABBTree2,
          typename OutputIterator >
void all_intersected_primitives(const AABBTree1 &tree1,
                                const AABBTree2 &tree2,
                                OutputIterator out)
{
  CGAL::internal::AABB_tree::Two_tree_listing_intersecting_primitives_traits traversal_traits(tree1.traits(), tree2.traits(), out);
  CGAL::internal::AABB_tree::two_tree_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
}

/// puts in `out` the ids of all pairs of intersected primitives in the given tree.
template< typename Concurrency_tag = Sequential_tag,
          typename AABBTree,
          typename OutputIterator>
void all_intersected_primitives(const AABBTree &tree,
                                OutputIterator out)
{
  all_intersected_primitives<Concurrency_tag>(tree, tree, out);
}


} // end namespace CGAL::AABB_tree

#endif
