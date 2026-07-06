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

#include <CGAL/AABB_tree/internal/AABB_two_trees_traversal.h>
#include <CGAL/AABB_tree/internal/AABB_two_trees_traversal_traits.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#endif

/// \file AABB_trees/intersection.h

namespace CGAL{
namespace AABB_trees {

  /// \ingroup PkgAABBTreeRef
  ///
  /// \brief Tests whether two AABB trees contain intersecting primitives.
  ///
  /// \cgalNamedParamsBegin
  ///   \cgalParamNBegin{transformation}
  ///     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm1` (`tm2`)}
  ///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
  ///                    as key type and `bool` as value type}
  ///     \cgalParamDefault{a constant property map returning `false` for any edge}
  ///   \cgalParamNEnd
  /// \cgalNamedParamsEnd
  ///
  /// Returns `true` if at least one primitive of `tree1` intersects
  /// a primitive of `tree2`, and `false` otherwise.
  template< typename Concurrency_tag = Sequential_tag,
            typename AABBTree1,
            typename AABBTree2,
            typename NamedParameters1 = parameters::Default_named_parameters,
            typename NamedParameters2 = parameters::Default_named_parameters>
  bool do_intersect(const AABBTree1 &tree1,
                    const AABBTree2 &tree2,
                    const NamedParameters1& np1 = parameters::default_values(),
                    const NamedParameters2& np2 = parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    if constexpr(is_default_parameter<NamedParameters1, internal_np::transformation_t>::value &&
                 is_default_parameter<NamedParameters2, internal_np::transformation_t>::value)
    {
      CGAL::internal::AABB_tree::Two_trees_do_intersect_traits traversal_traits(tree1.traits(), tree2.traits());
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
      return traversal_traits.is_intersection_found();
    }
    else
    {
      using Kernel = typename Kernel_traits<typename AABBTree1::AABBTraits::Point>::Kernel::type;
      using Aff_tr = Aff_transformation_3<Kernel>;
      const Aff_tr& tr1 = choose_parameter(get_parameter(np1, internal_np::transformation), Aff_tr(Identity_transformation()));
      const Aff_tr& tr2 = choose_parameter(get_parameter(np2, internal_np::transformation), Aff_tr(Identity_transformation()));
      CGAL::internal::AABB_tree::Two_trees_do_intersect_traits_with_transformation traversal_traits(tree1.traits(), tree2.traits(), tr1, tr2);
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
      return traversal_traits.is_intersection_found();
    }
  }

  /// \ingroup PkgAABBTreeRef
  ///
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
            typename OutputIterator,
            typename NamedParameters1 = parameters::Default_named_parameters,
            typename NamedParameters2 = parameters::Default_named_parameters>
  void all_pairs_of_intersecting_primitives(const AABBTree1 &tree1,
                                            const AABBTree2 &tree2,
                                            OutputIterator out,
                                            const NamedParameters1& np1 = parameters::default_values(),
                                            const NamedParameters2& np2 = parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;
    using parameters::is_default_parameter;
    if constexpr(is_default_parameter<NamedParameters1, internal_np::transformation_t>::value &&
                 is_default_parameter<NamedParameters2, internal_np::transformation_t>::value)
    {
      CGAL::internal::AABB_tree::Two_trees_listing_intersecting_primitives_traits traversal_traits(tree1.traits(), tree2.traits(), out);
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
    }
    else
    {
      using Kernel = typename Kernel_traits<typename AABBTree1::AABB_traits::Point>::Kernel;
      using Aff_tr = Aff_transformation_3<Kernel>;
      Aff_tr tr1 = choose_parameter(get_parameter(np1, internal_np::transformation), Aff_tr(Identity_transformation()));
      Aff_tr tr2 = choose_parameter(get_parameter(np2, internal_np::transformation), Aff_tr(Identity_transformation()));
      CGAL::internal::AABB_tree::Two_trees_listing_intersecting_primitives_traits_with_transformation traversal_traits(tree1.traits(), tree2.traits(), out, tr1, tr2);
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
    }
  }

  /// \ingroup PkgAABBTreeRef
  ///
  /// \brief Computes all self-intersecting primitive pairs in a single AABB tree.
  ///
  /// Intersections of a primitive with itself are not reported, and each intersecting
  /// pair of distinct primitives is reported only once.
  template< typename Concurrency_tag = Sequential_tag,
            typename AABBTree,
            typename OutputIterator>
  void all_pairs_of_intersecting_primitives(const AABBTree &tree,
                                            OutputIterator out)
  {
    CGAL::internal::AABB_tree::Listing_self_intersecting_primitives_traits traversal_traits(tree.traits(), out);
    CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree, tree, traversal_traits);
  }

}} // end namespace CGAL::AABB_trees

#endif
