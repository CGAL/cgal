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

#ifndef CGAL_AABB_TREES_INTERSECTIONS_H
#define CGAL_AABB_TREES_INTERSECTIONS_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/AABB_tree/internal/AABB_two_trees_traversal.h>
#include <CGAL/AABB_tree/internal/AABB_two_trees_traversal_traits.h>

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/FPU.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#endif

/// \file AABB_trees/intersection.h

namespace CGAL{
namespace AABB_trees {
  /// \ingroup PkgAABBTreeRef
  ///
  /// \brief tests if at least two primitives each from an AABB tree intersect.
  ///
  /// \cgalNamedParamsBegin
  ///   \cgalParamNBegin{concurrency_tag}
  ///     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
  ///     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
  ///     \cgalParamDefault{`CGAL::Sequential_tag`}
  ///     \cgalParamExtra{`np1` only}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{transformation}
  ///     \cgalParamDescription{An affine transformation apply to `tree1` (`tree2`)}
  ///     \cgalParamType{`CGAL::Aff_transformation_3<Kernel>` where `Kernel` is the kernel associated with `AABBTree1::AABB_traits::Point` (`AABBTree2::AABB_traits::Point`)}
  ///     \cgalParamDefault{An identity transformation}
  ///   \cgalParamNEnd
  /// \cgalNamedParamsEnd
  ///
  /// \return `true` if at least one primitive of `tree1` intersects
  /// a primitive of `tree2`, and `false` otherwise.
  template< typename AABBTree1,
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
    using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters1,
                                          Sequential_tag
                                        > ::type;

    if constexpr(is_default_parameter<NamedParameters1, internal_np::transformation_t>::value &&
                 is_default_parameter<NamedParameters2, internal_np::transformation_t>::value)
    {
      CGAL::internal::AABB_tree::Two_trees_do_intersect_traits traversal_traits(tree1.traits(), tree2.traits());
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
      return traversal_traits.is_intersection_found();
    }
    else
    {
      using Kernel = typename Kernel_traits<typename AABBTree1::AABB_traits::Point>::Kernel;
      // Get the dimension of the AABBTraits through the bbox and get the appropriate Aff_transformation
      using Aff_tr = std::conditional_t< std::is_same_v<typename AABBTree1::Bounding_box, Bbox_2>,
                                         Aff_transformation_2<Kernel>,
                                         Aff_transformation_3<Kernel>>;
      const Aff_tr& tr1 = choose_parameter(get_parameter(np1, internal_np::transformation), Aff_tr(Identity_transformation()));
      const Aff_tr& tr2 = choose_parameter(get_parameter(np2, internal_np::transformation), Aff_tr(Identity_transformation()));
      CGAL::internal::AABB_tree::Two_trees_do_intersect_traits_with_transformation<typename AABBTree1::AABB_traits,
                                                                                   typename AABBTree2::AABB_traits,
                                                                                   Aff_tr, false>
                                                        traversal_traits(tree1.traits(), tree2.traits(), tr1, tr2);
      CGAL::internal::AABB_tree::two_trees_traversal(tree1, tree2, traversal_traits);
      return traversal_traits.is_intersection_found();
    }
  }

  /// \ingroup PkgAABBTreeRef
  ///
  /// \brief Computes all pairs of intersecting primitive from two AABB trees.
  ///
  /// Traverses both trees and outputs all pairs of primitives that intersect.
  /// Each output element is a pair `(id1, id2)` where:
  ///       - `id1` is the ID of a primitive from `tree1`
  ///       - `id2` is the ID of a primitive from `tree2`
  ///
  /// \tparam AABBTree1        Type of the first AABB tree.
  /// \tparam AABBTree2        Type of the second AABB tree.
  /// \tparam OutputIterator   Output iterator storing std::pair<AABBTree1::Primitive::Id, AABBTree2::Primitive::Id>.
  /// \tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
  /// \tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
  ///
  /// \cgalNamedParamsBegin
  ///   \cgalParamNBegin{concurrency_tag}
  ///     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
  ///     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
  ///     \cgalParamDefault{`CGAL::Sequential_tag`}
  ///     \cgalParamExtra{`np1` only}
  ///   \cgalParamNEnd
  ///   \cgalParamNBegin{transformation}
  ///     \cgalParamDescription{An affine transformation apply to `tree1` (`tree2`)}
  ///     \cgalParamType{`CGAL::Aff_transformation_3<Kernel>` where `Kernel` is the kernel associated with `AABBTree1::AABB_traits::Point` (`AABBTree2::AABB_traits::Point`)}
  ///     \cgalParamDefault{An identity transformation}
  ///   \cgalParamNEnd
  /// \cgalNamedParamsEnd
  ///
  /// \see do_intersect()
  template< typename AABBTree1,
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
    using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                          internal_np::concurrency_tag_t,
                                          NamedParameters1,
                                          Sequential_tag
                                        > ::type;

    if constexpr(is_default_parameter<NamedParameters1, internal_np::transformation_t>::value &&
                 is_default_parameter<NamedParameters2, internal_np::transformation_t>::value)
    {
      CGAL::internal::AABB_tree::Two_trees_listing_intersecting_primitives_traits traversal_traits(tree1.traits(), tree2.traits(), out);
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
    }
    else
    {
      using Kernel = typename Kernel_traits<typename AABBTree1::AABB_traits::Point>::Kernel;
      // Get the dimension of the AABBTraits through the bbox and get the appropriate Aff_transformation
      using Aff_tr = std::conditional_t< std::is_same_v<typename AABBTree1::Bounding_box, Bbox_2>,
                                         Aff_transformation_2<Kernel>,
                                         Aff_transformation_3<Kernel>>;
      Aff_tr tr1 = choose_parameter(get_parameter(np1, internal_np::transformation), Aff_tr(Identity_transformation()));
      Aff_tr tr2 = choose_parameter(get_parameter(np2, internal_np::transformation), Aff_tr(Identity_transformation()));
      CGAL::internal::AABB_tree::Two_trees_listing_intersecting_primitives_traits_with_transformation<typename AABBTree1::AABB_traits,
                                                                                                      typename AABBTree2::AABB_traits,
                                                                                                      OutputIterator, Aff_tr, false>
                                                        traversal_traits(tree1.traits(), tree2.traits(), out, tr1, tr2);
      CGAL::internal::AABB_tree::two_trees_traversal<Concurrency_tag>(tree1, tree2, traversal_traits);
    }
  }

  /// \ingroup PkgAABBTreeRef
  ///
  /// \brief Computes all pairs of primitives from a single AABB tree that are intersecting.
  ///
  /// \note Whether two objects that share a common subfeature (e.g., two triangles sharing an edge) are considered to intersect depends on the primitive type used.
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
