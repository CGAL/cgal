// Copyright (c) 2026  Geometry Factory (France).
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

#ifndef CGAL_AABB_TWO_TREES_TRAVERSAL_H
#define CGAL_AABB_TWO_TREES_TRAVERSAL_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>

namespace CGAL {

namespace internal { namespace AABB_tree {

template <bool in_order = true,
          typename ConcurrencyTag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename TwoTreeTraversalTraits>
void two_trees_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                        const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                        const std::size_t nb_primitives_A,
                        const std::size_t nb_primitives_B,
                        TwoTreeTraversalTraits& traversal_traits)
{
#if CGAL_LINKED_WITH_TBB
  const std::size_t cutoff_parallel_call = 3000;
#endif
  auto recursive_call = [](const auto &node_A, const auto &node_B, const std::size_t nb_primitives_A, const std::size_t nb_primitives_B, auto &traversal_traits){
    if(traversal_traits.prefer_A_for_next_step(node_A, node_B, nb_primitives_A, nb_primitives_B))
      two_trees_traversal< in_order, ConcurrencyTag>(node_A, node_B, nb_primitives_A, nb_primitives_B, traversal_traits);
    else
      two_trees_traversal<!in_order, ConcurrencyTag>(node_B, node_A, nb_primitives_B, nb_primitives_A, traversal_traits);
  };
  switch(nb_primitives_A)
  {
  case 2:
  {
    if constexpr(in_order){
      traversal_traits.intersection(node_A.left_data(), node_B, nb_primitives_B);
      traversal_traits.intersection(node_A.right_data(), node_B, nb_primitives_B);
    } else {
      traversal_traits.intersection(node_B, nb_primitives_B, node_A.left_data());
      traversal_traits.intersection(node_B, nb_primitives_B, node_A.right_data());
    }
    break;
  }
  case 3:
  {
    if constexpr(in_order)
      traversal_traits.intersection(node_A.left_data(), node_B, nb_primitives_B);
    else
      traversal_traits.intersection(node_B, nb_primitives_B, node_A.left_data());

    bool do_intersect_right;
    if constexpr(in_order)
      do_intersect_right = traversal_traits.do_intersect(node_A.right_child(), node_B);
    else
      do_intersect_right = traversal_traits.do_intersect(node_B, node_A.right_child());
    if( do_intersect_right )
      two_trees_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B, 2, traversal_traits);
    break;
  }
  default:
  {
    bool do_intersect_left, do_intersect_right;
    if constexpr(in_order){
      do_intersect_left  = traversal_traits.do_intersect(node_A.left_child(), node_B);
      do_intersect_right = traversal_traits.do_intersect(node_A.right_child(), node_B);
    } else {
      do_intersect_left  = traversal_traits.do_intersect(node_B, node_A.left_child());
      do_intersect_right = traversal_traits.do_intersect(node_B, node_A.right_child());
    }
#if CGAL_LINKED_WITH_TBB
    if constexpr(ConcurrencyTag::is_parallel)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                recursive_call(node_A.left_child(), node_B, nb_primitives_A/2, nb_primitives_B, traversal_traits);}
              );
        recursive_call(node_A.right_child(), node_B, nb_primitives_A - nb_primitives_A/2, nb_primitives_B, traversal_traits);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          recursive_call(node_A.left_child(), node_B, nb_primitives_A/2, nb_primitives_B, traversal_traits);
        if( traversal_traits.go_further() && do_intersect_right )
          recursive_call(node_A.right_child(), node_B, nb_primitives_A - nb_primitives_A/2, nb_primitives_B, traversal_traits);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        recursive_call(node_A.left_child(), node_B, nb_primitives_A/2, nb_primitives_B, traversal_traits);
      if( traversal_traits.go_further() && do_intersect_right )
        recursive_call(node_A.right_child(), node_B, nb_primitives_A - nb_primitives_A/2, nb_primitives_B, traversal_traits);
    }
  }} // switch end
}

template<typename ConcurrencyTag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename TwoTreeTraversalTraits>
void two_trees_traversal(const Tree_A& tree_A,
                         const Tree_B& tree_B,
                         TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_trees_traversal<true, ConcurrencyTag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), traits);
}

// An optimization to traverse a tree with itself
template <typename ConcurrencyTag = Sequential_tag,
          typename AABBTraits,
          typename TwoTreeTraversalTraits>
void one_tree_traversal(const ::CGAL::AABB_node<AABBTraits>& node,
                        const std::size_t nb_primitives,
                        TwoTreeTraversalTraits& traversal_traits)
{
#if CGAL_LINKED_WITH_TBB
  const std::size_t cutoff_parallel_call = 3000;
#endif
  switch(nb_primitives)
  {
  case 2:
  {
    traversal_traits.intersection(node.left_data(), node.right_data());
    break;
  }
  case 3:
  {
    traversal_traits.intersection(node.left_data(), node.right_child(), 2);
    one_tree_traversal(node.right_child(), 2, traversal_traits);
    break;
  }
  default:
  {
#if CGAL_LINKED_WITH_TBB
    if constexpr(ConcurrencyTag::is_parallel)
    {
      if(nb_primitives > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        if(traversal_traits.do_intersect(node.left_child(), node.right_child())){
          tg.run([&]{
                  two_trees_traversal<true, ConcurrencyTag>(node.left_child(), node.right_child(), nb_primitives/2, nb_primitives - nb_primitives/2, traversal_traits);}
          );
        }
        tg.run([&]{
                one_tree_traversal<ConcurrencyTag>(node.left_child(), nb_primitives/2, traversal_traits);}
              );
        one_tree_traversal<ConcurrencyTag>(node.right_child(), nb_primitives - nb_primitives/2, traversal_traits);
        tg.wait();
      }
      else
      {
        if(traversal_traits.do_intersect(node.left_child(), node.right_child()))
          two_trees_traversal(node.left_child(), node.right_child(), nb_primitives/2, nb_primitives - nb_primitives/2, traversal_traits);
        one_tree_traversal(node.left_child(), nb_primitives/2, traversal_traits);
        one_tree_traversal(node.right_child(), nb_primitives - nb_primitives/2, traversal_traits);
      }
    }
    else
#endif
    {
      if(traversal_traits.do_intersect(node.left_child(), node.right_child()))
        two_trees_traversal(node.left_child(), node.right_child(), nb_primitives/2, nb_primitives - nb_primitives/2, traversal_traits);
      one_tree_traversal(node.left_child(), nb_primitives/2, traversal_traits);
      one_tree_traversal(node.right_child(), nb_primitives - nb_primitives/2, traversal_traits);
    }
  }} // switch end
}

template<typename ConcurrencyTag = Sequential_tag,
         typename Tree,
         typename TwoTreeTraversalTraits>
void one_tree_traversal(const Tree& tree,
                        TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree.size() != 0);
  one_tree_traversal<ConcurrencyTag>(*tree.root_node(), tree.size(), traits);
}

namespace experimental{

template <bool in_order = true,
          typename ConcurrencyTag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename TwoTreeTraversalTraits>
void two_trees_partial_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                                const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                                const std::size_t nb_primitives_A,
                                const std::size_t nb_primitives_B,
                                const std::size_t cutoff,
                                TwoTreeTraversalTraits& traversal_traits)
{
#if CGAL_LINKED_WITH_TBB
  const std::size_t cutoff_parallel_call = 3000;
#endif
  auto recursive_call = [&](const auto &node_A, const auto &node_B, const std::size_t nb_primitives_A, const std::size_t nb_primitives_B, auto &traversal_traits){
    if(traversal_traits.prefer_A_for_next_step(node_A, node_B, nb_primitives_A, nb_primitives_B))
      two_trees_partial_traversal< in_order, ConcurrencyTag>(node_A, node_B, nb_primitives_A, nb_primitives_B, cutoff, traversal_traits);
    else
      two_trees_partial_traversal<!in_order, ConcurrencyTag>(node_B, node_A, nb_primitives_B, nb_primitives_A, cutoff, traversal_traits);
  };
  if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    if constexpr(in_order)
      traversal_traits.intersection(node_A, node_B);
    else
      traversal_traits.intersection(node_B, node_A);
  }
  else if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    two_trees_partial_traversal<!in_order>(node_B, node_A, nb_primitives_B, nb_primitives_A, cutoff, traversal_traits);
  }
  else
  {
    bool do_intersect_left, do_intersect_right;
    if constexpr(in_order){
      do_intersect_left  = traversal_traits.do_intersect(node_A.left_child(), node_B);
      do_intersect_right = traversal_traits.do_intersect(node_A.right_child(), node_B);
    } else {
      do_intersect_left  = traversal_traits.do_intersect(node_B, node_A.left_child());
      do_intersect_right = traversal_traits.do_intersect(node_B, node_A.right_child());
    }
#if CGAL_LINKED_WITH_TBB
    if constexpr(ConcurrencyTag::is_parallel)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                recursive_call(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);
              });
        recursive_call(node_B, node_A.right_child(), nb_primitives_B, nb_primitives_A-nb_primitives_A/2, traversal_traits);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          recursive_call(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);
        if( traversal_traits.go_further() && do_intersect_right )
          recursive_call(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, traversal_traits);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        recursive_call(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);
      if( traversal_traits.go_further() && do_intersect_right )
        recursive_call(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, traversal_traits);
    }
  }
}

template<typename ConcurrencyTag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename TwoTreeTraversalTraits>
void two_trees_partial_traversal(const Tree_A& tree_A,
                                const Tree_B& tree_B,
                                const std::size_t cutoff,
                                TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_trees_partial_traversal<true, ConcurrencyTag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), cutoff, traits);
}

} // end of namespace experimental

}}} // end of namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
