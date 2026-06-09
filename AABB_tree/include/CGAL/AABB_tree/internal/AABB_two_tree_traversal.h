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

#ifndef CGAL_AABB_TWO_TREE_TRAVERSAL_H
#define CGAL_AABB_TWO_TREE_TRAVERSAL_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <optional>

namespace CGAL {

namespace internal { namespace AABB_tree {

template <bool in_order = true,
          typename Concurrency_tag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename TwoTreeTraversalTraits>
void two_tree_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                        const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                        const std::size_t nb_primitives_A,
                        const std::size_t nb_primitives_B,
                        TwoTreeTraversalTraits& traversal_traits)
{
  const std::size_t cutoff_parallel_call = 100000;
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

    if( traversal_traits.do_intersect(node_A.right_child(), node_B) )
      two_tree_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B, 2, traversal_traits);
    break;
  }
  default:
  {
    // TODO current strategy is we swap at each step. Try largest tree first or largest bbox first
    bool do_intersect_left  = traversal_traits.do_intersect(node_A.left_child(), node_B);
    bool do_intersect_right = traversal_traits.do_intersect(node_A.right_child(), node_B);
#if CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                two_tree_traversal<!in_order, Concurrency_tag>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);}
              );
        two_tree_traversal<!in_order, Concurrency_tag>(node_B, node_A.right_child(), nb_primitives_B, nb_primitives_A-nb_primitives_A/2, traversal_traits);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          two_tree_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);
        if( traversal_traits.go_further() && do_intersect_right )
          two_tree_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B, nb_primitives_A-nb_primitives_A/2, traversal_traits);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        two_tree_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traversal_traits);
      if( traversal_traits.go_further() && do_intersect_right )
        two_tree_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, traversal_traits);
    }
  }} // switch end
}

template<typename Concurrency_tag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename TwoTreeTraversalTraits>
void two_tree_traversal(const Tree_A& tree_A,
                        const Tree_B& tree_B,
                        TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_tree_traversal<true, Concurrency_tag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), traits);
}

namespace experimental{

template <bool in_order = true,
          typename Concurrency_tag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename TwoTreeTraversalTraits>
void two_tree_partial_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                                const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                                const std::size_t nb_primitives_A,
                                const std::size_t nb_primitives_B,
                                const std::size_t cutoff,
                                TwoTreeTraversalTraits& traversal_traits)
{
  const std::size_t cutoff_parallel_call = 100000;
  if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    if constexpr(in_order)
      traversal_traits.intersection(node_A, node_B);
    else
      traversal_traits.intersection(node_B, node_A);
  }
  else if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    two_tree_partial_traversal<!in_order>(node_B, node_A, nb_primitives_B, nb_primitives_A, cutoff, traversal_traits);
  }
  else
  {
    // TODO current strategy is we swap at each step. Try largest tree first or largest bbox first
    bool do_intersect_left  = traversal_traits.do_intersect(node_A.left_child(), node_B);
    bool do_intersect_right = traversal_traits.do_intersect(node_A.right_child(), node_B);
#if CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                two_tree_partial_traversal<!in_order, Concurrency_tag>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traversal_traits);}
              );
        two_tree_partial_traversal<!in_order, Concurrency_tag>(node_B, node_A.right_child(), nb_primitives_B, nb_primitives_A-nb_primitives_A/2, cutoff, traversal_traits);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          two_tree_partial_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traversal_traits);
        if( traversal_traits.go_further() && do_intersect_right )
          two_tree_partial_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B, nb_primitives_A-nb_primitives_A/2, cutoff, traversal_traits);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        two_tree_partial_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traversal_traits);
      if( traversal_traits.go_further() && do_intersect_right )
        two_tree_partial_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, cutoff, traversal_traits);
    }
  }
}

template<typename Concurrency_tag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename TwoTreeTraversalTraits>
void two_tree_partial_traversal(const Tree_A& tree_A,
                                const Tree_B& tree_B,
                                const std::size_t cutoff,
                                TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_tree_partial_traversal<true, Concurrency_tag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), cutoff, traits);
}

} // end of namespace experimental

}}} // end of namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
