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

#if CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#endif

namespace CGAL {

namespace internal { namespace AABB_tree {

template <typename AABBTraits_A,
          typename AABBTraits_B,
          typename TwoTreeTraversalTraits>
void sequential_two_trees_traversal(
  const CGAL::AABB_node<AABBTraits_A>& node_A,
  const CGAL::AABB_node<AABBTraits_B>& node_B,
  TwoTreeTraversalTraits& traversal_traits,
  bool split_A = true)
{
  // Mixed iterative/recursive version
  // Keep iterative while only one child do intersect
  // Did one recursive call if two children do_intersect and continue the current while
  const auto* current_A = &node_A;
  const auto* current_B = &node_B;
  while(traversal_traits.go_further())
  {
    if(split_A)
    {
      if(current_A->range() == 2)
      {
        traversal_traits.intersection(current_A->left_data(),  *current_B, current_B->range());
        traversal_traits.intersection(current_A->right_data(), *current_B, current_B->range());
        return;
      }
      else if(current_A->range() == 3)
      {
        // Treat directly the three primitives
        traversal_traits.intersection(current_A->left_data(), *current_B, current_B->range());
        traversal_traits.intersection(current_A->first_child().left_data() , *current_B, current_B->range());
        traversal_traits.intersection(current_A->first_child().right_data(), *current_B, current_B->range());
        return;
      }
      else
      {
        bool do_intersect_left  = traversal_traits.do_intersect(current_A->left_child(), *current_B);
        bool do_intersect_right = traversal_traits.do_intersect(current_A->right_child(), *current_B);
        if(do_intersect_left && do_intersect_right)
        {
          sequential_two_trees_traversal(current_A->left_child(), *current_B, traversal_traits,
                    traversal_traits.prefer_A_for_next_step(current_A->left_child(), *current_B, current_A->left_child().range(), current_B->range()));
          current_A = &current_A->right_child();
        }
        else if(!do_intersect_left && !do_intersect_right)
          return;
        else
          current_A = do_intersect_left ? &current_A->left_child() : &current_A->right_child();
      }
    }
    else
    {
      if(current_B->range() == 2)
      {
        traversal_traits.intersection(*current_A, current_A->range(), current_B->left_data());
        traversal_traits.intersection(*current_A, current_A->range(), current_B->right_data());
        return;
      }
      else if(current_B->range() == 3)
      {
        // Treat directly the three primitives
        traversal_traits.intersection(*current_A, current_A->range(), current_B->left_data());
        traversal_traits.intersection(*current_A, current_A->range(), current_B->first_child().left_data());
        traversal_traits.intersection(*current_A, current_A->range(), current_B->first_child().right_data());
        return;
      }
      else
      {
        bool do_intersect_left  = traversal_traits.do_intersect(*current_A, current_B->left_child());
        bool do_intersect_right = traversal_traits.do_intersect(*current_A, current_B->right_child());
        if(do_intersect_left && do_intersect_right)
        {
          sequential_two_trees_traversal(*current_A, current_B->left_child(), traversal_traits,
                    traversal_traits.prefer_A_for_next_step(*current_A, current_B->left_child(), current_A->range(), current_B->left_child().range()));
          current_B = &current_B->right_child();
        }
        else if(!do_intersect_left && !do_intersect_right)
          return;
        else
          current_B = do_intersect_left ? &current_B->left_child() : &current_B->right_child();
      }
    }
    split_A = traversal_traits.prefer_A_for_next_step(*current_A, *current_B, current_A->range(), current_B->range());
  }
}

#ifdef CGAL_LINKED_WITH_TBB

template<typename Node_A,
         typename Node_B,
         typename TwoTreeTraversalTraits>
class Two_trees_traversal_range
{
  const Node_A* m_node_A;
  const Node_B* m_node_B;

  // Cached splitting decision.
  bool m_empty = false;
  bool m_divisible = false;
  bool m_split_A;

  TwoTreeTraversalTraits& traits;

  void traversal_until_split_or_cutoff()
  {
    m_divisible = false;
    while(traits.go_further() &&
          m_node_A->range() >= grain_size &&
          m_node_B->range() >= grain_size)
    {
      if(m_split_A)
      {
        bool do_intersect_left  = traits.do_intersect(m_node_A->left_child(),  *m_node_B);
        bool do_intersect_right = traits.do_intersect(m_node_A->right_child(), *m_node_B);
        if(do_intersect_left && do_intersect_right){
          // Both intersect, trigger split
          m_divisible = true;
          return;
        }
        else if(!do_intersect_left && !do_intersect_right){
          // No intersection, set empty
          m_empty = true;
          return;
        }
        m_node_A = do_intersect_left ? &m_node_A->left_child() : &m_node_A->right_child();
      }
      else
      {
        bool do_intersect_left  = traits.do_intersect(*m_node_A, m_node_B->left_child());
        bool do_intersect_right = traits.do_intersect(*m_node_A, m_node_B->right_child());
        if(do_intersect_left && do_intersect_right){
          // Both intersect, trigger split
          m_divisible = true;
          return;
        }
        else if(!do_intersect_left && !do_intersect_right){
          // No intersection, set empty
          m_empty = true;
          return;
        }
        m_node_B = do_intersect_left ? &m_node_B->left_child() : &m_node_B->right_child();
      }
      m_split_A = traits.prefer_A_for_next_step(*m_node_A, *m_node_B, m_node_A->range(), m_node_B->range());
    }
  }

public:
  using size_type = std::size_t;

  // TODO finetune this value
  static constexpr size_type grain_size = 3000;

  Two_trees_traversal_range(
    const Node_A& node_A,
    const Node_B& node_B,
    TwoTreeTraversalTraits& traits)
    : m_node_A(&node_A),
      m_node_B(&node_B),
      m_split_A(true),
      traits(traits)
  {
    traversal_until_split_or_cutoff();
  }

  bool empty() const { return m_empty; }
  bool is_divisible() const { return m_divisible; }

  Two_trees_traversal_range(
    Two_trees_traversal_range& other,
    oneapi::tbb::split)
    : m_node_A(other.m_node_A),
      m_node_B(other.m_node_B),
      traits(other.traits)
  {
    CGAL_assertion(other.m_divisible);

    if(other.m_split_A)
    {
      // This range receives the right child and the original keeps the left child
      m_node_A = &other.m_node_A->right_child();
      other.m_node_A = &other.m_node_A->left_child();
    }
    else
    {
      // This range receives the right child and the original keeps the left child
      m_node_B = &other.m_node_B->right_child();
      other.m_node_B = &other.m_node_B->left_child();
    }
    m_split_A = traits.prefer_A_for_next_step(*m_node_A, *m_node_B, m_node_A->range(), m_node_B->range());
    other.m_split_A = traits.prefer_A_for_next_step(*other.m_node_A, *other.m_node_B, other.m_node_A->range(), other.m_node_B->range());

    // Do traversals until next splits
    other.traversal_until_split_or_cutoff();
    traversal_until_split_or_cutoff();
  }

  void run_in_sequential()
  {
    sequential_two_trees_traversal(*m_node_A, *m_node_B, traits, m_split_A);
  }
};
#endif

template<typename ConcurrencyTag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename TwoTreeTraversalTraits>
void two_trees_traversal(const Tree_A& tree_A,
                         const Tree_B& tree_B,
                         TwoTreeTraversalTraits &traits)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  if(tree_A.size() == 1){
    if(tree_B.size() == 1)
      traits.intersection(tree_A.singleton_data(), tree_B.singleton_data());
    else
      traits.intersection(tree_A.singleton_data(), *tree_B.root_node(), tree_B.size());
  }
  else if(tree_B.size() == 1)
  {
    traits.intersection(*tree_A.root_node(), tree_A.size(), tree_B.singleton_data());
  }
  else
  {
#if CGAL_LINKED_WITH_TBB
    using Range = Two_trees_traversal_range<typename Tree_A::Node,
                                            typename Tree_B::Node,
                                            TwoTreeTraversalTraits>;
    if constexpr(ConcurrencyTag::is_parallel)
    {
      Range range(*tree_A.root_node(), *tree_B.root_node(), traits);
      oneapi::tbb::parallel_for(range,
                                [](Range& r){ r.run_in_sequential(); },
                                oneapi::tbb::auto_partitioner{});
    }
    else
#endif
    {
      sequential_two_trees_traversal(*tree_A.root_node(), *tree_B.root_node(), traits, false);
    }
  }
}

template<typename ConcurrencyTag = Sequential_tag,
         typename Node_A,
         typename Node_B,
         typename TwoTreeTraversalTraits>
void two_trees_traversal(const Node_A& node_A,
                         const Node_B& node_B,
                         const std::size_t nb_primitives_A,
                         const std::size_t nb_primitives_B,
                         TwoTreeTraversalTraits &traits)
{
#if CGAL_LINKED_WITH_TBB
  using Range = Two_trees_traversal_range<Node_A, Node_B, TwoTreeTraversalTraits>;
  if(ConcurrencyTag::is_parallel && nb_primitives_A > Range::grain_size && nb_primitives_B > Range::grain_size)
  {
    Range range(node_A, node_B, traits);
    oneapi::tbb::parallel_for(range,
                              [](Range& r){ r.run_in_sequential(); },
                              oneapi::tbb::auto_partitioner{});
  }
  else
#endif
  {
    sequential_two_trees_traversal(node_A, node_B, traits);
  }
}

// An optimization to traverse a tree with itself
template <typename AABBTraits,
          typename TwoTreeTraversalTraits>
void sequential_one_tree_traversal(const ::CGAL::AABB_node<AABBTraits>& node,
                                   TwoTreeTraversalTraits& traversal_traits)
{
  // Mixed iterative/recursive version
  // Did recursive calls only for one child and if two children do intersect
  const auto* current = &node;
  while(traversal_traits.go_further())
  {
    if(current->range() == 2)
    {
      traversal_traits.intersection(current->left_data(), current->right_data());
      return;
    }
    else if(current->range() == 3)
    {
      // Treat directly the three primitives
      const auto& p1 = current->left_data(),
                  p2 = current->right_child().left_data(),
                  p3 = current->right_child().right_data();
      traversal_traits.intersection(p1, p2);
      traversal_traits.intersection(p2, p3);
      traversal_traits.intersection(p1, p3);
      return;
    }
    else
    {
      sequential_one_tree_traversal(current->left_child(), traversal_traits);
      if(traversal_traits.go_further() && traversal_traits.do_intersect(current->left_child(), current->right_child()))
        sequential_two_trees_traversal(current->left_child(), current->right_child(), traversal_traits);
      current = &current->right_child();
    }
  }
}

// An optimization to traverse a tree with itself
template <typename ConcurrencyTag = Sequential_tag,
          typename AABBTraits,
          typename TwoTreeTraversalTraits>
void one_tree_traversal(const ::CGAL::AABB_node<AABBTraits>& node,
                        TwoTreeTraversalTraits& traversal_traits)
{
  // TODO finetune and move in a functor
  const std::size_t cutoff = 3000;
#if CGAL_LINKED_WITH_TBB
  if constexpr(ConcurrencyTag::is_parallel)
  {
    if(node.range() > cutoff)
    {
      oneapi::tbb::task_group tg;
      if(traversal_traits.do_intersect(node.left_child(), node.right_child())){
        tg.run([&]{
          two_trees_traversal<ConcurrencyTag>(node.left_child(), node.right_child(), node.left_child().range(), node.right_child().range(), traversal_traits);}
        );
      }
      tg.run([&]{
              one_tree_traversal<ConcurrencyTag>(node.left_child(), traversal_traits);}
            );
      one_tree_traversal<ConcurrencyTag>(node.right_child(), traversal_traits);
      tg.wait();
    }
    else
    {
      sequential_one_tree_traversal(node, traversal_traits);
    }
  }
  else
#endif
  {
    sequential_one_tree_traversal(node, traversal_traits);
  }
}

template<typename ConcurrencyTag = Sequential_tag,
         typename Tree,
         typename TwoTreeTraversalTraits>
void one_tree_traversal(const Tree& tree,
                        TwoTreeTraversalTraits &traits)
{
  // Nothing to do if the tree is a singleton
  if(tree.size() < 2) return;
  one_tree_traversal<ConcurrencyTag>(*tree.root_node(), traits);
}

}}} // end of namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
