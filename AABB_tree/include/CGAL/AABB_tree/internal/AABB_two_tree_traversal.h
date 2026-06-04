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

#ifndef CGAL_AABB_TWO_TREE_TRAVERSAL_TRAITS_H
#define CGAL_AABB_TWO_TREE_TRAVERSAL_TRAITS_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_tree/internal/AABB_node.h>
#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <optional>

namespace CGAL {

namespace internal { namespace AABB_tree {

template< typename Value, typename OutputIterator >
class WrapOutputIterator
{
  Value first;
  OutputIterator out;
  bool in_order;
public:
  WrapOutputIterator(Value first_, OutputIterator out_, bool in_order_ = true): first(first_), out(out_), in_order(in_order_){}

  WrapOutputIterator& operator=(Value second){
    if (in_order)
      *out = std::make_pair(first, second);
    else
      *out = std::make_pair(second, first);
    return *this;
  }

  WrapOutputIterator& operator*(){ return *this; }
  WrapOutputIterator& operator++(){ ++out; return *this; }
  WrapOutputIterator& operator++(int){ ++out; return *this; }
  WrapOutputIterator& operator+(int d){ auto tmp = *this; ++out; return *this; }
};

template <bool in_order = true,
          typename Two_tree_traversal_traits,
          typename Concurrency_tag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename OutputIterator>
void two_tree_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                        const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                        const std::size_t nb_primitives_A,
                        const std::size_t nb_primitives_B,
                        const AABBTraits_A& traits_A,
                        const AABBTraits_B& traits_B,
                        OutputIterator out)
{
  using Primitive_id = typename AABBTraits_A::Primitive::Id;
  using WrapOutputIterator = WrapOutputIterator<Primitive_id, OutputIterator>;

  const std::size_t cutoff_parallel_call = 100000;
  switch(nb_primitives_A)
  {
  case 2:
  {
    WrapOutputIterator wrap_out_left(node_A.left_data().id(), out, in_order);
    Two_tree_traversal_traits traits_left(wrap_out_left, traits_B);
    node_B.traversal( internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.left_data(), traits_A), traits_left, nb_primitives_B);

    WrapOutputIterator wrap_out_right(node_A.right_data().id(), out, in_order);
    Two_tree_traversal_traits traits_right(wrap_out_right, traits_B);
    node_B.traversal( internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.right_data(), traits_A), traits_left, nb_primitives_B);
    break;
  }
  case 3:
  {
    WrapOutputIterator wrap_out_left(node_A.left_data().id(), out, in_order);
    Two_tree_traversal_traits traits_left(wrap_out_left, traits_B);
    node_B.traversal( internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.left_data(), traits_A), traits_left, nb_primitives_B);

    if( do_overlap(node_A.right_child().bbox(), node_B.bbox()) )
      two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.right_child(), nb_primitives_B, 2, traits_B, traits_A, out);
    break;
  }
  default:
  {
    // TODO current strategy is we swap at each step. Try largest tree first or largest bbox first
    bool do_intersect_left  = do_overlap(node_A.left_child().bbox(), node_B.bbox());
    bool do_intersect_right = do_overlap(node_A.right_child().bbox(), node_B.bbox());
#if CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                two_tree_traversal<!in_order, Two_tree_traversal_traits, Concurrency_tag>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traits_B, traits_A, out);}
              );
        two_tree_traversal<!in_order, Two_tree_traversal_traits, Concurrency_tag>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2,  traits_B, traits_A, out);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traits_B, traits_A, out);
        if( do_intersect_right )
          two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2,  traits_B, traits_A, out);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, traits_B, traits_A, out);
      if( do_intersect_right )
        two_tree_traversal<!in_order, Two_tree_traversal_traits>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2,  traits_B, traits_A, out);
    }
  }} // switch end
}

template <bool in_order = true,
          typename Two_tree_traversal_traits,
          typename Concurrency_tag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename OutputIterator>
void two_tree_traversal_iterative(const ::CGAL::AABB_node<AABBTraits_A>& root_A,
                                  const ::CGAL::AABB_node<AABBTraits_B>& root_B,
                                  std::size_t nb_primitives_A,
                                  std::size_t nb_primitives_B,
                                  const AABBTraits_A& traits_A,
                                  const AABBTraits_B& traits_B,
                                  OutputIterator out)
{
    using Primitive_id = typename AABBTraits_A::Primitive::Id;
    using WrapOutputIterator = WrapOutputIterator<Primitive_id, OutputIterator>;

    struct WorkItem {
        const ::CGAL::AABB_node<AABBTraits_A>* node_A;
        const ::CGAL::AABB_node<AABBTraits_B>* node_B;
        std::size_t nb_A;
        std::size_t nb_B;
        bool in_order_flag;
    };

    std::vector<WorkItem> stack;
    stack.push_back({&root_A, &root_B, nb_primitives_A, nb_primitives_B, in_order});

    const std::size_t cutoff_parallel_call = 100000;

    while (!stack.empty()) {
        WorkItem w = stack.back();
        stack.pop_back();

        const auto& node_A = *w.node_A;
        const auto& node_B = *w.node_B;

        switch (w.nb_A) {
            case 2: {
                WrapOutputIterator wrap_left(node_A.left_data().id(), out, w.in_order_flag);
                Two_tree_traversal_traits traits_left(wrap_left, traits_B);
                node_B.traversal(internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.left_data(), traits_A),
                                 traits_left, w.nb_B);

                WrapOutputIterator wrap_right(node_A.right_data().id(), out, w.in_order_flag);
                Two_tree_traversal_traits traits_right(wrap_right, traits_B);
                node_B.traversal(internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.right_data(), traits_A),
                                 traits_right, w.nb_B);
                break;
            }
            case 3: {
                WrapOutputIterator wrap_left(node_A.left_data().id(), out, w.in_order_flag);
                Two_tree_traversal_traits traits_left(wrap_left, traits_B);
                node_B.traversal(internal::Primitive_helper<AABBTraits_A>::get_datum(node_A.left_data(), traits_A),
                                 traits_left, w.nb_B);

                if (do_overlap(node_A.right_child().bbox(), node_B.bbox())) {
                    stack.push_back({&node_A.right_child(), &node_B, 2, w.nb_B, !w.in_order_flag});
                }
                break;
            }
            default: {
                bool do_left = do_overlap(node_A.left_child().bbox(), node_B.bbox());
                bool do_right = do_overlap(node_A.right_child().bbox(), node_B.bbox());

                if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>) {
                    if (do_left && do_right && w.nb_A > cutoff_parallel_call && w.nb_B > cutoff_parallel_call) {
                        oneapi::tbb::task_group tg;
                        tg.run([&] {
                            stack.push_back({&node_A.left_child(), &node_B, w.nb_A/2, w.nb_B, !w.in_order_flag});
                        });
                        stack.push_back({&node_A.right_child(), &node_B, w.nb_A - w.nb_A/2, w.nb_B, !w.in_order_flag});
                        tg.wait();
                        continue; // skip rest of loop iteration
                    }
                }

                if (do_right)
                    stack.push_back({&node_A.right_child(), &node_B, w.nb_A - w.nb_A/2, w.nb_B, !w.in_order_flag});
                if (do_left)
                    stack.push_back({&node_A.left_child(), &node_B, w.nb_A/2, w.nb_B, !w.in_order_flag});
            }
        }
    }
}

template<typename Concurrency_tag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename OutputIterator>
void two_tree_listing_intersecting_primitives(const Tree_A& tree_A,
                                              const Tree_B& tree_B,
                                              OutputIterator out)
{
  using Primitive_id = typename Tree_A::AABB_traits::Primitive::Id;
  using WrapOutputIterator = WrapOutputIterator<Primitive_id, OutputIterator>;
  using Two_tree_Listing_primitive_traits = CGAL::internal::AABB_tree::Listing_primitive_traits<typename Tree_B::AABB_traits, typename Tree_A::AABB_traits::Primitive::Datum, WrapOutputIterator>;
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_tree_traversal<true, Two_tree_Listing_primitive_traits, Concurrency_tag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), tree_A.traits(), tree_B.traits(), out);
}

template <bool in_order = true,
          typename Concurrency_tag = Sequential_tag,
          typename AABBTraits_A,
          typename AABBTraits_B,
          typename OutputIterator,
          typename Bbox_map_A,
          typename Bbox_map_B>
void two_tree_partial_traversal(const ::CGAL::AABB_node<AABBTraits_A>& node_A,
                                const ::CGAL::AABB_node<AABBTraits_B>& node_B,
                                const std::size_t nb_primitives_A,
                                const std::size_t nb_primitives_B,
                                const std::size_t cutoff,
                                const AABBTraits_A& traits_A,
                                const AABBTraits_B& traits_B,
                                const Bbox_map_A &bb_A,
                                const Bbox_map_B &bb_B,
                                OutputIterator out)
{
  const std::size_t cutoff_parallel_call = 100000;
  if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    if constexpr(in_order)
      *out++ = std::make_pair(std::addressof(node_A), std::addressof(node_B));
    else
      *out++ = std::make_pair(std::addressof(node_B), std::addressof(node_A));
  }
  else if(nb_primitives_A < cutoff && nb_primitives_B < cutoff)
  {
    two_tree_partial_traversal<!in_order>(node_B, node_A, nb_primitives_B, nb_primitives_A, cutoff, traits_B, traits_A, bb_B, bb_A, out);
  }
  else
  {
    // TODO current strategy is we swap at each step. Try largest tree first or largest bbox first
    bool do_intersect_left  = do_overlap(node_A.left_child().bbox(), node_B.bbox());
    bool do_intersect_right = do_overlap(node_A.right_child().bbox(), node_B.bbox());
#if CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      if(do_intersect_left && do_intersect_right && nb_primitives_A > cutoff_parallel_call && nb_primitives_B > cutoff_parallel_call)
      {
        oneapi::tbb::task_group tg;
        tg.run([&]{
                two_tree_partial_traversal<!in_order, Concurrency_tag>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traits_B, traits_A, bb_B, bb_A, out); }
              );
        two_tree_partial_traversal<!in_order, Concurrency_tag>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, cutoff,  traits_B, traits_A, bb_B, bb_A, out);
        tg.wait();
      }
      else
      {
        if( do_intersect_left )
          two_tree_partial_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traits_B, traits_A, bb_B, bb_A, out);
        if( do_intersect_right )
          two_tree_partial_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, cutoff,  traits_B, traits_A, bb_B, bb_A, out);
      }
    }
    else
#endif
    {
      if( do_intersect_left )
        two_tree_partial_traversal<!in_order>(node_B, node_A.left_child(), nb_primitives_B, nb_primitives_A/2, cutoff, traits_B, traits_A, bb_B, bb_A, out);
      if( do_intersect_right )
        two_tree_partial_traversal<!in_order>(node_B, node_A.right_child(), nb_primitives_B,  nb_primitives_A-nb_primitives_A/2, cutoff,  traits_B, traits_A, bb_B, bb_A, out);
    }
  }
}

template<typename Concurrency_tag = Sequential_tag,
         typename Tree_A,
         typename Tree_B,
         typename OutputIterator,
         typename Bbox_map_A,
         typename Bbox_map_B>
void two_tree_listing_intersecting_patches(const Tree_A& tree_A,
                                           const Tree_B& tree_B,
                                           OutputIterator out,
                                           const std::size_t cutoff,
                                           const Bbox_map_A &bb_A,
                                           const Bbox_map_B &bb_B)
{
  CGAL_precondition(tree_A.size() != 0 && tree_B.size() != 0);
  two_tree_partial_traversal<true, Concurrency_tag>(*tree_A.root_node(), *tree_B.root_node(), tree_A.size(), tree_B.size(), cutoff, tree_A.traits(), tree_B.traits(), bb_A, bb_B, out);
}


}}} // end namespace CGAL::internal::AABB_tree

#endif // CGAL_AABB_TRAVERSAL_TRAITS_H
